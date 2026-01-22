#!/usr/bin/env python3
"""
Revert Box files to their earliest versions.

This script reads a CSV file (in the format produced by box-file-versions.py)
and reverts each file to its earliest version using the Box "promote version" API.

The promote API creates a copy of the previous version and makes it current,
without needing to download and re-upload the file content.

Uses OAuth 2.0 authentication with a browser-based login flow. This works with
enterprise SSO accounts without requiring admin approval.

Setup:
    1. Go to https://developer.box.com and create a new app
    2. Choose "Custom App" with "User Authentication (OAuth 2.0)"
    3. Under "Configuration", add redirect URI: http://localhost:8080/callback
    4. Note your Client ID and Client Secret
    5. Either create a config file or pass credentials via CLI (see below)

Config file option (e.g., box_oauth_config.json):
    {
        "client_id": "YOUR_CLIENT_ID",
        "client_secret": "YOUR_CLIENT_SECRET"
    }

Usage (with config file):
    python box-revert-versions.py --input files_to_revert.csv --config box_oauth_config.json

Usage (with CLI credentials):
    python box-revert-versions.py --input files_to_revert.csv \\
        --client-id YOUR_CLIENT_ID --client-secret YOUR_CLIENT_SECRET

Dry run (preview without making changes):
    python box-revert-versions.py --input files_to_revert.csv --config box_oauth_config.json --dry-run

Headless operation:
    python box-revert-versions.py --input files_to_revert.csv \\
        --client-id YOUR_CLIENT_ID --client-secret YOUR_CLIENT_SECRET \\
        --access-token 'ACCESS_TOKEN' --refresh-token 'REFRESH_TOKEN'

CSV Format (as produced by box-file-versions.py):
    path,name,modified_at,version_count,size,id
    folder/file.jpg,file.jpg,2024-01-15 10:30:00,3,12345,123456789

Only the 'id' column is required; other columns are used for display purposes.
"""

import argparse
import csv
import json
import os
import sys
import threading
import time
import webbrowser
from datetime import datetime
from http.server import HTTPServer, BaseHTTPRequestHandler
from typing import Optional
from urllib.parse import urlparse, parse_qs

try:
    from box_sdk_gen import (
        BoxClient,
        BoxOAuth,
        OAuthConfig,
        GetAuthorizeUrlOptions,
        AccessToken,
        BoxAPIError,
        PromoteFileVersionType,
    )
except ImportError:
    print("Error: box-sdk-gen is not installed. Install it with: pip install box-sdk-gen")
    sys.exit(1)


TOKEN_CACHE_PATH = os.path.expanduser("~/.box_tokens.json")
REDIRECT_URI = "http://localhost:8080/callback"

# Rate limiting and retry settings
MAX_RETRIES = 10
BASE_BACKOFF_SECONDS = 15.0
MAX_BACKOFF_SECONDS = 300.0

# Global to capture auth code from callback
_auth_code = None
_auth_error = None

# Progress tracking
_files_processed = 0
_files_reverted = 0
_files_skipped = 0
_files_failed = 0
_files_already_reverted = 0
_last_progress_time = 0

# API call tracking
_api_call_timestamps = []
_start_time = None


def api_call_with_retry(func, *args, **kwargs):
    """
    Execute a Box API call with retry logic for rate limits and transient errors.

    Handles:
    - 429 Too Many Requests (rate limit)
    - 500, 502, 503, 504 (server errors)
    - Connection errors

    Uses exponential backoff with jitter.
    """
    global _api_call_timestamps
    last_exception = None

    for attempt in range(MAX_RETRIES):
        try:
            _api_call_timestamps.append(time.time())
            return func(*args, **kwargs)
        except BoxAPIError as e:
            last_exception = e

            # Get status code from the new SDK's error structure
            status = getattr(e, 'response_info', None)
            status_code = getattr(status, 'status_code', None) if status else None

            # Check if retryable
            if status_code in (429, 500, 502, 503, 504):
                # Exponential backoff with jitter
                wait_time = min(
                    BASE_BACKOFF_SECONDS * (2 ** attempt) + (time.time() % 1),
                    MAX_BACKOFF_SECONDS
                )

                print(f"Rate limited or server error (HTTP {status_code}). "
                      f"Retrying in {wait_time:.1f}s (attempt {attempt + 1}/{MAX_RETRIES})...",
                      file=sys.stderr)
                time.sleep(wait_time)
            else:
                # Non-retryable error
                raise
        except (ConnectionError, TimeoutError, OSError) as e:
            last_exception = e
            wait_time = min(
                BASE_BACKOFF_SECONDS * (2 ** attempt),
                MAX_BACKOFF_SECONDS
            )
            print(f"Connection error: {e}. "
                  f"Retrying in {wait_time:.1f}s (attempt {attempt + 1}/{MAX_RETRIES})...",
                  file=sys.stderr)
            time.sleep(wait_time)

    # Exhausted retries
    raise last_exception


def print_progress(force: bool = False):
    """Print progress update if enough time has passed."""
    global _last_progress_time

    now = time.time()
    # Print progress every 10 seconds, or if forced
    if force or (now - _last_progress_time) >= 10:
        # Calculate API call stats
        total_api_calls = len(_api_call_timestamps)
        one_minute_ago = now - 60
        calls_last_minute = sum(1 for t in _api_call_timestamps if t > one_minute_ago)

        # Calculate elapsed minutes
        elapsed_minutes = 0.0
        if _start_time:
            elapsed_minutes = (now - _start_time) / 60.0

        print(f"Progress: {_files_processed} processed, {_files_reverted} reverted, "
              f"{_files_already_reverted} already reverted, {_files_skipped} skipped, {_files_failed} failed | "
              f"API calls: {total_api_calls} total, {calls_last_minute} in last min, "
              f"{elapsed_minutes:.1f} min elapsed",
              file=sys.stderr)
        _last_progress_time = now


class OAuthCallbackHandler(BaseHTTPRequestHandler):
    """HTTP handler to capture OAuth callback."""

    def do_GET(self):
        global _auth_code, _auth_error

        parsed = urlparse(self.path)
        params = parse_qs(parsed.query)

        if 'code' in params:
            _auth_code = params['code'][0]
            self.send_response(200)
            self.send_header('Content-type', 'text/html')
            self.end_headers()
            self.wfile.write(b"""
                <html><body>
                <h1>Authorization successful!</h1>
                <p>You can close this window and return to the terminal.</p>
                </body></html>
            """)
        elif 'error' in params:
            _auth_error = params.get('error_description', params['error'])[0]
            self.send_response(400)
            self.send_header('Content-type', 'text/html')
            self.end_headers()
            self.wfile.write(f"""
                <html><body>
                <h1>Authorization failed</h1>
                <p>{_auth_error}</p>
                </body></html>
            """.encode())
        else:
            self.send_response(400)
            self.end_headers()

    def log_message(self, format, *args):
        pass  # Suppress HTTP server logs


def save_tokens(access_token: str, refresh_token: str):
    """Save tokens to cache file."""
    with open(TOKEN_CACHE_PATH, 'w') as f:
        json.dump({
            'access_token': access_token,
            'refresh_token': refresh_token
        }, f)
    os.chmod(TOKEN_CACHE_PATH, 0o600)


def load_tokens() -> Optional[dict]:
    """Load tokens from cache file."""
    if os.path.exists(TOKEN_CACHE_PATH):
        try:
            with open(TOKEN_CACHE_PATH, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return None


def authenticate(
    config_path: Optional[str] = None,
    client_id: Optional[str] = None,
    client_secret: Optional[str] = None,
    access_token: Optional[str] = None,
    refresh_token: Optional[str] = None
) -> BoxClient:
    """
    Authenticate with Box using OAuth 2.0.

    If access_token/refresh_token are provided, uses those directly (no storage).
    Otherwise, checks cached tokens, then falls back to browser auth flow.

    Args:
        config_path: Path to JSON file with client_id and client_secret
        client_id: Box OAuth client ID (alternative to config_path)
        client_secret: Box OAuth client secret (alternative to config_path)
        access_token: Optional access token to use directly
        refresh_token: Optional refresh token to use directly

    Returns:
        Authenticated Box client
    """
    global _auth_code, _auth_error

    # Load client credentials from config file or use provided values
    if config_path:
        with open(config_path, 'r') as f:
            config = json.load(f)
        client_id = config['client_id']
        client_secret = config['client_secret']
    elif not client_id or not client_secret:
        raise ValueError("Must provide either --config or both --client-id and --client-secret")

    # Create OAuth config
    oauth_config = OAuthConfig(client_id=client_id, client_secret=client_secret)
    auth = BoxOAuth(oauth_config)

    # If tokens provided via CLI, use them directly
    if access_token:
        token = AccessToken(access_token=access_token, refresh_token=refresh_token)
        auth.token_storage.store(token)
        client = BoxClient(auth=auth)
        # Test the connection
        client.users.get_user_me()
        print("Using provided tokens.", file=sys.stderr)
        return client

    # Check for cached tokens
    cached = load_tokens()
    if cached:
        try:
            token = AccessToken(
                access_token=cached.get('access_token'),
                refresh_token=cached.get('refresh_token')
            )
            auth.token_storage.store(token)
            client = BoxClient(auth=auth)
            # Test the connection
            client.users.get_user_me()
            print("Using cached authentication.", file=sys.stderr)
            return client
        except Exception:
            print("Cached tokens expired, re-authenticating...", file=sys.stderr)
            # Create fresh auth object
            auth = BoxOAuth(oauth_config)

    # Need to do fresh OAuth flow
    auth_url = auth.get_authorize_url(
        options=GetAuthorizeUrlOptions(redirect_uri=REDIRECT_URI)
    )

    # Start local server to receive callback
    server = HTTPServer(('localhost', 8080), OAuthCallbackHandler)
    server_thread = threading.Thread(target=server.handle_request)
    server_thread.start()

    # Open browser for user authorization
    print("Opening browser for Box authorization...", file=sys.stderr)
    print(f"If browser doesn't open, visit: {auth_url}", file=sys.stderr)
    webbrowser.open(auth_url)

    # Wait for callback
    server_thread.join(timeout=120)
    server.server_close()

    if _auth_error:
        raise Exception(f"Authorization failed: {_auth_error}")
    if not _auth_code:
        raise Exception("Authorization timed out. Please try again.")

    # Exchange code for tokens
    auth.get_tokens_authorization_code_grant(_auth_code)

    # Save tokens from storage
    token = auth.token_storage.get()
    if token:
        save_tokens(token.access_token, token.refresh_token)

    return BoxClient(auth=auth)


def get_version_count(client: BoxClient, file_id: str) -> int:
    """
    Get the total number of versions for a file.

    Returns the count including the current version (so a file with 1 previous
    version returns 2).
    """
    try:
        def fetch_versions():
            return client.file_versions.get_file_versions(file_id)

        versions_response = api_call_with_retry(fetch_versions)

        # The API returns previous versions (not including current),
        # so add 1 for the current version
        return len(versions_response.entries) + 1 if versions_response.entries else 1
    except BoxAPIError as e:
        status = getattr(e, 'response_info', None)
        status_code = getattr(status, 'status_code', None) if status else None
        if status_code == 404:
            return 0
        raise


def get_earliest_version(client: BoxClient, file_id: str) -> Optional[dict]:
    """
    Get the earliest version of a file.

    Returns the oldest version (not the current one).
    Returns None if there is no previous version.
    """
    try:
        def fetch_versions():
            return client.file_versions.get_file_versions(file_id)

        versions_response = api_call_with_retry(fetch_versions)

        # Debug: print what we're getting
        print(f"DEBUG: versions_response type={type(versions_response)}, entries={versions_response.entries}", file=sys.stderr)

        if not versions_response.entries:
            return None

        # Find the earliest version by created_at timestamp
        earliest_version = min(
            versions_response.entries,
            key=lambda v: getattr(v, 'created_at', '') or ''
        )
        return {
            'id': earliest_version.id,
            'sha1': getattr(earliest_version, 'sha1', None),
            'created_at': getattr(earliest_version, 'created_at', None),
            'modified_at': getattr(earliest_version, 'modified_at', None),
            'size': getattr(earliest_version, 'size', None),
        }
    except BoxAPIError as e:
        status = getattr(e, 'response_info', None)
        status_code = getattr(status, 'status_code', None) if status else None
        if status_code == 404:
            return None
        raise


def promote_version(client: BoxClient, file_id: str, version_id: str) -> dict:
    """
    Promote a previous version to become the current version.

    This creates a copy of the old version and makes it the current version.
    No download/re-upload is required.

    Args:
        client: Authenticated Box client
        file_id: The ID of the file
        version_id: The ID of the version to promote

    Returns:
        Dict with the new version info
    """
    def do_promote():
        return client.file_versions.promote_file_version(
            file_id,
            id=version_id,
            type=PromoteFileVersionType.FILE_VERSION
        )

    result = api_call_with_retry(do_promote)
    return {
        'id': result.id,
        'sha1': getattr(result, 'sha1', None),
        'version_number': getattr(result, 'version_number', None),
    }


def format_size(size_bytes) -> str:
    """Format file size in human-readable form."""
    if size_bytes is None:
        return "N/A"

    try:
        size_bytes = int(size_bytes)
    except (ValueError, TypeError):
        return "N/A"

    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} PB"


def read_input_csv(input_path: str) -> list[dict]:
    """
    Read the input CSV file and return list of file records.

    Args:
        input_path: Path to the CSV file

    Returns:
        List of dicts with file info (at minimum 'id' key)
    """
    files = []
    with open(input_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)

        # Verify required column exists
        if 'id' not in reader.fieldnames:
            raise ValueError("CSV file must have an 'id' column")

        for row in reader:
            files.append(row)

    return files


def main():
    parser = argparse.ArgumentParser(
        description="Revert Box files to their previous versions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Required arguments
    parser.add_argument(
        '--input', '-i',
        required=True,
        help="Path to CSV file with files to revert (must have 'id' column)"
    )

    # Authentication options (config file OR client-id/client-secret)
    parser.add_argument(
        '--config',
        help="Path to Box OAuth config JSON file (with client_id and client_secret)"
    )
    parser.add_argument(
        '--client-id',
        help="Box OAuth client ID (alternative to --config)"
    )
    parser.add_argument(
        '--client-secret',
        help="Box OAuth client secret (alternative to --config)"
    )

    # Operation options
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help="Preview what would be reverted without making changes"
    )
    parser.add_argument(
        '--output', '-o',
        help="Output file path for results log (CSV format)"
    )

    # Token options for headless operation
    parser.add_argument(
        '--access-token',
        help="Box access token (for headless operation, bypasses stored tokens)"
    )
    parser.add_argument(
        '--refresh-token',
        help="Box refresh token (for headless operation, use with --access-token)"
    )
    parser.add_argument(
        '--clear-tokens',
        action='store_true',
        help="Clear cached tokens and re-authenticate"
    )
    parser.add_argument(
        '--print-tokens',
        action='store_true',
        help="Print current tokens after auth (for copying to headless machine)"
    )

    args = parser.parse_args()

    # Handle token clearing
    if args.clear_tokens and os.path.exists(TOKEN_CACHE_PATH):
        os.remove(TOKEN_CACHE_PATH)
        print("Cleared cached tokens.", file=sys.stderr)

    # Validate auth options
    if not args.config and not (args.client_id and args.client_secret):
        print("Error: Must provide either --config or both --client-id and --client-secret",
              file=sys.stderr)
        sys.exit(1)

    # Read input file
    print(f"Reading input file: {args.input}", file=sys.stderr)
    try:
        files_to_revert = read_input_csv(args.input)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files_to_revert)} files to process", file=sys.stderr)

    if not files_to_revert:
        print("No files to process. Exiting.", file=sys.stderr)
        sys.exit(0)

    # Authenticate
    print("Authenticating with Box...", file=sys.stderr)
    try:
        client = authenticate(
            config_path=args.config,
            client_id=args.client_id,
            client_secret=args.client_secret,
            access_token=args.access_token,
            refresh_token=args.refresh_token
        )
    except Exception as e:
        print(f"Error: Failed to authenticate: {e}", file=sys.stderr)
        sys.exit(1)

    # Print tokens if requested (for headless setup)
    if args.print_tokens:
        tokens = load_tokens()
        if tokens:
            print("\n=== Tokens (copy these for headless use) ===", file=sys.stderr)
            print(f"--access-token '{tokens.get('access_token')}'", file=sys.stderr)
            print(f"--refresh-token '{tokens.get('refresh_token')}'", file=sys.stderr)
            print("", file=sys.stderr)

    if args.dry_run:
        print("\n*** DRY RUN MODE - No changes will be made ***\n", file=sys.stderr)

    # Setup output file if specified
    output_file = None
    csv_writer = None
    if args.output:
        output_file = open(args.output, 'w', newline='', encoding='utf-8')
        csv_writer = csv.DictWriter(output_file, fieldnames=[
            'file_id', 'path', 'name', 'status', 'previous_version_id',
            'new_version_id', 'error'
        ])
        csv_writer.writeheader()

    # Process files
    global _files_processed, _files_reverted, _files_skipped, _files_failed, _files_already_reverted, _start_time
    _start_time = time.time()

    results = []

    print("\nProcessing files...", file=sys.stderr)
    print("-" * 80, file=sys.stderr)

    try:
        for file_info in files_to_revert:
            _files_processed += 1
            file_id = file_info['id']
            file_path = file_info.get('path', '')
            file_name = file_info.get('name', f'file_{file_id}')

            result = {
                'file_id': file_id,
                'path': file_path,
                'name': file_name,
                'status': '',
                'previous_version_id': '',
                'new_version_id': '',
                'error': ''
            }

            try:
                # Check version count first - if > 2, file was already reverted
                version_count = get_version_count(client, file_id)

                if version_count == 3:
                    result['status'] = 'already_reverted'
                    result['error'] = 'File has 3 versions (already reverted)'
                    _files_already_reverted += 1
                    print_progress()
                    results.append(result)
                    if csv_writer:
                        csv_writer.writerow(result)
                        output_file.flush()
                    continue

                # Get the earliest version
                earliest_version = get_earliest_version(client, file_id)

                if not earliest_version:
                    result['status'] = 'skipped'
                    result['error'] = 'No previous version available'
                    _files_skipped += 1
                    print(f"SKIP: {file_path or file_name} - No previous version", file=sys.stderr)
                else:
                    result['previous_version_id'] = earliest_version['id']

                    if args.dry_run:
                        result['status'] = 'dry_run'
                        print(f"WOULD REVERT: {file_path or file_name} -> version {earliest_version['id']}", file=sys.stderr)
                        _files_reverted += 1
                    else:
                        # Promote the earliest version
                        new_version = promote_version(client, file_id, earliest_version['id'])
                        result['status'] = 'reverted'
                        result['new_version_id'] = new_version['id']
                        _files_reverted += 1
                        print(f"REVERTED: {file_path or file_name} -> version {earliest_version['id']} (new: {new_version['id']})", file=sys.stderr)

            except BoxAPIError as e:
                status = getattr(e, 'response_info', None)
                status_code = getattr(status, 'status_code', None) if status else None
                error_msg = f"API error (HTTP {status_code}): {e}"
                result['status'] = 'failed'
                result['error'] = error_msg
                _files_failed += 1
                print(f"FAILED: {file_path or file_name} - {error_msg}", file=sys.stderr)

            except Exception as e:
                result['status'] = 'failed'
                result['error'] = str(e)
                _files_failed += 1
                print(f"FAILED: {file_path or file_name} - {e}", file=sys.stderr)

            results.append(result)

            # Write to output file immediately
            if csv_writer:
                csv_writer.writerow(result)
                output_file.flush()

            print_progress()

    finally:
        if output_file:
            output_file.close()

    # Final progress and summary
    print("-" * 80, file=sys.stderr)
    print_progress(force=True)

    print("\n=== Summary ===", file=sys.stderr)
    print(f"Total files processed: {_files_processed}", file=sys.stderr)
    print(f"Successfully reverted: {_files_reverted}", file=sys.stderr)
    print(f"Already reverted (3 versions): {_files_already_reverted}", file=sys.stderr)
    print(f"Skipped (no previous version): {_files_skipped}", file=sys.stderr)
    print(f"Failed: {_files_failed}", file=sys.stderr)

    if args.output:
        print(f"\nResults written to: {args.output}", file=sys.stderr)

    if args.dry_run:
        print("\n*** This was a DRY RUN - no changes were made ***", file=sys.stderr)

    # Exit with error code if any failures
    if _files_failed > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
