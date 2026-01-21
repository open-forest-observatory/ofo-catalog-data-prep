#!/usr/bin/env python3
"""
List Box files matching date range and/or version criteria.

This script recursively searches a Box folder and lists files based on:
- Modification date range (optional)
- File version number (optional)

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
    python box-file-versions.py --folder-id FOLDER_ID --config box_oauth_config.json \\
        --start-date 2024-01-01 --end-date 2024-12-31

Usage (with CLI credentials):
    python box-file-versions.py --folder-id FOLDER_ID \\
        --client-id YOUR_CLIENT_ID --client-secret YOUR_CLIENT_SECRET \\
        --start-date 2024-01-01 --end-date 2024-12-31

Headless operation:
    1. On a machine with a browser, authenticate and get tokens:
       python box-file-versions.py --folder-id 0 --config box_oauth_config.json --print-tokens

    2. Copy the printed tokens to your headless machine and run:
       python box-file-versions.py --folder-id FOLDER_ID \\
           --client-id YOUR_CLIENT_ID --client-secret YOUR_CLIENT_SECRET \\
           --access-token 'ACCESS_TOKEN' --refresh-token 'REFRESH_TOKEN' \\
           --start-date 2024-01-01 --output results.csv

On first run (without --access-token), your browser will open for you to log in
via SSO. The token is cached in ~/.box_tokens.json for subsequent runs.
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
from typing import Generator, Optional
from urllib.parse import urlparse, parse_qs

try:
    from boxsdk import OAuth2, Client
    from boxsdk.object.folder import Folder
    from boxsdk.exception import BoxAPIException
except ImportError:
    print("Error: boxsdk is not installed. Install it with: pip install boxsdk")
    sys.exit(1)

try:
    from dateutil import parser as date_parser
except ImportError:
    print("Error: python-dateutil is not installed. Install it with: pip install python-dateutil")
    sys.exit(1)


TOKEN_CACHE_PATH = os.path.expanduser("~/.box_tokens.json")
REDIRECT_URI = "http://localhost:8080/callback"

# Rate limiting and retry settings
MAX_RETRIES = 5
BASE_BACKOFF_SECONDS = 1.0
MAX_BACKOFF_SECONDS = 60.0

# Global to capture auth code from callback
_auth_code = None
_auth_error = None

# Progress tracking
_files_scanned = 0
_folders_scanned = 0
_files_matched = 0
_last_progress_time = 0


def api_call_with_retry(func, *args, **kwargs):
    """
    Execute a Box API call with retry logic for rate limits and transient errors.

    Handles:
    - 429 Too Many Requests (rate limit)
    - 500, 502, 503, 504 (server errors)
    - Connection errors

    Uses exponential backoff with jitter.
    """
    last_exception = None

    for attempt in range(MAX_RETRIES):
        try:
            return func(*args, **kwargs)
        except BoxAPIException as e:
            last_exception = e

            # Check if retryable
            if e.status in (429, 500, 502, 503, 504):
                # Get retry-after header if available (for 429)
                retry_after = None
                if e.status == 429 and hasattr(e, 'headers'):
                    retry_after = e.headers.get('Retry-After')

                if retry_after:
                    wait_time = float(retry_after)
                else:
                    # Exponential backoff with jitter
                    wait_time = min(
                        BASE_BACKOFF_SECONDS * (2 ** attempt) + (time.time() % 1),
                        MAX_BACKOFF_SECONDS
                    )

                print(f"Rate limited or server error (HTTP {e.status}). "
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
        print(f"Progress: {_folders_scanned} folders, {_files_scanned} files scanned, "
              f"{_files_matched} files matched so far...",
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
) -> Client:
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

    # If tokens provided via CLI, use them directly (no storage)
    if access_token:
        oauth = OAuth2(
            client_id=client_id,
            client_secret=client_secret,
            access_token=access_token,
            refresh_token=refresh_token,
            store_tokens=lambda a, r: None  # Don't store CLI-provided tokens
        )
        client = Client(oauth)
        # Test the connection
        client.user().get()
        print("Using provided tokens.", file=sys.stderr)
        return client

    # Token storage callback for SDK
    def store_tokens(access_token, refresh_token):
        save_tokens(access_token, refresh_token)

    # Check for cached tokens
    cached = load_tokens()
    if cached:
        try:
            oauth = OAuth2(
                client_id=client_id,
                client_secret=client_secret,
                access_token=cached.get('access_token'),
                refresh_token=cached.get('refresh_token'),
                store_tokens=store_tokens
            )
            client = Client(oauth)
            # Test the connection
            client.user().get()
            print("Using cached authentication.", file=sys.stderr)
            return client
        except Exception:
            print("Cached tokens expired, re-authenticating...", file=sys.stderr)

    # Need to do fresh OAuth flow
    oauth = OAuth2(
        client_id=client_id,
        client_secret=client_secret,
        store_tokens=store_tokens
    )

    auth_url, csrf_token = oauth.get_authorization_url(REDIRECT_URI)

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
    access_token, refresh_token = oauth.authenticate(_auth_code)
    save_tokens(access_token, refresh_token)

    return Client(oauth)


def get_file_version_count(client: Client, file_id: str) -> int:
    """
    Get the number of versions for a file.

    Box counts the current version plus all previous versions.
    A file with no previous versions has version count of 1.
    """
    try:
        file_obj = client.file(file_id)

        def fetch_versions():
            return list(file_obj.get_previous_versions())

        versions = api_call_with_retry(fetch_versions)
        return len(versions) + 1  # +1 for current version
    except Exception as e:
        print(f"Warning: Could not get versions for file {file_id}: {e}", file=sys.stderr)
        return 1


def list_files_recursively(
    client: Client,
    folder: Folder,
    path_prefix: str = ""
) -> Generator[dict, None, None]:
    """
    Recursively list all files in a folder and its subfolders.

    Args:
        client: Authenticated Box client
        folder: Box folder object to start from
        path_prefix: Current path for building relative paths

    Yields:
        Dict with file information (id, name, path, modified_at, size)
    """
    global _files_scanned, _folders_scanned

    try:
        def fetch_items():
            return list(folder.get_items(
                limit=1000,
                fields=['id', 'name', 'type', 'modified_at', 'size', 'content_modified_at']
            ))

        items = api_call_with_retry(fetch_items)
        _folders_scanned += 1

        for item in items:
            current_path = f"{path_prefix}/{item.name}" if path_prefix else item.name

            if item.type == 'file':
                _files_scanned += 1
                print_progress()
                yield {
                    'id': item.id,
                    'name': item.name,
                    'path': current_path,
                    'size': getattr(item, 'size', 0),
                    'modified_at': getattr(item, 'modified_at', None),
                    'content_modified_at': getattr(item, 'content_modified_at', None),
                }
            elif item.type == 'folder':
                subfolder = client.folder(item.id)
                yield from list_files_recursively(client, subfolder, current_path)

    except Exception as e:
        print(f"Error listing folder '{path_prefix}': {e}", file=sys.stderr)


def parse_date(date_str: str) -> datetime:
    """Parse a date string in various formats."""
    try:
        # Try ISO format first
        return datetime.fromisoformat(date_str.replace('Z', '+00:00'))
    except ValueError:
        pass

    # Try common formats
    for fmt in ['%Y-%m-%d', '%Y/%m/%d', '%m/%d/%Y', '%d-%m-%Y']:
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue

    # Fall back to dateutil parser
    return date_parser.parse(date_str)


def filter_files(
    client: Client,
    files: Generator[dict, None, None],
    start_date: Optional[datetime] = None,
    end_date: Optional[datetime] = None,
    min_version: Optional[int] = None,
    max_version: Optional[int] = None,
) -> Generator[dict, None, None]:
    """
    Filter files by date range and/or version count.

    Args:
        client: Authenticated Box client (needed for version queries)
        files: Generator of file dicts from list_files_recursively
        start_date: Include files modified on or after this date
        end_date: Include files modified on or before this date
        min_version: Include files with at least this many versions
        max_version: Include files with at most this many versions

    Yields:
        Filtered file dicts with version_count added
    """
    global _files_matched

    need_versions = min_version is not None or max_version is not None

    for file_info in files:
        # Check date filter
        if start_date or end_date:
            # Use content_modified_at if available, fall back to modified_at
            mod_date_str = file_info.get('content_modified_at') or file_info.get('modified_at')
            if mod_date_str:
                try:
                    mod_date = parse_date(mod_date_str)
                    # Make naive datetime for comparison if needed
                    if mod_date.tzinfo is not None and start_date and start_date.tzinfo is None:
                        mod_date = mod_date.replace(tzinfo=None)

                    if start_date and mod_date < start_date:
                        continue
                    if end_date and mod_date > end_date:
                        continue
                except Exception as e:
                    print(f"Warning: Could not parse date '{mod_date_str}': {e}", file=sys.stderr)

        # Check version filter
        if need_versions:
            version_count = get_file_version_count(client, file_info['id'])
            file_info['version_count'] = version_count

            if min_version and version_count < min_version:
                continue
            if max_version and version_count > max_version:
                continue
        else:
            file_info['version_count'] = None

        _files_matched += 1
        yield file_info


def format_size(size_bytes: int) -> str:
    """Format file size in human-readable form."""
    if size_bytes is None:
        return "N/A"

    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} PB"


def main():
    parser = argparse.ArgumentParser(
        description="List Box files by modification date and/or version number",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Required arguments
    parser.add_argument(
        '--folder-id',
        required=True,
        help="Box folder ID to search (use '0' for root folder)"
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

    # Optional date filters
    parser.add_argument(
        '--start-date',
        help="Include files modified on or after this date (YYYY-MM-DD)"
    )
    parser.add_argument(
        '--end-date',
        help="Include files modified on or before this date (YYYY-MM-DD)"
    )

    # Optional version filters
    parser.add_argument(
        '--min-version',
        type=int,
        help="Include files with at least this many versions"
    )
    parser.add_argument(
        '--max-version',
        type=int,
        help="Include files with at most this many versions"
    )
    parser.add_argument(
        '--version',
        type=int,
        help="Include files with exactly this version count (shorthand for --min-version N --max-version N)"
    )

    # Output options
    parser.add_argument(
        '--output', '-o',
        help="Output file path (CSV format). If not specified, prints to stdout"
    )
    parser.add_argument(
        '--clear-tokens',
        action='store_true',
        help="Clear cached tokens and re-authenticate"
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
        '--print-tokens',
        action='store_true',
        help="Print current tokens after auth (for copying to headless machine)"
    )

    args = parser.parse_args()

    # Handle token clearing
    if args.clear_tokens and os.path.exists(TOKEN_CACHE_PATH):
        os.remove(TOKEN_CACHE_PATH)
        print("Cleared cached tokens.", file=sys.stderr)

    # Handle --version shorthand
    min_version = args.min_version
    max_version = args.max_version
    if args.version is not None:
        min_version = args.version
        max_version = args.version

    # Parse dates
    start_date = None
    end_date = None
    if args.start_date:
        start_date = parse_date(args.start_date)
    if args.end_date:
        end_date = parse_date(args.end_date)
        # If end_date has no time, set to end of day
        if end_date.hour == 0 and end_date.minute == 0 and end_date.second == 0:
            end_date = end_date.replace(hour=23, minute=59, second=59)

    # Validate that at least one filter is specified
    if not any([start_date, end_date, min_version, max_version]):
        print("Warning: No filters specified. Listing all files.", file=sys.stderr)

    # Validate auth options
    if not args.config and not (args.client_id and args.client_secret):
        print("Error: Must provide either --config or both --client-id and --client-secret",
              file=sys.stderr)
        sys.exit(1)

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

    # Get the folder
    print(f"Searching folder {args.folder_id}...", file=sys.stderr)
    try:
        folder = client.folder(args.folder_id)
        folder_info = folder.get(fields=['name'])
        folder_name = folder_info.name
    except Exception as e:
        print(f"Error: Could not access folder {args.folder_id}: {e}", file=sys.stderr)
        sys.exit(1)

    # Build filter description for output
    filters = []
    if start_date:
        filters.append(f"from {start_date.strftime('%Y-%m-%d')}")
    if end_date:
        filters.append(f"to {end_date.strftime('%Y-%m-%d')}")
    if min_version:
        filters.append(f"min {min_version} versions")
    if max_version:
        filters.append(f"max {max_version} versions")
    filter_desc = ", ".join(filters) if filters else "no filters"

    print(f"Folder: {folder_name}", file=sys.stderr)
    print(f"Filters: {filter_desc}", file=sys.stderr)
    print("", file=sys.stderr)

    # List and filter files
    print("Scanning folders...", file=sys.stderr)
    files = list_files_recursively(client, folder)
    filtered_files = filter_files(
        client, files,
        start_date=start_date,
        end_date=end_date,
        min_version=min_version,
        max_version=max_version
    )

    # Collect results
    results = list(filtered_files)

    # Final progress
    print_progress(force=True)
    print(f"Scan complete.", file=sys.stderr)

    # Output results
    if args.output:
        # Write to CSV
        with open(args.output, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'path', 'name', 'modified_at', 'version_count', 'size', 'id'
            ])
            writer.writeheader()
            for file_info in results:
                writer.writerow({
                    'path': file_info['path'],
                    'name': file_info['name'],
                    'modified_at': file_info.get('content_modified_at') or file_info.get('modified_at', ''),
                    'version_count': file_info.get('version_count', ''),
                    'size': file_info.get('size', ''),
                    'id': file_info['id']
                })
        print(f"Wrote {len(results)} files to {args.output}", file=sys.stderr)
    else:
        # Print human-readable table
        if not results:
            print("No files found matching the criteria.")
        else:
            # Calculate column widths
            path_width = max(len(f['path']) for f in results)
            path_width = min(path_width, 80)  # Cap at 80 chars

            # Print header
            header = f"{'Path':<{path_width}}  {'Modified':<20}  {'Ver':>4}  {'Size':>10}"
            print(header)
            print("-" * len(header))

            # Print files
            for file_info in results:
                path = file_info['path']
                if len(path) > path_width:
                    path = "..." + path[-(path_width-3):]

                mod_date = file_info.get('content_modified_at') or file_info.get('modified_at', '')
                if mod_date:
                    try:
                        dt = parse_date(mod_date)
                        mod_date = dt.strftime('%Y-%m-%d %H:%M:%S')
                    except:
                        pass

                version = file_info.get('version_count')
                version_str = str(version) if version is not None else '-'

                size_str = format_size(file_info.get('size'))

                print(f"{path:<{path_width}}  {mod_date:<20}  {version_str:>4}  {size_str:>10}")

            print("")
            print(f"Total: {len(results)} files")


if __name__ == '__main__':
    main()
