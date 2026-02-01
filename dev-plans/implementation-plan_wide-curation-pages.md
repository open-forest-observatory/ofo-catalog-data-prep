# Implementation Plan: Wide Layout for Curation Pages

## Goal
Make the drone mission curation pages wider to accommodate two side-by-side mission views, without affecting the width of other pages on the site.

## Overview
We will add a custom CSS class to identify curation pages, then add CSS rules that override the default container width only when that class is present.

---

## Step 1: Add a wrapper div with a custom class to the template

**File to edit:** `deploy/drone-imagery-ingestion/10_drone-mission-web-catalog/templates/drone-mission-details_curation.md`

**What to do:**
1. Open the template file
2. Immediately after the closing `---` of the front matter (line 7), add an opening `<div>` tag with a custom class
3. At the very end of the file, add the closing `</div>` tag

**Before (lines 1-9):**
```markdown
---
title: Mission {* dataset_id *} (Curation View)

date:
show_date: false
profile: false
---

<div class="container">
```

**After (lines 1-11):**
```markdown
---
title: Mission {* dataset_id *} (Curation View)

date:
show_date: false
profile: false
---

<div class="wide-page-content">

<div class="container">
```

4. At the very end of the file (after the closing `</script>` tag on line 151), add:
```html
</div>
```

---

## Step 2: Add CSS rules to widen the container on curation pages

**File to edit:** `ofo-website-3/assets/scss/custom.scss`

**What to do:**
1. Open the custom.scss file
2. Add the following CSS rules at the end of the file:

```scss
/* Wide layout for curation pages */
.wide-page-content {
  /* Override the parent container's max-width when inside wide-page-content */
}

.wide-page-content ~ .universal-wrapper,
.wide-page-content + .universal-wrapper,
.article-container:has(.wide-page-content),
body:has(.wide-page-content) .universal-wrapper {
  max-width: 95%;
}

/* Also widen any Bootstrap containers within the wide content */
.wide-page-content .container,
.wide-page-content .container-sm,
.wide-page-content .container-md,
.wide-page-content .container-lg,
.wide-page-content .container-xl {
  max-width: 100%;
}
```

---


## Step 4: If Step 2 doesn't work, try alternative CSS approach

Hugo Blox wraps page content in specific containers. If the CSS selectors in Step 2 don't work, try this alternative approach:

**Replace the CSS in Step 2 with:**
```scss
/* Wide layout for curation pages - Alternative approach */
/* This uses the :has() selector to detect pages with the wide-page-content class */
.universal-wrapper:has(.wide-page-content) {
  max-width: 95% !important;
  width: 95% !important;
}

/* Fallback for browsers that don't support :has() */
.wide-page-content {
  margin-left: -10%;
  margin-right: -10%;
  width: 120%;
}
```
