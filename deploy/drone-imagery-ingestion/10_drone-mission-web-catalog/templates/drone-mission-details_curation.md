---
title: Mission {* dataset_id *} (Curation View)

date:
show_date: false
profile: false
---

<div class="container">
    <div class="row">
        <div class="col-sm">
            <div class="text-center">
            {{< cta cta_text="Previous mission" cta_link="{* previous_dataset_page_path *}" cta_new_tab="false" >}}
            </div>
        </div>
        <div class="col-sm">
            <div class="text-center">
            {{< cta cta_text="Next mission" cta_link="{* next_dataset_page_path *}" cta_new_tab="false" >}}
            </div>
        </div>
    </div>
</div>

<p style="line-height: 125%; text-align:center;"><b>Internal curation view.</b> This page shows mission-level metadata and image locations for data curation purposes.</p>

<!-- View Toggle UI -->
<div style="text-align: center; margin: 20px 0;">
    <div class="btn-group" role="group" aria-label="View toggle">
        <button type="button" class="btn btn-outline-primary view-toggle-btn" onclick="setView('pre')">Pre-curation</button>
        <button type="button" class="btn btn-outline-primary view-toggle-btn" onclick="setView('post')">Post-curation</button>
        <button type="button" class="btn btn-primary view-toggle-btn active" onclick="setView('both')">Both</button>
    </div>
</div>

<style>
.no-data-message {
    color: #6c757d;
    font-style: italic;
    text-align: center;
    padding: 50px 20px;
    background-color: #f8f9fa;
    border: 1px dashed #dee2e6;
    border-radius: 4px;
    min-height: 200px;
    display: flex;
    align-items: center;
    justify-content: center;
}
.pre-curation-panel, .post-curation-panel {
    padding: 10px;
}
.panel-header {
    font-size: 1.1em;
    font-weight: 600;
    margin-bottom: 10px;
    text-align: center;
}
.view-toggle-btn.active {
    font-weight: bold;
}
</style>

## Image locations and attributes

<div class="row">
<div class="col-sm-6 pre-curation-panel">
<div class="panel-header">Pre-curation</div>
<iframe src="{* map_html_path *}" frameborder="0" scrolling="yes" seamless="seamless" style="display:block; width:100%; height:65vh; background: rgba(0,0,0,0);"></iframe>
</div>
<div class="col-sm-6 post-curation-panel">
<div class="panel-header">Post-curation</div>
{% if has_post_curation_data -%}
<iframe src="{* post_curation_map_html_path *}" frameborder="0" scrolling="yes" seamless="seamless" style="display:block; width:100%; height:65vh; background: rgba(0,0,0,0);"></iframe>
{%- else -%}
<p class="no-data-message">No post-curation data available</p>
{%- endif -%}
</div>
</div>

<br>

## Mission metadata

<div class="row">
<div class="col-sm-6 pre-curation-panel">
<div class="panel-header">Pre-curation</div>
<iframe src="{* datatable_html_path *}" onload='javascript:(function(o){o.style.height=o.contentWindow.document.body.scrollHeight+"px";}(this));' style="height:200px;width:100%;border:none;overflow:hidden;padding:0;"></iframe>
</div>
<div class="col-sm-6 post-curation-panel">
<div class="panel-header">Post-curation</div>
{% if has_post_curation_data -%}
<iframe src="{* post_curation_datatable_html_path *}" onload='javascript:(function(o){o.style.height=o.contentWindow.document.body.scrollHeight+"px";}(this));' style="height:200px;width:100%;border:none;overflow:hidden;padding:0;"></iframe>
{%- else -%}
<p class="no-data-message">No post-curation data available</p>
{%- endif -%}
</div>
</div>

<br>

<!-- JavaScript Toggle Logic -->
<script>
function setView(view) {
    const prePanels = document.querySelectorAll('.pre-curation-panel');
    const postPanels = document.querySelectorAll('.post-curation-panel');
    const buttons = document.querySelectorAll('.view-toggle-btn');

    // Update button states
    buttons.forEach(btn => {
        btn.classList.remove('active', 'btn-primary');
        btn.classList.add('btn-outline-primary');
    });
    // Find and activate the clicked button
    buttons.forEach(btn => {
        if ((view === 'pre' && btn.textContent === 'Pre-curation') ||
            (view === 'post' && btn.textContent === 'Post-curation') ||
            (view === 'both' && btn.textContent === 'Both')) {
            btn.classList.remove('btn-outline-primary');
            btn.classList.add('active', 'btn-primary');
        }
    });

    if (view === 'pre') {
        prePanels.forEach(p => {
            p.style.display = 'block';
            p.className = 'col-12 pre-curation-panel';
        });
        postPanels.forEach(p => p.style.display = 'none');
    } else if (view === 'post') {
        prePanels.forEach(p => p.style.display = 'none');
        postPanels.forEach(p => {
            p.style.display = 'block';
            p.className = 'col-12 post-curation-panel';
        });
    } else { // both
        prePanels.forEach(p => {
            p.style.display = 'block';
            p.className = 'col-sm-6 pre-curation-panel';
        });
        postPanels.forEach(p => {
            p.style.display = 'block';
            p.className = 'col-sm-6 post-curation-panel';
        });
    }
}

// Initialize to 'both' view on page load
document.addEventListener('DOMContentLoaded', function() {
    setView('both');
});
</script>
