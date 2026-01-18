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

## Image locations and attributes

<iframe src="{* map_html_path *}" frameborder="0" scrolling="yes" seamless="seamless" style="display:block; width:100%; height:75vh; background: rgba(0,0,0,0);" class="tester"></iframe>

<br>

## Mission metadata

<iframe src="{* datatable_html_path *}" onload='javascript:(function(o){o.style.height=o.contentWindow.document.body.scrollHeight+"px";}(this));' style="height:200px;width:100%;border:none;overflow:hidden;padding:0;"></iframe>

<br>

<!-- Script to make the datatable the height to fit the data -->
<script type="application/javascript">
    var iframe = document.getElementById("myIframe");

    iframe.onload = function(){
    iframe.contentWindow.document.body.scrollHeight + 'px';
    }
</script>
