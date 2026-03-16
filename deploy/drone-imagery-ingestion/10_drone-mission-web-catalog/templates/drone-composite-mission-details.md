---
title: Composite Mission {* composite_id *}

date:
show_date: false
profile: false
---

<div class="container">
    <div class="row">
        <div class="col-sm">
            <div class="text-center">
            {{< cta cta_text="Previous composite" cta_link="{* previous_dataset_page_path *}" cta_new_tab="false" >}}
            </div>
        </div>
        <div class="col-sm">
            <div class="text-center">
            {{< cta cta_text="Next composite" cta_link="{* next_dataset_page_path *}" cta_new_tab="false" >}}
            </div>
        </div>
    </div>
</div>

<p style="line-height: 125%; text-align:center;"><b>This page is in beta.</b> The drone data web catalog is under active development and will continue to improve. Feel free to <a href="/about/#contact-us">contact us</a> with feedback!</p>

{% if withhold_from_training -%}
<p style="color: red; line-height: 125%; text-align:center;"><b>DATASET WITHHELD FROM ML TRAINING.</b> We have designated a stratified subset of approximately 20% of OFO datasets to be excluded from training ML models that are intended for use by the forest mapping/ML community. This is a withheld plot. For developing ML models intended for broad use, this plot is reserved for use as test data for model evaluation and intercomparison.</p>

{% endif -%}

## Constituent missions

This composite mission combines two overlapping drone missions flown at different altitudes:

- [Mission {* mission_id_a *}]({* individual_mission_page_path_a *}) (higher altitude, nadir camera pitch)
- [Mission {* mission_id_b *}]({* individual_mission_page_path_b *}) (lower altitude, oblique camera pitch)

<br>

{% if ttops_exists -%}

## Detected trees

<iframe src="{* itd_map_path *}" frameborder="0" scrolling="yes" seamless="seamless" style="display:block; width:100%; height:75vh; background: rgba(0,0,0,0);" class="tester"></iframe>

[Download tree points]({* ttops_url *})

<br>

{% endif -%}

{% if ortho_exists -%}

## Orthomosaic

{{< figure src="{* ortho_url_thumb *}" caption="[Download full orthomosaic]({* ortho_url_full *})" >}}

<br>

{% endif -%}

{% if chm_exists -%}

## Canopy height model

{{< figure src="{* chm_url_thumb *}" caption="[Download full CHM]({* chm_url_full *})" >}}

<br>

{% endif -%}

{% if dsm_exists -%}

## Digital surface model

{{< figure src="{* dsm_url_thumb *}" caption="[Download full DSM]({* dsm_url_full *})" >}}

<br>

{% endif -%}

{% if dtm_exists -%}

## Digital terrain model

{{< figure src="{* dtm_url_thumb *}" caption="[Download full DTM]({* dtm_url_full *})" >}}

<br>

{% endif -%}

{% if pc_exists -%}

## Point cloud

Preview in development. For now, you can paste [this url]({* pc_url_full *}) into a point cloud viewer like [Eptium](https://viewer.copc.io/).

[Download full point cloud]({* pc_url_full *})

<br>

{% endif -%}

{% if mesh_exists -%}

## Mesh model

Preview in development.

[Download full mesh model]({* mesh_url_full *})

<br>

{% endif -%}

## Raw drone images

Access the raw drone images via the constituent missions:

- [Mission {* mission_id_a *}]({* individual_mission_page_path_a *}) (higher altitude, nadir camera pitch)
- [Mission {* mission_id_b *}]({* individual_mission_page_path_b *}) (lower altitude, oblique camera pitch)

Note that mission composites may use a subset of the images from the each of the constituent missions (focusing on the area where the two missions overlap). {% if image_metadata_exists -%} The images that are retained from each consituent mission are indicated in the [image metadata file]({* image_metadata_url *}).{% endif %}

<br>

## Mission attributes

<iframe src="{* composite_details_datatable_path *}" onload='javascript:(function(o){o.style.height=o.contentWindow.document.body.scrollHeight+"px";}(this));' style="height:200px;width:100%;border:none;overflow:hidden;padding:0;"></iframe>

<br>

## Image locations

<iframe src="{* composite_details_map_path *}" frameborder="0" scrolling="yes" seamless="seamless" style="display:block; width:100%; height:75vh; background: rgba(0,0,0,0);" class="tester"></iframe>

<br>

{% if image_metadata_exists or footprint_exists or cameras_exists or log_exists -%}

## Other data

{% if image_metadata_exists -%}
[Composite mission image metadata]({* image_metadata_url *})

{% endif -%}

{% if footprint_exists -%}
[Composite mission geospatial footprint]({* footprint_url *})

{% endif -%}

{% if cameras_exists -%}
[Photogrammetry camera locations and parameters]({* cameras_url *})

{% endif -%}

{% if log_exists -%}
[Photogrammetry processing log]({* log_url *})

{% endif -%}

{% endif -%}


<!-- Script to make the datatable the height to fit the data -->
<script type="application/javascript">
    var iframe = document.getElementById("myIframe");
    iframe.onload = function(){
        iframe.style.height = iframe.contentWindow.document.body.scrollHeight + 'px';
    }
</script>
