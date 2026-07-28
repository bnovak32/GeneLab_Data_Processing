process GENERATE_PROTOCOL {
  tag "${ params.gldsAccession }"
  publishDir "${ publishdir }/GeneLab",
    mode: params.publish_dir_mode,
    pattern: "*.txt"

  input:
    val(publishdir)
    val(ch_meta)
    path(software_versions_yaml)
    tuple val(ensemblVersion), val(ensemblSource)
    val(bioconductor_annotations)
    path(annotations_db_info)
    val(skipDE)
  
  output:
        path("protocol_GLmicroarray.txt")
  
  script:
    def skipDGE = skipDE ? "--skip-DGE" : ''
    def custom_annot_file = params.annotation_config_path ? "--custom_annot_config ${params.annotation_config_path}" : ''
    def annot_db_info_file = annotations_db_info ? "--annotations_db_info ${annotations_db_info}" : ''
    
    """
    generate_protocol.py \
        --outdir . \
        --software_table ${software_versions_yaml} \
        --assay_suffix "_GLmicroarray" \
        --workflow_version ${workflow.manifest.version} \
        --organism "${ch_meta.organism}" \
        --reference_source ${ensemblSource} \
        --reference_version ${ensemblVersion} \
        --biomart_attribute "${ch_meta.biomart_id}" \
        --bioconductor_annotations ${bioconductor_annotations} \
        $annot_db_info_file \
        $custom_annot_file \
        $skipDGE
    """
}