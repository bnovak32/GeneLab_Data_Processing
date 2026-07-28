process PROCESS_AGILE1CH {
  publishDir "${publishdir}/GeneLab",
    pattern: "NF_MAAgilent1ch_v${workflow.manifest.version}_GLmicroarray.html",
    mode: params.publish_dir_mode
  stageInMode 'copy'

  input:
    val(publishdir)
    path(qmd) // quarto qmd file to render
    path(runsheet_csv) // runsheet to supply as parameter
    path(array_data_files) // staged, locally-named, decompressed raw array data files
    path(annotation_file_path) // gene annotation file
    tuple val(ensemblVersion), val(ensemblSource)
    path(referenceStorePath) // path to custom annotation references
    path(annotation_config_path) // path to custom annotation config file
    val(skipDE) // whether to skip DE

  output:
    path("NF_MAAgilent1ch_v${workflow.manifest.version}_GLmicroarray.html"), emit: report

    tuple path("02-limma_DGE"),
          path("01-limma_NormExp"),
          path("00-RawData"), emit: de

    path("versions.yml"), emit: versions // Note: Quarto version captured in script body.  R versions captured during render (part of qmd code).

  script:
    def run_DE = skipDE ? "-P run_DE:'false'" : ''
    """
        export HOME=\$PWD;
        
        quarto render \$PWD/${qmd} \
            -P 'workflow_version:${workflow.manifest.version}' \
            -P 'runsheet:${runsheet_csv}' \
            -P 'annotation_file_path:${annotation_file_path}' \
            -P 'ensembl_version:${ensemblVersion}' \
            -P 'local_annotation_dir:${referenceStorePath}' \
            -P 'annotation_config_path:${annotation_config_path}' \
            ${run_DE}

        # Rename report
        mv Agile1CMP.html NF_MAAgilent1ch_v${workflow.manifest.version}_GLmicroarray.html

        cat >> versions.yml <<END_OF_VERSIONS
        - name: quarto
          version: \$(quarto --version)
          homepage: https://quarto.org/
          workflow task: ${task.process}
        END_OF_VERSIONS
    """
}