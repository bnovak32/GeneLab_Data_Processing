/*
 * DESeq2 Differential Gene Expression Analysis 
 */
// ERCC counts are removed before normalization

process DGE_DESEQ2 {

    input:
        val(meta)
        path(runsheet_path)
        path(gene_counts)
        val(output_prefix)

    output:
        tuple path("Normalized_Counts${params.assay_suffix}.csv"),
              path(params.mode == "microbes" ? "FeatureCounts_Unnormalized_Counts${params.assay_suffix}.csv" : 
                   "RSEM_Unnormalized_Counts${params.assay_suffix}.csv"),                emit: norm_counts
        path("contrasts${params.assay_suffix}.csv"),                                                     emit: contrasts
        path("SampleTable${params.assay_suffix}.csv"),                                                   emit: sample_table      
        path("differential_expression_no_annotations${params.assay_suffix}.csv"),        emit: dge_table
        path("VST_Normalized_Counts${params.assay_suffix}.csv"),                         emit: vst_norm_counts
        path("summary.txt"),                                                                             emit: summary
        path("versions2.txt"),                                                                           emit: versions

    script:
        def output_filename_prefix = output_prefix ?: ""
        def output_filename_suffix = params.assay_suffix ?: ""
        def microbes = params.mode == 'microbes' ? 'TRUE' : 'FALSE'
        def dge_rmd_file = "${projectDir}/bin/dge_deseq2.Rmd"
        def debug_dummy_counts = params.use_dummy_gene_counts ? 'TRUE'  : 'FALSE'

        """
        Rscript -e "rmarkdown::render('${dge_rmd_file}', 
            output_file = 'DGE_DESeq2.html',
            output_dir = '\${PWD}',
            params = list(
                cpus = ${task.cpus},
                work_dir = '\${PWD}',
                output_directory = '\${PWD}',
                output_filename_prefix = '${output_filename_prefix}',
                output_filename_suffix = '${output_filename_suffix}',
                runsheet_path = '${runsheet_path}',
                microbes = ${microbes},
                gene_id_type = '${meta.gene_id_type}',
                input_counts = '${gene_counts}',
                DEBUG_MODE_LIMIT_GENES = FALSE,
                DEBUG_MODE_ADD_DUMMY_COUNTS = ${debug_dummy_counts}
            ))"

        Rscript -e "versions <- c(); 
                    versions['R'] <- gsub(' .*', '', gsub('R version ', '', R.version\\\$version.string));
                    versions['BioConductor'] <- as.character(BiocManager::version()); 
                    pkg_list <- c('BiocParallel', 'DESeq2', 'tidyverse', 'dplyr', 'knitr', 'stringr', 'yaml');
                    for(pkg in pkg_list) {
                        versions[pkg] <- as.character(packageVersion(pkg))
                    };
                    cat('"RNASEQ_DGE_DESEQ2":\\n', 
                        paste0('    ', names(versions), ': ', versions, collapse='\\n'), 
                        '\\n', sep='', file='versions2.txt')"
        """
}