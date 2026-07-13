nextflow.enable.dsl=2

include { paramsHelp } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { paramsSummaryLog } from 'plugin/nf-schema'

include { PARSE_ANNOTATION_TABLE } from './modules/PARSE_ANNOTATION_TABLE.nf'
include { VV_AGILE1CH } from './modules/VV_AGILE1CH.nf'
include { PROCESS_AGILE1CH } from './modules/PROCESS_AGILE1CH.nf'
include { RUNSHEET_FROM_GLDS } from './modules/RUNSHEET_FROM_GLDS.nf'
include { RUNSHEET_FROM_ISA } from './modules/RUNSHEET_FROM_ISA.nf'
include { GENERATE_SOFTWARE_TABLE } from './modules/GENERATE_SOFTWARE_TABLE'
include { DUMP_META } from './modules/DUMP_META'
include { GENERATE_PROTOCOL } from './modules/POST_PROCESSING/GENERATE_PROTOCOL'

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/

workflow {
	main:

    // color defs
    c_bright_green = "\u001b[32;1m";
    c_reset = "\033[0m";

    /**************************************************
    * HELP MENU  **************************************
    **************************************************/
    if ( params.help ) {
      before_text = """
************************************************
* Microarray Agilent 1 Channel Pipeline: ${workflow.manifest.version} *
************************************************

Usage example 1: Processing OSDR datasets
    > nextflow run ./main.nf --osdAccession OSD-548 --gldsAccession GLDS-548

Usage example 2: Processing Other datasets (requires a user-created runsheet)
    > nextflow run ./main.nf --runsheetPath </path/to/runsheet>


"""
      log.info paramsHelp(
        beforeText: before_text, 
        afterText: "For more information, please see the README.md file in the workflow code directory.",
        fullHelp: true)
      exit(0)
    }

    // validate parameters and print parameter summary log (includes only parameters that are not set to default values)
    validateParameters(cast_cli_params: true)
    log.info paramsSummaryLog(workflow)

    /**************************************************
    * CHECK REQUIRED PARAMS AND LOAD  *****************
    **************************************************/
    log.info("Resolved output directory: \"${ params.resultsDir }")

    // Catch case where only one of OSD/GLDS is set or where neither is set and runsheetPath is also not set
    if (params.gldsAccession && !params.osdAccession) {
      error("ERROR: GLDS accession set but OSD accession is not set.  Please set both or neither.")
    } else if (params.gldsAccession && !params.osdAccession) {
      error("ERROR: OSD accession set but GLDS accession is not set.  Please set both or neither.")
    } else if ((!params.gldsAccession && !params.osdAccession) && (!params.runsheetPath && !params.isaArchivePath)) {
      error("ERROR: Neither OSD/GLDS accessions nor runsheetPath or isaArchivePath are set.  Please set either OSD/GLDS accessions, runsheetPath, or isaArchivePath.")
    }

    if ( !params.runsheetPath && !params.isaArchivePath) {
        RUNSHEET_FROM_GLDS( 
          params.osdAccession,
          params.gldsAccession,
          "${ projectDir }/bin/dp_tools__agilent_1_channel" // dp_tools plugin location
        ) 
        RUNSHEET_FROM_GLDS.out.runsheet | set{ ch_runsheet }
    } else if ( !params.runsheetPath && params.isaArchivePath ) {
        RUNSHEET_FROM_ISA( 
          params.osdAccession,
          params.gldsAccession,
          params.isaArchivePath,
          "${ projectDir }/bin/dp_tools__agilent_1_channel" // dp_tools plugin location
        )
        RUNSHEET_FROM_ISA.out.runsheet | set{ ch_runsheet }
    } else if ( params.runsheetPath && !params.isaArchivePath ) {
        ch_runsheet = channel.fromPath( params.runsheetPath )
    } else if ( params.runsheetPath && params.isaArchivePath ) {
        error("Error: User supplied both runsheetPath and isaArchivePath.  Only one or neither is allowed to be supplied!")
    }

    ch_runsheet | splitCsv(header: true) | first |  set{ ch_meta }

    PARSE_ANNOTATION_TABLE(
        params.annotation_file_path,
        ch_meta | map { it -> it.organism }
    )

    PROCESS_AGILE1CH(
      channel.fromPath( "${ projectDir }/bin/Agile1CMP.qmd" ),
      ch_runsheet,
      PARSE_ANNOTATION_TABLE.out.annotations_db_url,
      PARSE_ANNOTATION_TABLE.out.reference_version_and_source,
      params.limit_biomart_query,
      params.skipDE
    )

    VV_AGILE1CH( 
      ch_runsheet, 
      PROCESS_AGILE1CH.out.de,
      params.skipVV,
      "${ projectDir }/bin/${ params.skipDE ? 'dp_tools__agilent_1_channel_skipDE' : 'dp_tools__agilent_1_channel' }" // dp_tools plugin
      )

    // Software Version Capturing
    nf_version = "- name: nextflow\n  ".concat(
"""
  version: ${nextflow.version}
  homepage: https://www.nextflow.io
  workflow task: N/A
""")
    ch_software_versions = channel.value(nf_version)
    PROCESS_AGILE1CH.out.versions | map{ it -> it.text } | mix(ch_software_versions) | set{ch_software_versions}
    VV_AGILE1CH.out.versions | map{ it -> it.text } | mix(ch_software_versions) | set{ch_software_versions}

    GENERATE_SOFTWARE_TABLE(
      ch_software_versions | unique | collectFile(newLine: true, sort: true, cache: false),
      ch_runsheet | splitCsv(header: true, quote: '"') | first | map{ row -> row['Array Data File Name'] },
      params.skipDE
    )

    // export meta for post processing usage
    ch_meta | DUMP_META

    GENERATE_PROTOCOL(
      ch_meta,
      ch_software_versions | unique | collectFile(newLine: true, sort: true, cache: false),
      PARSE_ANNOTATION_TABLE.out.reference_version_and_source,
      PARSE_ANNOTATION_TABLE.out.bioconductor_annotations,
      PARSE_ANNOTATION_TABLE.out.annotations_db_info_url,
      params.skipDE
    )


  workflow.onComplete = { 
    println "${c_bright_green}Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if ( workflow.success ) {
      println "Raw and Processed data location: ${ params.resultsDir }"
      println "V&V logs location: ${ params.resultsDir }/VV_Logs"
      println "Pipeline tracing/visualization files location:  ${ params.resultsDir }/Resource_Usage${c_reset}"
    }
  }

}