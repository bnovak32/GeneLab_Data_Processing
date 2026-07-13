nextflow.enable.dsl=2

include { GENERATE_MD5SUMS } from './modules/GENERATE_MD5SUMS.nf'
include { UPDATE_ISA_TABLES } from './modules/UPDATE_ISA_TABLES.nf'

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/

workflow {

  
  main:

    // color defs
    c_back_bright_red = "\u001b[41;1m";
    c_reset = "\033[0m";


    /**************************************************
    * HELP MENU  **************************************
    **************************************************/
    if (params.help) {
      println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
      println("┇ Microarray Agilent 1 Channel Post Processing Pipeline: $workflow.manifest.version  ┇")
      println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
      println("Post processing workflow. Generates md5sum of output files and updates ISA archive tables. Help menu refinements to come")
      exit 0
      }

    println "PARAMS: $params"
    println "\n"

    /**************************************************
    * CHECK REQUIRED PARAMS AND LOAD  *****************
    **************************************************/
    println("Resolved output directory: ${ params.resultsDir }")

    ch_processed_directory = channel.fromPath("${ params.resultsDir }", checkIfExists: true)
    ch_runsheet = channel.fromPath("${ params.resultsDir }/Metadata/*_runsheet.csv", checkIfExists: true)

    GENERATE_MD5SUMS(      
      ch_processed_directory, 
      ch_runsheet,       
      "${ projectDir }/bin/${ params.skipDE ? 'dp_tools__agilent_1_channel_skipDE' : 'dp_tools__agilent_1_channel' }" // dp_tools plugin
    )

    def isa_file = file("${ params.resultsDir }/Metadata/*ISA*.zip")
    if ( isa_file ) {
      ch_isa = channel.fromPath("${ params.resultsDir }/Metadata/*ISA*.zip")
      UPDATE_ISA_TABLES(
        ch_processed_directory, 
        ch_runsheet,
        ch_isa,
        "${ projectDir }/bin/${ params.skipDE ? 'dp_tools__agilent_1_channel_skipDE' : 'dp_tools__agilent_1_channel' }"
      )
    } else {
      println "${ c_back_bright_red }WARNING: No ISA archive found in ${ params.resultsDir }/Metadata/ -- skipping UPDATE_ISA_TABLES${ c_reset }"
    }
    
}