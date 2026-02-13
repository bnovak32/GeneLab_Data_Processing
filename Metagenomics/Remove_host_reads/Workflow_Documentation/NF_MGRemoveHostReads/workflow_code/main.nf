// main.nf
nextflow.enable.dsl=2

include { paramsHelp         } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

include { KRAKEN2_DB         } from './modules/kraken2_db.nf'
include { KRAKEN_2           } from './modules/kraken2.nf'
include { SUMMARY            } from './modules/summary.nf'
include { COMPILE_SUMMARY    } from './modules/summary.nf'

include { SOFTWARE_VERSIONS  } from './modules/utils.nf'
include { GENERATE_PROTOCOL  } from './modules/generate_protocol.nf'

workflow {
    
    main:

    // check input parameters
    validateParameters()

    // Capture software versions
    software_versions_ch = channel.empty()

    // Get host info
    host_info = channel
        .fromPath(params.hosts_table)
        .splitCsv(header:true)
        .filter { row -> row.name.toLowerCase() == params.host.toLowerCase() }  // match host
        .map { row ->
            def host_id   = row.name.replaceAll(' ', '_').toLowerCase() 
            tuple(row.name, host_id, row.species, row.refseq, row.genome, row.fasta) }

    host_info
    .ifEmpty { error("INPUT ERROR: Host '${params.host}' not found in hosts table '${params.hosts_table}'") }
    
    // Check if kraken2 database already exists or needs to be built
    def host_id = params.host.replaceAll(' ', '_').toLowerCase()
    def host_db = file("${params.ref_dbs_dir}/kraken2-${host_id}-db")
    def db_exists = host_db.exists()

    if (db_exists)
        database_ch = channel.value(host_db)
    else {
        build_ch = host_info.map { name, hostID, species, refseq, genome, fasta -> tuple(name, host_id, fasta) }
        KRAKEN2_DB(build_ch, params.ref_dbs_dir)
        database_ch = KRAKEN2_DB.out.first()
    }
    
    channel
        .fromPath(params.sample_id_list)
        .splitText()
        .map { it.trim() }
        .filter { it } //Ignore blank lines
        .map { sample_id ->
            def meta = [id: sample_id, paired_end: !params.is_single]
                
            def reads = params.is_single ?
                file("$params.reads_dir/${sample_id}$params.single_suffix") :
                [file("$params.reads_dir/${sample_id}$params.R1_suffix"), file("$params.reads_dir/${sample_id}$params.R2_suffix")]

            return tuple(meta, reads)
        }
        .set {generated_reads_ch}

    KRAKEN_2(database_ch, generated_reads_ch, params.out_suffix)
    KRAKEN_2.out.version | mix(software_versions_ch) | set{software_versions_ch}
    
    // Generate summary and compile into one file
    SUMMARY(KRAKEN_2.out.output, KRAKEN_2.out.report)
    COMPILE_SUMMARY(SUMMARY.out.collect(), channel.fromPath(params.sample_id_list), params.host) 

    // Software Version Capturing - combining all captured software versions
    nf_version = "Nextflow Version ".concat("${nextflow.version}")
    nextflow_version_ch = channel.value(nf_version)

    //  Write software versions to file
    software_versions_ch | map { it -> it.text.strip() }
                         | unique
                         | mix(nextflow_version_ch)
                         | collectFile({it -> it}, newLine: true, cache: false)
                         | SOFTWARE_VERSIONS

    // Protocol always needs name, refseq ID, and genome build
    protocol_ch = host_info.map { name, hostID, species, refseq, genome, fasta -> tuple(name, refseq, genome) }
    
    def protocol = host_db.resolve('read-removal-protocol-text.txt')
    protocol_out = GENERATE_PROTOCOL(protocol_ch, SOFTWARE_VERSIONS.out, channel.value(protocol))
    
    publish:
    protocol_out = protocol_out
    software_versions = SOFTWARE_VERSIONS.out
    fastq_out = KRAKEN_2.out.host_removed
    kraken2_out = KRAKEN_2.out.output
    kraken2_report = KRAKEN_2.out.report
    summary_stats = COMPILE_SUMMARY.out.summary_file
    
}

output {
    protocol_out {
        path "processing_info"
    }

    software_versions {
        path "processing_info"
    }

    fastq_out {
        path "${params.reads_outdir}"
    }

    kraken2_out {
        path "results/kraken2-output"
    }

    kraken2_report {
        path "results/kraken2-output"
    }

    summary_stats {
        path "results"
    }
}