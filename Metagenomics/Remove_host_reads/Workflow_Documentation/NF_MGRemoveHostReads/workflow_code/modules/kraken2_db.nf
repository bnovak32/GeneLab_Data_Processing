process KRAKEN2_DB {
    tag "Creating host reads database in ${ref_dbs_dir}"
    publishDir "${ref_dbs_dir}", mode: 'copy'
        
    input:
    tuple val(host_name), val(host_url) , path(host_fasta)
    path ref_dbs_dir

    output:
    path("kraken2-${host_name}-db/"), emit: krakendb_dir
    path("versions.txt"), emit: version

    script:
    if (host_url != null)
        """
        echo "Downloading and unpacking database from ${host_url}
        wget -O ${host_name}.tar.gz --timeout=3600 --tries=0 --continue ${host_url}

        mkdir kraken2-${host_name}-db/ && tar -zxvf -C kraken2-${host_name}-db/

        # Cleaning up
        [ -f  ${host_name}.tar.gz ] && rm -rf ${host_name}.tar.gz

        echo "Kraken2 \$(kraken2 -version | head -n 1 | awk '{print \$3}')" >> versions.txt            
        """
    else if (host_fasta != null)
        """
        echo "Attempting to build a custom ${host_name} reference database from ${host_fasta}"

        # install taxonomy
        k2 download-taxonomy --db kraken2-${host_name}-db/

        # add sequence to database's genomic library
        k2 add-to-library --db kraken2-${host_name}-db/ --threads ${task.cpus} \
                          --files ${host_fasta} --no-masking

        # build the kraken2 database
        k2 build --db kraken2-${host_name}-db/ --threads ${task.cpus} \
                 --kmer-len 35 --minimizer-len 31

        # remove intermediate files
        k2 clean --db kraken2-${host_name}-db/

        echo "Kraken2 \$(kraken2 -version | head -n 1 | awk '{print \$3}')" >> versions.txt            
        """
    else if (host_name != null)
        """
        echo "Download and build kraken reference for named host: ${host_name}"

        # download genomic sequences
        k2 download-library --db kraken2-${host_name}-db/ --threads ${task.cpus} \
                            --library ${host_name} --no-masking

        # install taxonomy
        k2 download-taxonomy --db kraken2-${host_name}-db/

        # build the kraken2 database
        k2 build --db kraken2-${host_name}-db/ --threads ${task.cpus} \
                 --kmer-len 35 --minimizer-len 31

        # remove intermediate files
        k2 clean --db kraken2-${host_name}-db/ 
            
        echo "Kraken2 \$(kraken2 -version | head -n 1 | awk '{print \$3}')" >> versions.txt            
        """
    else
        error "Input error, host_name, host_url, and host_fasta are all set to null. Please supply at least one valid parameter for database creation"


}