process GENERATE_PROTOCOL {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating analysis protocol text..."

    input:
        tuple val(host), val(refSeq_ID), val(genome)
        path software_versions
        path protocol

    output:
        path("protocol.txt")
    
    script:
        if (protocol.exists())
        """
        cp ${protocol} protocol.txt
        """
        else
        """
        generate_protocol.sh ${software_versions} ${host} "${refSeq_ID}" ${genome} > protocol.txt
        """
}