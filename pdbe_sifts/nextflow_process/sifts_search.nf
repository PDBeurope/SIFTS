process RUN_SIFTS_SEARCH {

    tag { fasta_file.getName() }

    input:
    path fasta_file

    output:
    path "DONE.txt"

    cpus params.tool.threads
    memory params.tool.memory

    script:
    """
    mkdir -p ${params.work.outdir}
    python /hps/software/users/pdbe/user/adamb/pdbe_sifts/pdbe_sifts/sifts_search_nf.py \
        -i ${params.work.fasta_file} \
        -od ${params.work.outdir} \
        -db ${params.work.db_file} \
        -t ${params.tool.name} \
        -threads ${params.tool.threads}
    
    echo "OK" > DONE.txt
    """
}

workflow {
    Channel.fromPath(params.work.fasta_file).set { fasta_ch }

    RUN_SIFTS_SEARCH(fasta_ch)
}
