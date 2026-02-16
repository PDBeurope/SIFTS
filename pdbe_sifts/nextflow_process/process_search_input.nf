nextflow.enable.dsl=2

params.mmcif_list = "test_entries.txt"
params.out_fasta  = "test_entries.fasta"

process MMCIF_TO_FASTA {

    tag { pdb_id }

    cpus 30
    memory '64 GB'

    input:
    path cif

    output:
    path "${pdb_id}.fasta"

    script:
    pdb_id = cif.baseName

    """
    python /hps/software/users/pdbe/user/adamb/pdbe_sifts/pdbe_sifts/sifts_mmcif2fasta.py \
        --cif ${cif} \
        --out ${pdb_id}.fasta
    """
}

process MERGE_FASTA {

    cpus 1
    memory '24 GB'

    input:
    path fastas

    output:
    path params.out_fasta

    script:
    """
    cat ${fastas.join(' ')} > ${params.out_fasta}
    """
}

workflow {

    Channel
        .fromPath(params.mmcif_list)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { file(it) }
        .set { mmcifs }

    mmcifs | MMCIF_TO_FASTA | collect | MERGE_FASTA
}
