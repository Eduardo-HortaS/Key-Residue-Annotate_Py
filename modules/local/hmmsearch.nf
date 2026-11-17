process HMMSEARCH {
    tag "hmmsearch"
    label 'process_medium'

    conda "${projectDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-kra:latest' :
        'quay.io/nf-core/kra:latest' }"

    input:
    path fasta
    path hmm

    output:
    path "hmmsearch_per_domain.json"  , emit: hmmsearch_per_domain_json
    path "hmmsearch_sequences.json"   , emit: hmmsearch_sequences_json
    path "hmmsearch_sequences.txt"    , emit: hmmsearch_sequences_txt
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def nucleotide_flag = params.nucleotide ? '--nucleotide' : ''
    def bit_cutoffs = params.bit_cutoffs ?: 'gathering'
    """
    mkdir -p logs
    python ${projectDir}/run_hmmsearch.py \\
        --fasta ${fasta} \\
        --hmm ${hmm} \\
        --output-dir . \\
        ${nucleotide_flag} \\
        --bit-cutoffs ${bit_cutoffs} \\
        --log logs/hmmsearch.log \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        pyhmmer: \$(python -c "import pyhmmer; print(pyhmmer.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo '{"PF00005": {"tr|A7K9T9|A7K9T9_9PHYC": [{"hmm_name": "ABC_tran", "target_seq_name": "tr|A7K9T9|A7K9T9_9PHYC", "bitscore": 86.0, "ali_from": 410, "ali_to": 540, "ali_range": "/410-540", "subseq": "LNSTRLNIKRGNHYGLIG"}]}}' > hmmsearch_per_domain.json
    echo '{"sequences": ["tr|A7K9T9|A7K9T9_9PHYC", "tr|A7K8Q5|A7K8Q5_9PHYC"]}' > hmmsearch_sequences.json
    echo 'tr|A7K9T9|A7K9T9_9PHYC' > hmmsearch_sequences.txt
    echo 'tr|A7K8Q5|A7K8Q5_9PHYC' >> hmmsearch_sequences.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d' ' -f2 || echo "stub")
        pyhmmer: "stub"
    END_VERSIONS
    """
}
