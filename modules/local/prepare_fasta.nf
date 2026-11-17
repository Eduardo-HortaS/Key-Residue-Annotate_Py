process PREPARE_FASTA {
    tag "prepare_fasta"
    label 'process_low'

    conda "${projectDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-kra:latest' :
        'quay.io/nf-core/kra:latest' }"

    input:
    path hmmsearch_per_domain_json
    path fasta
    path resources

    output:
    path "PF*", type: 'dir'           , emit: domain_dirs
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p logs

    # Extract domain list from hmmsearch_per_domain.json
    domains=\$(python -c "import json; print(' '.join(json.load(open('${hmmsearch_per_domain_json}')).keys()))")

    # Process each domain
    for domain in \$domains; do
        python ${projectDir}/prepare_fasta_per_domain.py \\
            --per-dom-json ${hmmsearch_per_domain_json} \\
            --domain-accession \$domain \\
            --resource-dir ${resources} \\
            --output-dir . \\
            --log logs/prepare_fasta_\${domain}.log \\
            ${args}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    # Create stub domain directories with expected files
    mkdir -p PF00005
    echo '{"can_align": true, "hmm_file": "resources/PF00005/domain.hmm", "seed_alignment": "resources/PF00005/alignment.seed", "conservations": "resources/PF00005/conservations.json", "annotations": "resources/PF00005/annotations.json", "pfam_id_hmmaligned": "PF00005/PF00005_hmmalign.sth", "dom_fasta": "PF00005/PF00005_hits.fasta"}' > PF00005/domain_info.json
    echo '>tr|A7K9T9|A7K9T9_9PHYC/target/410-540' > PF00005/PF00005_hits.fasta
    echo 'LNSTRLNIKRGNHYGLIG' >> PF00005/PF00005_hits.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d' ' -f2 || echo "stub")
    END_VERSIONS
    """
}
