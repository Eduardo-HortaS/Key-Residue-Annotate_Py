process HMMALIGN {
    tag "${domain_dir.baseName}"
    label 'process_low'

    conda "${projectDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-kra:latest' :
        'quay.io/nf-core/kra:latest' }"

    input:
    path domain_dir
    path resources

    output:
    path "${domain_dir}"              , emit: aligned_domains
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def domain_accession = domain_dir.baseName
    def trim_flag = params.trim ? '--trim' : ''
    """
    mkdir -p logs

    # Run hmmalign for this domain
    python ${projectDir}/run_hmmalign.py \\
        --dom-info ${domain_dir}/domain_info.json \\
        --domain-accession ${domain_accession} \\
        ${trim_flag} \\
        --log logs/hmmalign_${domain_accession}.log \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        hmmer: \$(hmmalign -h | head -2 | tail -1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def domain_accession = domain_dir.baseName
    """
    # Create stub alignment file
    cat <<-EOF > ${domain_dir}/${domain_accession}_hmmalign.sth
# STOCKHOLM 1.0
#=GF ID   ${domain_accession}
#=GF AC   ${domain_accession}
tr|A7K9T9|A7K9T9_9PHYC/target/410-540    LNSTRLNIKRGNHYGLIG
//
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d' ' -f2 || echo "stub")
        hmmer: "stub"
    END_VERSIONS
    """
}
