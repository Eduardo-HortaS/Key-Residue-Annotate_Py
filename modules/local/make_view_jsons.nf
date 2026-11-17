process MAKE_VIEW_JSONS {
    tag "make_view_jsons"
    label 'process_medium'

    conda "${projectDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-kra:latest' :
        'quay.io/nf-core/kra:latest' }"

    input:
    path sequence_dirs, stageAs: 'sequences/*'

    output:
    path "sequences/*", includeInputs: true , emit: sequence_reports
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p logs

    # Run make_view_jsons for each sequence
    for seq_dir in sequences/*/; do
        seq_name=\$(basename \$seq_dir)
        python ${projectDir}/make_view_jsons.py \\
            --sequence-dir \$seq_dir \\
            --sequence "\$seq_name" \\
            --log logs/make_views_\${seq_name}.log \\
            ${args} || true
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    # Stub mode - create sequence directories with expected outputs
    mkdir -p sequences/tr-A7K9T9-A7K9T9_9PHYC
    echo '{}' > sequences/tr-A7K9T9-A7K9T9_9PHYC/aggregated_report.json
    echo '[]' > sequences/tr-A7K9T9-A7K9T9_9PHYC/PF00005_ranges.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d' ' -f2 || echo "stub")
    END_VERSIONS
    """
}
