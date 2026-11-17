process RUN_IPRSCAN {
    tag "batch_${batch_index}"
    label 'process_high'

    // InterProScan uses local installation - no container for this process
    // Container support for InterProScan will be added in Phase 2

    input:
    tuple val(batch_index), val(batch_sequences)
    path sequence_fastas

    output:
    path "iprscan_by_sequence"        , emit: iprscan_results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def iprscan_path = params.iprscan_path ?: 'interproscan.sh'
    def analyses = params.iprscan_analyses ? "-iA ${params.iprscan_analyses}" : ''
    def output_format = params.iprscan_output_format ?: 'TSV'
    def enable_precalc = params.iprscan_enable_precalc ? '-iDpc' : ''
    def disable_res = params.iprscan_disable_res ? '-iDr' : ''
    def cpu_cores = params.cpu_cores_iprscan ?: task.cpus
    def batch_dash = batch_sequences.collect { it.replace('|', '-') }.join(',')
    """
    mkdir -p iprscan_by_sequence logs

    # run_iprscan.py builds the batch from <parent>/<sequence>/sequence.fasta and writes the
    # split output to <parent>/<sequence>/iprscan.tsv, exactly as executor.py runs it.
    for seq in \$(echo "${batch_dash}" | tr ',' ' '); do
        mkdir -p "iprscan_by_sequence/\$seq"
        cp "${sequence_fastas}/\$seq/sequence.fasta" "iprscan_by_sequence/\$seq/"
    done

    python ${projectDir}/run_iprscan.py \\
        -iPr ${iprscan_path} \\
        -iOf ${output_format} \\
        -iCc ${cpu_cores} \\
        -l logs/iprscan_batch_${batch_index}.log \\
        ${analyses} \\
        ${enable_precalc} \\
        ${disable_res} \\
        ${args} \\
        batch \\
        -sB "${batch_dash}" \\
        -sBi ${batch_index} \\
        -sPd iprscan_by_sequence

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(${iprscan_path} --version 2>&1 | head -1 | grep -oP '[0-9]+\\.[0-9]+-[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    # Create stub per-sequence InterProScan output
    mkdir -p iprscan_by_sequence/tr-A7K9T9-A7K9T9_9PHYC
    echo -e "tr|A7K9T9|A7K9T9_9PHYC\t1234567890abcdef\t450\tPfam\tPF00005\tABC_tran\t410\t540\t1.2E-25\tT\t20-01-2026\tIPR003439\tABC transporter-like\tGO:0005524|GO:0016887" > iprscan_by_sequence/tr-A7K9T9-A7K9T9_9PHYC/iprscan.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: "stub"
    END_VERSIONS
    """
}
