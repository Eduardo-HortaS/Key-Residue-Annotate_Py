process COLLECT_IPRSCAN {
    tag "collect_iprscan"
    label 'process_single'

    input:
    path batch_dirs, stageAs: 'batch_?/*'

    output:
    path "iprscan_by_sequence"        , emit: iprscan_by_sequence

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Gather every batch's <sequence>/iprscan.tsv into one directory, the layout
    # transfer_annotations.py expects.
    mkdir -p iprscan_by_sequence
    for seq_dir in batch_*/iprscan_by_sequence/*/; do
        [ -f "\${seq_dir}iprscan.tsv" ] || continue
        seq=\$(basename "\$seq_dir")
        mkdir -p "iprscan_by_sequence/\$seq"
        cp "\${seq_dir}iprscan.tsv" "iprscan_by_sequence/\$seq/"
    done
    """

    stub:
    """
    mkdir -p iprscan_by_sequence
    """
}
