process MERGE_REPORTS {
    tag "merge_reports"
    label 'process_low'

    conda "${projectDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-kra:latest' :
        'quay.io/nf-core/kra:latest' }"

    input:
    path sequence_report_dirs, stageAs: 'reports_?/*'
    path all_sequences_json
    path fasta
    path iprscan_results, stageAs: 'iprscan/*'

    output:
    path "sequences/*"                , emit: sequence_dirs
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Create sequence directories structure
    mkdir -p sequences logs

    # The sequence_report_dirs are per-sequence directories (sp-*, tr-* or any other
    # header) containing PF*_report.json files. Each is staged under its own reports_N/
    # parent, because a sequence with domains in several families arrives once per family.
    # Copy all reports to sequences/
    echo "Processing sequence report directories..."
    for src_dir in reports_*/*/; do
        if [ -d "\$src_dir" ]; then
            seq_name=\$(basename "\$src_dir")
            # Skip placeholder directories
            if [ "\$seq_name" = "tr-placeholder-empty" ]; then
                continue
            fi
            mkdir -p "sequences/\$seq_name"
            # Copy all report files
            cp "\$src_dir"/*_report.json "sequences/\$seq_name/" 2>/dev/null || true
        fi
    done

    # Extract sequence list from all_sequences.json and ensure all sequences have directories
    python -c "
import json
import os
import sys

# Load all sequences
with open('${all_sequences_json}', 'r') as f:
    all_seqs = json.load(f)

# Flatten batches to get all sequence IDs
all_sequence_ids = []
for batch, seqs in all_seqs.items():
    all_sequence_ids.extend(seqs)

# Create sequence directories (if they don't exist)
for seq_id in all_sequence_ids:
    dir_name = seq_id.replace('|', '-')
    seq_dir = os.path.join('sequences', dir_name)
    os.makedirs(seq_dir, exist_ok=True)

print(f'Ensured {len(all_sequence_ids)} sequence directories exist', file=sys.stderr)
"

    # Run merge_reports for each sequence
    for seq_dir in sequences/*/; do
        seq_name=\$(basename \$seq_dir)
        python ${projectDir}/merge_reports_in_sequences.py \\
            --sequence "\$seq_name" \\
            --sequence-dir \$seq_dir \\
            --log logs/merge_\${seq_name}.log \\
            ${args} || true
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p sequences/tr-A7K9T9-A7K9T9_9PHYC
    echo '{}' > sequences/tr-A7K9T9-A7K9T9_9PHYC/aggregated_report.json
    echo '>tr|A7K9T9|A7K9T9_9PHYC' > sequences/tr-A7K9T9-A7K9T9_9PHYC/sequence.fasta
    echo 'MSEQVENCE' >> sequences/tr-A7K9T9-A7K9T9_9PHYC/sequence.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d' ' -f2 || echo "stub")
    END_VERSIONS
    """
}
