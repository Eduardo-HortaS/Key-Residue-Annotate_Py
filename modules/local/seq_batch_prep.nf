process SEQ_BATCH_PREP {
    tag "seq_batch_prep"
    label 'process_single'

    conda "${projectDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-kra:latest' :
        'quay.io/nf-core/kra:latest' }"

    input:
    path fasta
    path hmmsearch_sequences_json

    output:
    path "all_sequences.json"         , emit: all_sequences_json
    path "batches/"                   , emit: batches
    path "sequence_fastas"            , emit: sequence_fastas
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def nucleotide_flag = params.nucleotide ? '--nucleotide' : ''
    def batch_size = params.seq_batch_size ?: 2000
    """
    # Create output structure
    mkdir -p batches logs

    python ${projectDir}/seq_and_batch_prep.py \\
        --fasta ${fasta} \\
        --output-dir . \\
        --batch-size ${batch_size} \\
        ${nucleotide_flag} \\
        --log logs/seq_batch_prep.log \\
        ${args}

    # Move batch files to batches directory
    mv batch_*.fasta batches/ 2>/dev/null || true

    # seq_and_batch_prep.py writes <sequence>/sequence.fasta for every input sequence;
    # run_iprscan.py builds its batches from them, as in executor.py.
    mkdir -p sequence_fastas
    for d in */; do
        d=\${d%/}
        [ -d "\$d" ] || continue
        [ -L "\$d" ] && continue
        case "\$d" in batches|logs|sequence_fastas) continue ;; esac
        mv "\$d" sequence_fastas/
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p batches sequence_fastas/tr-A7K9T9-A7K9T9_9PHYC
    echo '>tr|A7K9T9|A7K9T9_9PHYC' > sequence_fastas/tr-A7K9T9-A7K9T9_9PHYC/sequence.fasta
    echo 'MSEQVENC' >> sequence_fastas/tr-A7K9T9-A7K9T9_9PHYC/sequence.fasta
    echo '{"batch_1": ["tr|A7K9T9|A7K9T9_9PHYC", "tr|A7K8Q5|A7K8Q5_9PHYC"]}' > all_sequences.json
    echo '>tr|A7K9T9|A7K9T9_9PHYC' > batches/batch_1_tr-A7K9T9-A7K9T9_9PHYC_to_tr-A7K8Q5-A7K8Q5_9PHYC.fasta
    echo 'MSEQVENC' >> batches/batch_1_tr-A7K9T9-A7K9T9_9PHYC_to_tr-A7K8Q5-A7K8Q5_9PHYC.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d' ' -f2 || echo "stub")
        biopython: "stub"
    END_VERSIONS
    """
}
