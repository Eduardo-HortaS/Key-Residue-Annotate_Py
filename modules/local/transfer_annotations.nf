process TRANSFER_ANNOTATIONS {
    tag "${domain_dir.baseName}"
    label 'process_medium'

    conda "${projectDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-kra:latest' :
        'quay.io/nf-core/kra:latest' }"

    input:
    path domain_dir
    path resources
    path iprscan_by_sequence, stageAs: 'iprscan_by_sequence'

    output:
    path "${domain_dir}"              , emit: domain_reports
    path "sequence_dirs/*", type: 'dir', optional: true , emit: sequence_dirs
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def domain_accession = domain_dir.baseName
    def eco_codes = params.eco_codes ? "--eco-codes ${params.eco_codes}" : ''
    def extra_donors_flag = params.extra_donors ? '--extra-donors' : ''
    """
    mkdir -p logs

    # transfer_annotations.py reads each target's GO terms from
    # <output-dir>/<sequence>/iprscan.tsv, which executor.py has run_iprscan.py create.
    # Link those files in for this domain's targets.
    if [ -d iprscan_by_sequence ]; then
        grep '^>' ${domain_dir}/${domain_accession}_hits.fasta | sed -e 's/^>//' -e 's/target\\/.*\$//' | tr '|' '-' | sort -u | while read -r seq; do
            if [ -f "iprscan_by_sequence/\$seq/iprscan.tsv" ]; then
                mkdir -p "\$seq"
                ln -s "\$PWD/iprscan_by_sequence/\$seq/iprscan.tsv" "\$seq/iprscan.tsv"
            fi
        done
    fi

    # Run transfer_annotations for this domain
    python ${projectDir}/transfer_annotations.py \\
        --dom-align ${domain_dir}/${domain_accession}_hmmalign.sth \\
        --resource-dir ${resources} \\
        --domain-accession ${domain_accession} \\
        --output-dir . \\
        ${eco_codes} \\
        ${extra_donors_flag} \\
        --log logs/transfer_${domain_accession}.log \\
        ${args}

    # Collect every per-sequence directory (sp-*, tr-* or any other header), not only
    # tr-*. Staged inputs are symlinks and are left where they are.
    mkdir -p sequence_dirs
    for d in */; do
        d=\${d%/}
        [ -d "\$d" ] || continue
        [ -L "\$d" ] && continue
        case "\$d" in logs|sequence_dirs|${domain_accession}) continue ;; esac
        mv "\$d" sequence_dirs/
    done

    # Ensure at least one directory exists (empty placeholder if no hits)
    if [ -z "\$(ls -A sequence_dirs)" ]; then
        mkdir -p sequence_dirs/tr-placeholder-empty
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        goatools: \$(python -c "import goatools; print(goatools.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    def domain_accession = domain_dir.baseName
    """
    # Create stub domain report
    cat <<-EOF > ${domain_dir}/${domain_accession}_report.json
{
    "domain_id": "${domain_accession}",
    "sequences": {
        "tr|A7K9T9|A7K9T9_9PHYC": {
            "hit_intervals": {
                "410-540": {
                    "sequence": "LNSTRLNIKRGNHYGLIG",
                    "length": 131,
                    "hit_start": 410,
                    "hit_end": 540,
                    "annotations": {
                        "positions": {}
                    }
                }
            }
        }
    }
}
EOF

    # Create stub sequence directory with report
    mkdir -p sequence_dirs/tr-A7K9T9-A7K9T9_9PHYC
    cat <<-EOF > sequence_dirs/tr-A7K9T9-A7K9T9_9PHYC/${domain_accession}_report.json
{
    "sequence_id": "tr-A7K9T9-A7K9T9_9PHYC",
    "domain": {
        "${domain_accession}": {
            "hit_intervals": {
                "410-540": {
                    "annotations": {}
                }
            }
        }
    }
}
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d' ' -f2 || echo "stub")
        goatools: "stub"
    END_VERSIONS
    """
}
