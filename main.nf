#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/kra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Key Residue Annotate (KRA) Pipeline
    A Nextflow DSL2 pipeline for genomic-scale protein function annotation using Pfam domains.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Eduardo-HortaS/Key-Residue-Annotate_Py
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HMMSEARCH            } from './modules/local/hmmsearch'
include { SEQ_BATCH_PREP       } from './modules/local/seq_batch_prep'
include { RUN_IPRSCAN          } from './modules/local/run_iprscan'
include { COLLECT_IPRSCAN      } from './modules/local/collect_iprscan'
include { PREPARE_FASTA        } from './modules/local/prepare_fasta'
include { HMMALIGN             } from './modules/local/hmmalign'
include { TRANSFER_ANNOTATIONS } from './modules/local/transfer_annotations'
include { MERGE_REPORTS        } from './modules/local/merge_reports'
include { MAKE_VIEW_JSONS      } from './modules/local/make_view_jsons'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow KRA {

    take:
    ch_fasta           // channel: [ path(fasta) ]
    ch_hmm             // channel: [ path(hmm) ]
    ch_resources       // channel: [ path(resources_dir) ]

    main:

    ch_versions = Channel.empty()

    // The resources directory is a single item reused by every domain task (HMMALIGN,
    // TRANSFER_ANNOTATIONS): as a queue channel it would let only one domain run.
    ch_resources = ch_resources.first()

    //
    // STEP 1: HMMSEARCH - Search input sequences against Pfam HMM database
    //
    HMMSEARCH (
        ch_fasta,
        ch_hmm
    )
    ch_versions = ch_versions.mix(HMMSEARCH.out.versions)

    //
    // STEP 2: SEQ_BATCH_PREP - Prepare sequence batches for InterProScan
    //
    SEQ_BATCH_PREP (
        ch_fasta,
        HMMSEARCH.out.hmmsearch_sequences_json
    )
    ch_versions = ch_versions.mix(SEQ_BATCH_PREP.out.versions)

    //
    // STEP 3: RUN_IPRSCAN - Run InterProScan on sequence batches
    // Processes batches in parallel, one per batch of all_sequences.json, through
    // run_iprscan.py as executor.py does. Per-sequence results can instead be supplied
    // with params.iprscan_precomputed (<dir>/<sequence>/iprscan.tsv), or skipped.
    //
    if (params.iprscan_precomputed) {
        ch_iprscan_by_sequence = Channel.fromPath(params.iprscan_precomputed, type: 'dir', checkIfExists: true).first()
    } else if (!params.skip_iprscan) {
        ch_batches = SEQ_BATCH_PREP.out.all_sequences_json.flatMap { json ->
            new groovy.json.JsonSlurper().parse(json).collect { name, seqs -> [ (name - 'batch_') as Integer, seqs ] }
        }
        RUN_IPRSCAN (
            ch_batches,
            SEQ_BATCH_PREP.out.sequence_fastas.first()
        )
        ch_versions = ch_versions.mix(RUN_IPRSCAN.out.versions.first())

        // Gather all per-sequence InterProScan results into one directory
        COLLECT_IPRSCAN (
            RUN_IPRSCAN.out.iprscan_results.collect()
        )
        ch_iprscan_by_sequence = COLLECT_IPRSCAN.out.iprscan_by_sequence.first()
    } else {
        // No InterProScan results when skipping
        ch_iprscan_by_sequence = Channel.value([])
    }
    ch_iprscan_results = ch_iprscan_by_sequence

    //
    // STEP 4: PREPARE_FASTA - Prepare per-domain FASTA files
    //
    PREPARE_FASTA (
        HMMSEARCH.out.hmmsearch_per_domain_json,
        ch_fasta,
        ch_resources
    )
    ch_versions = ch_versions.mix(PREPARE_FASTA.out.versions)

    //
    // STEP 5: HMMALIGN - Align sequences against domain HMMs
    // Processes domains in parallel
    //
    HMMALIGN (
        PREPARE_FASTA.out.domain_dirs.flatten(),
        ch_resources
    )
    ch_versions = ch_versions.mix(HMMALIGN.out.versions.first())

    // Collect all alignment results
    ch_aligned_domains = HMMALIGN.out.aligned_domains.collect()

    //
    // STEP 6: TRANSFER_ANNOTATIONS - Transfer annotations from seeds to targets
    // Processes domains in parallel
    //
    TRANSFER_ANNOTATIONS (
        HMMALIGN.out.aligned_domains.flatten(),
        ch_resources,
        ch_iprscan_by_sequence
    )
    ch_versions = ch_versions.mix(TRANSFER_ANNOTATIONS.out.versions.first())

    // Collect all domain reports and sequence directories with reports
    ch_domain_reports = TRANSFER_ANNOTATIONS.out.domain_reports.collect()
    ch_sequence_report_dirs = TRANSFER_ANNOTATIONS.out.sequence_dirs.flatten().collect()

    //
    // STEP 7: MERGE_REPORTS - Merge domain reports into sequence reports
    //
    MERGE_REPORTS (
        ch_sequence_report_dirs,
        SEQ_BATCH_PREP.out.all_sequences_json,
        ch_fasta,
        ch_iprscan_results
    )
    ch_versions = ch_versions.mix(MERGE_REPORTS.out.versions)

    //
    // STEP 8: MAKE_VIEW_JSONS - Create view JSONs for visualization
    //
    MAKE_VIEW_JSONS (
        MERGE_REPORTS.out.sequence_dirs.collect()
    )
    ch_versions = ch_versions.mix(MAKE_VIEW_JSONS.out.versions)

    emit:
    sequence_reports = MAKE_VIEW_JSONS.out.sequence_reports
    domain_reports   = ch_domain_reports
    versions         = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    //
    // VALIDATE INPUT PARAMETERS
    //
    if (!params.fasta) {
        error "Please provide an input FASTA file using --fasta"
    }
    if (!params.hmm) {
        error "Please provide a Pfam HMM file using --hmm"
    }
    if (!params.resources) {
        error "Please provide the resources directory using --resources"
    }

    //
    // CREATE INPUT CHANNELS
    //
    ch_fasta     = Channel.fromPath(params.fasta, checkIfExists: true)
    ch_hmm       = Channel.fromPath(params.hmm, checkIfExists: true)
    ch_resources = Channel.fromPath(params.resources, checkIfExists: true)

    //
    // RUN MAIN WORKFLOW
    //
    KRA (
        ch_fasta,
        ch_hmm,
        ch_resources
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
