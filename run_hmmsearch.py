"""
run_hmmsearch.py

Copyright 2025 Eduardo Horta Santos <GitHub: Eduardo-HortaS>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.

This script contains all functions dealing directly with hmmsearch and needs 3 arguments:
a FASTA file, a HMM targets file and the path to the output directory. Optionally,
a flag to indicate nucleotide sequences, bit score cutoffs for reporting hits, and a log file path.

1 - Runs pyHMMER hmmsearch on the FASTA file using the provided HMM database file,
generates a 'hmmsearch_per_domain.json' file in the output directory.
Also, generates a 'hmmsearch_sequences.txt' and a 'hmmsearch_sequences.json'
file in the output directory, these last two contain the sequence IDs with at least 1 domain hit.

    1.5 - Loads the sequence file, either as a SequenceFile or a DigitalSequenceBlock,
    depending on size and available memory. For nucleotide sequences, performs translation to protein sequences before searching.

This script assumes that the user will provide HMM profiles from a curated database,
where specific bit score thresholds for each profile should be present,
including both per-sequence and per-domain reporting and inclusion thresholds.
Particularly, the HMMs should come from Pfam, ensuring that the default choice - "gathering" - is available.
We also accept "trusted" and "noise" bit score cutoffs, defined as follows:
    - gathering: Used to define what gets included in Pfam Full alignments based on searches with Pfam seed models. Reliable curated thresholds for defining family membership.
    - trusted: Represents the score of the lowest-scoring known true positive in a family that is above all known false positives.
    - noise: Represents the score of the highest-scoring known false positive for the family.
"""

import os
import json
import argparse
import logging
import sys
from typing import Union
import psutil
import pyhmmer
from pyhmmer.easel import DigitalSequenceBlock, DigitalSequence
from utils import get_logger
# from modules.decorators import measure_time_and_memory

def parse_arguments():
    """Parse command-line arguments for running hmmsearch
    using a sequence database against target HMMs,
    aiming to generate a hmmsearch_output.json.
    Also includes an output directory and a log file path.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description=
    'Runs hmmsearch using a sequence database against target HMMs \
    aiming to generate a hmmsearch_per_domain.json file.')
    parser.add_argument("-iF", "--fasta", help="Path to fasta file",
                        required=True, type=str)
    parser.add_argument("-iH", "--hmm", help="Path to target HMMs database file",
                        required=True, type=str)
    parser.add_argument("-o", "--output-dir", help="Output dir path",
                        required=True, type=str)
    parser.add_argument("-n", "--nucleotide",
                        help="Flag to indicate nucleotide instead of default protein sequences",
                        action="store_true")
    parser.add_argument("-bc", "--bit-cutoffs",
                        help="Bit score cutoffs for reporting hits. Options: 'noise', 'gathering', 'trusted'",
                        required=False, type=str, default="gathering")
    parser.add_argument("-l", "--log",
                        help="Log path",
                        required=False, type=str, default="logs/run_hmmsearch.log")
    return parser.parse_args()

def load_and_translate_sequence_file(fasta_path: str, logger: logging.Logger, is_nucleotide: bool = False) -> Union[DigitalSequenceBlock, DigitalSequence]:
    """
    Loads a multifasta file either as a SequenceFile or a DigitalSequenceBlock,
    the latter only if the file loaded to memory would take less than 20% of the available memory.

    Args:
        fasta_path: Path to the multifasta file
        logger: Logger instance for tracking execution
        is_nucleotide: If True, treats input as nucleotide sequences (default: False)

    Returns:
        Union[DigitalSequenceBlock, DigitalSequence]: Loaded sequence file
    """
    available_memory = psutil.virtual_memory().available
    target_size= os.stat(fasta_path).st_size
    logger.info(f"RUN_HMMSEARCH --- LOAD_TRANSLATE --- Available memory: {available_memory/1024:.1f} KiB")
    logger.info(f"RUN_HMMSEARCH --- LOAD_TRANSLATE --- Database on-disk size: {target_size/1024:.1f} KiB")

    if is_nucleotide:
        alphabet = pyhmmer.easel.Alphabet.dna()
    else:
        alphabet = pyhmmer.easel.Alphabet.amino()

    with pyhmmer.easel.SequenceFile(fasta_path, digital=True, alphabet=alphabet) as seq_file:
        if target_size < available_memory * 0.2:
            logger.info("RUN_HMMSEARCH --- LOAD_TRANSLATE --- Pre-fetching targets into memory")
            targets = seq_file.read_block()
            if is_nucleotide:
                logger.info("RUN_HMMSEARCH --- LOAD_TRANSLATE --- Translating nucleotide sequences to protein sequences")
                targets = targets.translate()
            logger.info(f"RUN_HMMSEARCH --- LOAD_TRANSLATE --- Database in-memory size: {(sys.getsizeof(targets) + sum(sys.getsizeof(target) for target in targets))/1024:.1f} KiB")
        else:
            targets = seq_file
            if is_nucleotide:
                logger.info("RUN_HMMSEARCH --- LOAD_TRANSLATE --- Translating nucleotide sequences to protein sequences")
                targets = targets.translate()
    return targets

def run_hmmsearch(hmm: str, fasta_path: str, output_dir: str, logger: logging.Logger, bit_cutoffs: str = "gathering", is_nucleotide: bool = False) -> None:
    """Run HMMER search against target sequences and save results.

    Executes hmmsearch using HMM profiles as queries against target sequences.
    Processes hits to extract domain information and saves results in multiple formats.

    Args:
        hmm: Path to HMM profiles database file
        fasta_path: Path to target sequences FASTA file
        output_dir: Directory to save output files
        logger: Logger instance for tracking execution
        bit_cutoffs: Bit score cutoffs for reporting hits ("noise", "gathering", or "trusted")
        is_nucleotide: If True, treats input as nucleotide sequences (default: False)

    Returns:
        set[str]: Set of sequence IDs that had at least one domain hit

    Outputs:
        - hmmsearch_per_domain.json: JSON file containing detailed domain hits
          Structure: {pfam_id: {seq_id: [{seq_hits_data}]}}
        - hmmsearch_sequences.txt: Plain text file with hit sequence IDs
        - hmmsearch_sequences.json: JSON file with hit sequence IDs

    Domain data includes:
        - hmm_name: Name of the matching HMM profile
        - target_seq_name: ID of the target sequence
        - bitscore: HMMER bit score for the match
        - ali_from: Start position of alignment in target
        - ali_to: End position of alignment in target
        - ali_range: Formatted string of alignment range
        - subseq: Aligned subsequence without gaps
    """
    valid_cutoffs = ["noise", "gathering", "trusted"]
    if bit_cutoffs not in valid_cutoffs:
        logger.warning("RUN_HMMSEARCH --- RUN --- Invalid bit_cutoffs value: '%s'. Using default 'gathering'.", bit_cutoffs)
        bit_cutoffs = "gathering"

    os.makedirs(output_dir, exist_ok=True)

    with open(hmm, 'rb') as f:
        hmms = list(pyhmmer.plan7.HMMFile(f))

    targets = load_and_translate_sequence_file(fasta_path, logger, is_nucleotide)

    hits_per_domain = {}
    hit_sequences = set()

    for top_hits in pyhmmer.hmmsearch(hmms, targets, bit_cutoffs=bit_cutoffs):
        for hit in top_hits:
            target_seq = hit.name.decode("utf-8")
            hit_sequences.add(target_seq)
            for domain in hit.domains.included:
                alignment = domain.alignment
                hmm_name = alignment.hmm_name.decode("utf-8")
                # Accession is the Pfam ID
                accession = alignment.hmm_accession.decode("utf-8").split(".")[0] if alignment.hmm_accession.decode("utf-8") != hmm_name else None
                if accession is None:
                    continue
                if accession not in hits_per_domain:
                    hits_per_domain[accession] = {}
                ali_from_1 = alignment.target_from
                ali_to_1 = alignment.target_to
                dirty_subseq = alignment.target_sequence
                subseq = dirty_subseq.translate(str.maketrans('', '', '-_'))  # Remove gaps and insertions from the subsequence
                ali_range = f"/{ali_from_1}-{ali_to_1}"
                if target_seq not in hits_per_domain[accession] and target_seq not in [None, ""]:
                    hits_per_domain[accession][target_seq] = []
                hits_per_domain[accession][target_seq].append({
                    "hmm_name": hmm_name,
                    "target_seq_name": target_seq,
                    "bitscore": hit.score,
                    "ali_from": ali_from_1,
                    "ali_to": ali_to_1,
                    "ali_range": ali_range,
                    "subseq": subseq
                })

    sequences_txt_path = os.path.join(output_dir, "hmmsearch_sequences.txt")
    sequences_json_path = os.path.join(output_dir, "hmmsearch_sequences.json")
    per_domain_output = os.path.join(output_dir, "hmmsearch_per_domain.json")

    # Text file (human readable, grep and so on)
    with open(sequences_txt_path, "w", encoding='utf-8') as f:
        for seq in sorted(hit_sequences):
            f.write(f"{seq}\n")

    # JSON file (programmatic access)
    with open(sequences_json_path, "w", encoding='utf-8') as f:
        json.dump({"sequences": list(sorted(hit_sequences))}, f, indent=4)

    with open(per_domain_output, "w", encoding='utf-8') as f:
        json.dump(hits_per_domain, f, indent=4)

    logger.info(f"RUN_HMMSEARCH --- RUN --- HmmSearch hit sequences saved in text format - {sequences_txt_path}")
    logger.info(f"RUN_HMMSEARCH --- RUN --- HmmSearch hit sequences saved in JSON format - {sequences_json_path}")
    logger.info(f"RUN_HMMSEARCH --- RUN --- HmmSearch TopHits results saved per domain - {per_domain_output}")

def main():
    """Main function, initializes this script"""
    args = parse_arguments()
    logger, _ = get_logger(args.log, scope="main")
    input_fasta = args.fasta
    input_hmm = args.hmm
    output_dir = args.output_dir
    bit_cutoffs = args.bit_cutoffs
    is_nucleotide = args.nucleotide
    logger.info("RUN_HMMSEARCH --- MAIN --- Running hmmsearch with arguments: %s", args)

    # Run hmmsearch for all sequences
    run_hmmsearch(input_hmm, input_fasta, output_dir, logger, bit_cutoffs, is_nucleotide)

if __name__ == '__main__':
    main()
