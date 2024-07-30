"""
hmmsearch.py

Copyright 2024 Eduardo Horta Santos <GitHub: Eduardo-HortaS>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
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
a FASTA file, a HMM targets file and the path to the output directory.

# Moved to these general functions to utils.py, for better modularity:
    1 - Parses the FASTA file and stores each sequence's queryname (ex.: TTHY_XENLA) in a list,
    and generates a SeqRecord object to be used sequentially and *only once*.

    2 - For each sequence, creates a directory
    and writes the sequence to a 'sequence.fasta' file in the directory.

3 - Runs pyhmmer hmmsearch on the FASTA file using the provided HMM database file,
generates a 'hmmsearch_per_domain_output.json' file, one level above all sequence dirs.

    3.5 - Loads the sequence file, either as a SequenceFile or a DigitalSequenceBlock,
    depending on size and available memory.

"""

import os
import json
import argparse
import logging
import sys
import psutil
import pyhmmer
import utils
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
    parser.add_argument("-iF", "--fasta", help="Path to fasta file", required=True, type=str)
    parser.add_argument("-iH", "--hmm", help="Path to target HMMs database file", required=True, type=str)
    parser.add_argument("-o", "--output-dir", help="Output dir path", required=True, type=str)
    parser.add_argument("-l", "--log", help="Log path", \
        required=False, type=str, default="logs/run_hmmsearch.log")
    return parser.parse_args()

def setup_logging(log_path):
    """Set up logging for the script."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(filename=log_path, level=logging.DEBUG, filemode='w', \
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def load_sequence_file(fasta_path):
    """
    Loads a multifasta file either as a SequenceFile or a DigitalSequenceBlock,
    the latter only if the file loaded to memory would take less than 20% of the available memory.
    """
    available_memory = psutil.virtual_memory().available
    target_size= os.stat(fasta_path).st_size
    print(f"Available memory: {available_memory/1024:.1f} KiB")
    print(f"Database on-disk size: {target_size/1024:.1f} KiB")
    with pyhmmer.easel.SequenceFile(fasta_path, digital=True) as seq_file:
        if target_size < available_memory * 0.2:
            print("Pre-fetching targets into memory")
            targets = seq_file.read_block()
            print(f"Database in-memory size: {(sys.getsizeof(targets) + sum(sys.getsizeof(target) for target in targets))/1024:.1f} KiB")
        else:
            targets = seq_file
    return targets

def run_hmmsearch(hmm, fasta_path, output_dir):
    """Runs pyhmmer hmmsearch with HMMs as query and target (multi)FASTA file and
    returns a JSON with the domains found for the target sequences with structure
    hits_per_domain[accession][target_seq] = [{<domain_data>}]
    for use in dealing with hits on a per domain basis.
    Domain data contains the bitscore, alignment sequence start and end, alignment range and subsequence."""
    with open(hmm, 'rb') as f:
        hmms = list(pyhmmer.plan7.HMMFile(f))

    targets = load_sequence_file(fasta_path)

    hits_per_domain = {}
    per_domain_output = os.path.join(output_dir, "hmmsearch_per_domain_output.json")
    for top_hits in pyhmmer.hmmsearch(hmms, targets, bit_cutoffs="gathering"):
        for hit in top_hits:
            target_seq = hit.name.decode("utf-8")
            for domain in hit.domains.included:
                alignment = domain.alignment
                hmm_name = alignment.hmm_name.decode("utf-8")
                # Accession is the Pfam ID
                accession = alignment.hmm_accession.decode("utf-8") if alignment.hmm_accession.decode("utf-8") != hmm_name else None
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

    with open(per_domain_output, "w", encoding='utf-8') as f:
        json.dump(hits_per_domain, f, indent=4)
    print(f"HmmSearch TopHits results saved per domain - {per_domain_output}")

def main():
    """Main function, initializes this script"""
    args = parse_arguments()
    input_fasta = args.fasta
    input_hmm = args.hmm
    output_dir = args.output_dir
    setup_logging(args.log)

    yielded_sequences = utils.seqrecord_yielder(input_fasta)
    # May stop individualized fastas later,
    # if they end up not being useful in the next steps
    utils.make_dirs_and_write_fasta(yielded_sequences, output_dir)
    run_hmmsearch(input_hmm, input_fasta, output_dir)

if __name__ == '__main__':
    main()
