"""
seq_and_batch_prep.py

Copyright 2025 Eduardo Horta Santos <GitHub: Eduardo-HortaS>

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
a FASTA file, a HMM targets file and the path to the output directory. Optionally,
a flag to indicate nucleotide sequences and a log file path.

1 - Splits multifasta file into individual sequence fastas,
each inside its respective {sequence_name} subdir.

2 - B

# Moved these general functions to utils.py, for better modularity:
    2 - Parses the FASTA file and stores each sequence's queryname (ex.: TTHY_XENLA) in a list,
    and generates a SeqRecord object to be used sequentially and *only once*.

    3 - For each sequence with at least 1 domain hit, creates a directory
    and writes the sequence to a 'sequence.fasta' file in the directory.

"""

import os
import json
import argparse
import logging
from typing import Iterator
from Bio.SeqRecord import SeqRecord
from utils import get_logger, seqrecord_yielder, make_dirs_and_write_fasta

def parse_arguments():
    """Parse command-line arguments for sequence and batch preparation"""
    parser = argparse.ArgumentParser(description="Prepare sequences and batches for pipeline processing")
    parser.add_argument("-iF", "--fasta", help="Path to input FASTA file",
                    required=True, type=str)
    parser.add_argument("-o", "--output-dir", help="Output directory path",
                    required=True, type=str)
    parser.add_argument("-b", "--batch-size", help="Maximum sequences per batch",
                    type=int, default=2000)
    parser.add_argument("-n", "--nucleotide",
                    help="Flag to indicate nucleotide instead of default protein sequences",
                    action="store_true")
    parser.add_argument("-l", "--log",
                    help="Log path",
                    required=False, type=str, default="logs/seq_and_batch_prep.log")
    return parser.parse_args()


def create_sequence_batches_json(
    sequences: Iterator[SeqRecord],
    batch_size: int,
    output_dir: str,
    logger: logging.Logger) -> None:
    """Creates JSON file mapping batch numbers to sequence IDs.

    Args:
        sequences: Iterator of SeqRecord objects
        batch_size: Maximum sequences per batch
        output_dir: Output directory path
        logger: Logger instance for tracking execution
    """
    all_sequences = {}
    current_batch = []
    batch_num = 1

    for record in sequences:
        seq_id = record.id
        current_batch.append(seq_id)

        if len(current_batch) >= batch_size:
            all_sequences[f"batch_{batch_num}"] = current_batch
            current_batch = []
            batch_num += 1

    # Handle any remaining sequences
    if current_batch:
        all_sequences[f"batch_{batch_num}"] = current_batch

    output_path = os.path.join(output_dir, "all_sequences.json")
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(all_sequences, f, indent=4)

    logger.info("SEQ_BATCH_PREP --- Created all_sequences.json with %d batches", batch_num)

def main():
    """Main function to prepare sequences and batches"""
    args = parse_arguments()
    logger, _ = get_logger(args.log, scope="main")

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # First pass - create sequence batches JSON
    sequences = seqrecord_yielder(args.fasta, args.nucleotide, logger)
    create_sequence_batches_json(sequences, args.batch_size, args.output_dir, logger)

    # Second pass - create individual sequence directories/files
    sequences = seqrecord_yielder(args.fasta, args.nucleotide, logger)
    make_dirs_and_write_fasta(sequences, args.output_dir)

    logger.info("SEQ_BATCH_PREP --- Sequence and batch preparation completed")

if __name__ == "__main__":
    main()