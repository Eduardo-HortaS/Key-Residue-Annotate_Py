"""
utils.py

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

This script contains all utility functions for SSF-Predict, used in the main script
or by running Snakemake.

"""

import argparse
import logging
import os
import pyhmmer.easel
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
from typing import Iterator
from Bio import SeqIO
from typing import TypeVar, List, Callable


# T = input type (could be str, int, dict, etc.)
# R = return type (could be different from input)
T = TypeVar('T')
R = TypeVar('R')

def get_querynames(fasta: str) -> List[str]:
    """Parses a FASTA file and returns a list of sequence "query names"."""
    sequences = []
    seq = SeqIO.parse(fasta, "fasta")
    for record in seq:
        queryname = record.id.replace("|", "-")
        sequences.append(queryname)
    return sequences

def make_dirs_and_write_fasta(sequences: Iterator[SeqRecord], base_dir: str) -> None:
    """Creates a directory for each sequence and writes the sequence to a fasta file."""
    for record in sequences:
        used_queryname = record.id.replace("|", "-")
        os.makedirs(os.path.join(base_dir, used_queryname), exist_ok=True)
        SeqIO.write(record, os.path.join(base_dir, used_queryname, "sequence.fasta"), "fasta")

def translate_sequence(seq_record: SeqRecord, logger: logging.Logger) -> SeqRecord:
    """Translates a nucleotide sequence using pyHMMER's translation table."""
    # Create digital sequence using pyHMMER
    digital_seq = pyhmmer.easel.TextSequence(
        name=seq_record.id.encode(),
        sequence=str(seq_record.seq).encode()
    ).digitize(pyhmmer.easel.Alphabet.dna())

    # Translate using pyHMMER
    translated = digital_seq.translate()
    logger.info("Translated %s to amino acids.", seq_record.id)
    logger.info("Original sequence: %s", seq_record.seq)
    logger.info("Translated sequence: %s", translated.sequence.decode())

    # Create new SeqRecord with translated sequence
    return SeqRecord(
        Seq(translated.sequence.decode()),
        id=seq_record.id,
        description=seq_record.description
    )

def seqrecord_yielder(fasta: str, is_nucleotide: bool = False, logger: logging.Logger = None) -> Iterator[SeqRecord]:
    """Yields SeqRecord objects from a FASTA file, translating if nucleotide."""
    for record in SeqIO.parse(fasta, "fasta"):
        if is_nucleotide:
            yield translate_sequence(record, logger)
        else:
            yield record

def parallelize(func: Callable[[T], R], inputs: List[T], num_workers: int | None = None) -> List[R]:
    """
    Example usage of T and R:
    - T could be str if inputs is List[str]
    - R could be int if func returns integers

    parallelize(len, ["hello", "world"]) would have:
    - T = str (input strings)
    - R = int (output lengths)
    """

    if num_workers is None:
        num_cores = os.cpu_count()  # Automatically detect the number of CPU cores
        num_workers = max(1, num_cores - 1)  # Use num_cores - 1 to leave one core free

    with Pool(processes=num_workers) as pool:
        results = pool.map(func, inputs)
    return results

def parse_arguments():
    """
    Parse command line arguments for running hmmscan on a fasta file.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="""Runs hmmscan on a fasta file using an hmm file \
        and outputs the results to a specified output dir, also generates domains.json for a sequence's domains, \
        in the output dir's subdirectory relative to each sequence \
        and separates sequences in a multifasta into single fasta files \
        and stores each in the aforementioned subdirectory.""")
    parser.add_argument("-iF", "--input-fasta", help="Path to fasta file", required=True, type=str)
    parser.add_argument("-iH", "--input-hmm", help="Path to hmm file", required=True, type=str)
    parser.add_argument("-o", "--output", help="Output dir path", required=True, type=str)
    parser.add_argument("-l", "--log", help="Log path", \
        required=False, type=str, default="logs/hmmscan.log")
    return parser.parse_args()

def setup_logging(log_file: str ) -> None:
    """Configure logging settings."""
    log_dir = os.path.dirname(log_file)
    os.makedirs(log_dir, exist_ok=True)

    logging.basicConfig(level=logging.INFO, filename=log_file, filemode='w',
                        format='%(asctime)s - %(levelname)s - %(message)s')
