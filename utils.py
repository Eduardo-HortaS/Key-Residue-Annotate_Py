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
from multiprocessing import Pool
from Bio import SeqIO

def get_querynames(fasta):
    """Parses a FASTA file and returns a list of sequence "query names"."""
    sequences = []
    seq = SeqIO.parse(fasta, "fasta")
    for record in seq:
        queryname = record.id.replace("|", "-")
        sequences.append(queryname)
    return sequences

def make_dirs_and_write_fasta(sequences, base_dir):
    """Creates a directory for each sequence and writes the sequence to a fasta file."""
    for record in sequences:
        queryname = record.id.replace("|", "-")
        os.makedirs(os.path.join(base_dir, queryname), exist_ok=True)
        SeqIO.write(record, os.path.join(base_dir, queryname, "sequence.fasta"), "fasta")

def seqrecord_yielder(fasta):
    """Yields SeqRecord objects from a FASTA file."""
    for record in SeqIO.parse(fasta, "fasta"):
        yield record

def parallelize(func, inputs, num_workers=None):
    if num_workers is None:
        num_cores = os.cpu_count()  # Automatically detect the number of CPU cores
        num_workers = max(1, num_cores - 1)  # Use num_cores - 1 to leave one core free

    with Pool(processes=num_workers) as pool:
        results = pool.map(func, inputs)
    return results

# TODO: Change it to work with all necessary arguments for the whole pipeline
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

# TODO: Change it to work separately or only once for the whole pipeline.
def setup_logging(log_file):
    """Configure logging settings."""
    log_dir = os.path.dirname(log_file)
    os.makedirs(log_dir, exist_ok=True)

    logging.basicConfig(level=logging.INFO, filename=log_file, filemode='w',
                        format='%(asctime)s - %(levelname)s - %(message)s')
