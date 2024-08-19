"""
executor.py

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

This script is intended to be the main executor of the pipeline,
running all scripts in the proper order when called with the necessary arguments.

"""

import os
import json
import glob
import subprocess
import logging
import argparse
import sys
import traceback
from joblib import Parallel, delayed

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run the pipeline")
    parser.add_argument('-iF', '--fasta', required=True, help='Input fasta file')
    parser.add_argument('-iH', '--hmm', required=True, help='Input hmm file')
    parser.add_argument('-r', '--resource_dir', required=True, help='Resource directory')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-t', '--threads', required=False, type=int, default=1, help='Number of threads')
    parser.add_argument('-p', '--python', required=True, help='Path to the Python executable')
    parser.add_argument("-n", "--nucleotide", required=False, default=False, action="store_true", help="Flag to indicate nucleotide instead of default protein sequences")
    parser.add_argument("-e", "--eco-codes", required=False, type=str, default=[], nargs="*", help="List of ECO codes to filter annotations")
    parser.add_argument('-l', '--log', required=True, type=str, help='Log path')
    return parser.parse_args()

def configure_logging(log_path):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename=log_path)
    logger = logging.getLogger()
    return logger

def run_command(command, logger):
    try:
        subprocess.run(command, shell=True, check=True)
    except Exception as e:
        error_info = traceback.format_exc()
        logger.error("Command failed: %s \n Error: %s, Info: %s", command, e, error_info)
        sys.exit(1)

def main():
    """Main function for running the pipeline."""
    args = parse_arguments()
    input_fasta = args.fasta
    input_hmm = args.hmm
    resource_dir = args.resource_dir
    output_dir = args.output_dir
    eco_codes = args.eco_codes
    threads = args.threads
    python_executable = args.python
    logger = configure_logging(args.log)

    # run_hmmsearch.py
    per_dom_json = os.path.join(output_dir, "hmmsearch_per_domain.json")
    if os.path.exists(per_dom_json):
        logger.info("Output for hmmsearch step %s already exists. Skipping.", per_dom_json)
    else:
        run_hmmsearch_call = f"{python_executable} run_hmmsearch.py -iF {input_fasta} -iH {input_hmm} -o {output_dir} -l {args.log}"
        run_command(run_hmmsearch_call, logger)

    # prepare_fasta_per_domain.py
    prepare_fasta_done = os.path.join(output_dir, "prepare_fasta_per_domain.done")
    if os.path.exists(prepare_fasta_done):
        logger.info("PREPARE_FASTA_PER_DOMAIN.PY --- Skipping, output already exists")
    else:
        with open(per_dom_json, "r", encoding="utf-8") as f:
            hits_per_domain = json.load(f)

        prepare_fasta_tasks = [
            f"{python_executable} prepare_fasta_per_domain.py -iJ {per_dom_json} -iD {dom_accession} -r {resource_dir} -o {output_dir} -l {args.log}"
            for dom_accession in hits_per_domain
        ]
        Parallel(n_jobs=threads)(delayed(run_command)(task, logger) for task in prepare_fasta_tasks)
        with open(prepare_fasta_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("PREPARE_FASTA_PER_DOMAIN.PY --- Executed")

    # run_hmmalign.py
    run_hmmalign_done = os.path.join(output_dir, "run_hmmalign.done")
    if os.path.exists(run_hmmalign_done):
        logger.info("RUN_HMMALIGN.PY --- Skipping, output already exists")
    else:
        run_hmmalign_tasks = []
        for subdir in os.listdir(output_dir):
            subdir_path = os.path.join(output_dir, subdir)
            if os.path.isdir(subdir_path) and subdir.startswith("PF"):
                domain_info = os.path.join(subdir_path, "domain_info.json")
                if os.path.isfile(domain_info):
                    run_hmmalign_tasks.append(f"{python_executable} run_hmmalign.py -iDI {domain_info} -l {args.log}")
        Parallel(n_jobs=threads)(delayed(run_command)(task, logger) for task in run_hmmalign_tasks)
        with open(run_hmmalign_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("RUN_HMMALIGN.PY --- Executed.")

    # transfer_annotations.py
    transfer_annotations_done = os.path.join(output_dir, "transfer_annotations.done")
    if os.path.exists(transfer_annotations_done):
        logger.info("TRANSFER_ANNOTATIONS.PY --- Skipping, output already exists")
    else:
        transfer_annotations_tasks = []
        for subdir in os.listdir(output_dir):
            subdir_path = os.path.join(output_dir, subdir)
            if os.path.isdir(subdir_path) and subdir.startswith("PF"):
                dom_aligns = [dom_align for dom_align in glob.glob(os.path.join(subdir_path, "PF*_hmmalign.sth")) if os.path.isfile(dom_align)]
                for dom_align in dom_aligns:
                    transfer_annotations_tasks.append(f"{python_executable} transfer_annotations.py -iA {dom_align} -r {resource_dir} -o {output_dir} --eco-codes {' '.join(eco_codes)} -l {args.log}")
        logger.info("Transfer annotations tasks: %s", transfer_annotations_tasks)
        Parallel(n_jobs=threads)(delayed(run_command)(task, logger) for task in transfer_annotations_tasks)
        with open(transfer_annotations_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("TRANSFER_ANNOTATIONS.PY --- Executed")

    # merge_sequences.py
    merge_sequences_done = os.path.join(output_dir, "merge_sequences.done")
    if os.path.exists(merge_sequences_done):
        logger.info("MERGE_SEQUENCES --- Skipping, output already exists")
    else:
        merge_sequences_call = f"{python_executable} merge_sequences.py -o {output_dir} -l {args.log}"
        run_command(merge_sequences_call, logger)
        with open(merge_sequences_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("MERGE_SEQUENCES --- Executed")

    logger.info("Pipeline finished successfully")



if __name__ == "__main__":
    main()
