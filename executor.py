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
import argparse
import subprocess
import glob
from joblib import Parallel, delayed

def parse_arguments():
    """
    Parse command-line arguments for running all
    scripts in the pipeline.
    Neccessary arguments: proteome fasta, HMM database,
    output directory and resource directory.
    Optional: number of threads to use, default is 1.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """

    parser = argparse.ArgumentParser(description=
    'Runs hmmsearch using a sequence database against target HMMs \
    aiming to generate a hmmsearch_per_domain.json file.')
    parser.add_argument("-iF", "--fasta", help="Path to fasta file", required=True, type=str)
    parser.add_argument("-iH", "--hmm", help="Path to target HMMs database file", required=True, type=str)
    parser.add_argument("-r", "--resource-dir", help="Resource dir path", required=True, type=str)
    parser.add_argument("-o", "--output-dir", help="Output dir path", required=True, type=str)
    parser.add_argument("-t", "--threads", help="Number of threads to use", required=False, type=int, default=1)
    return parser.parse_args()

def main():
    """Main function for running the pipeline."""
    args = parse_arguments()
    input_fasta = args.fasta
    input_hmm = args.hmm
    resource_dir = args.resource_dir
    output_dir = args.output_dir
    threads = args.threads

    # run_hmmsearch.py
    run_hmmsearch_call = f"python3 run_hmmsearch.py -iF {input_fasta} -iH {input_hmm} -o {output_dir}"
    subprocess.run(run_hmmsearch_call, shell=True, check=True)

    # prepare_fasta_per_domain.py
    per_dom_json = os.path.join(output_dir, "hmmsearch_per_domain.json")
    with open(per_dom_json, "r", encoding="utf-8") as f:
        hits_per_domain = json.load(f)

    prepare_fasta_tasks = [
        f"python3 prepare_fasta_per_domain.py -iJ {per_dom_json} -iD {dom_accession} -r {resource_dir} -o {output_dir}"
        for dom_accession in hits_per_domain
    ]
    Parallel(n_jobs=threads)(delayed(subprocess.run)(task, shell=True, check=True) for task in prepare_fasta_tasks)

    # run_hmmalign.py
    run_hmmalign_tasks = []
    for subdir in os.listdir(output_dir):
        subdir_path = os.path.join(output_dir, subdir)
        if os.path.isdir(subdir_path) and subdir.startswith("PF"):
            domain_info = os.path.join(subdir_path, "domain_info.json")
            if os.path.isfile(domain_info):
                run_hmmalign_tasks.append(f"python3 run_hmmalign.py -iDI {domain_info}")
    Parallel(n_jobs=threads)(delayed(subprocess.run)(task, shell=True, check=True) for task in run_hmmalign_tasks)

    # transfer_annotations.py
    transfer_annotations_tasks = []
    for subdir in os.listdir(output_dir):
        subdir_path = os.path.join(output_dir, subdir)
        if os.path.isdir(subdir_path) and subdir.startswith("PF"):
            dom_aligns = [dom_align for dom_align in glob.glob(os.path.join(subdir_path, "*_hmmalign.sth")) if os.path.isfile(dom_align)]
            for dom_align in dom_aligns:
                transfer_annotations_tasks.append(f"python3 transfer_annotations.py -iA {dom_align} -r {resource_dir} -o {output_dir}")
    Parallel(n_jobs=threads)(delayed(subprocess.run)(task, shell=True, check=True) for task in transfer_annotations_tasks)

    # merge_sequences.py
    merge_sequences_call = f"python3 merge_sequences.py -o {output_dir}"
    subprocess.run(merge_sequences_call, shell=True, check=True)
    print("Pipeline finished successfully.")

if __name__ == "__main__":
    main()
