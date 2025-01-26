"""
executor.py

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
from configparser import ConfigParser
from joblib import Parallel, delayed

def load_config(config_file=None):
    """Load configuration from INI file"""
    config = ConfigParser()
    if config_file and os.path.exists(config_file):
        config.read(config_file)
        return {
            "fasta": config.get("Inputs", "fasta", fallback=None),
            "hmm": config.get("Inputs", "hmm", fallback=None),
            "iprscan_path": config.get("Paths", "iprscan_path", fallback=None),
            "output_format_iprscan": config.get("Parameters", "output_format_iprscan", fallback="TSV, XML, GFF3"),
            "databases": config.get("Parameters", "databases", fallback=""),
            "resource_dir": config.get("Paths", "resource_dir", fallback=None),
            "output_dir": config.get("Paths", "output_dir", fallback=None),
            "threads": config.getint("Parameters", "threads", fallback=1),
            "python": config.get("Paths", "python", fallback=sys.executable),
            "nucleotide": config.getboolean("Parameters", "nucleotide", fallback=False),
            "eco_codes": config.get("Parameters", "eco_codes", fallback="").split(),
            "log": config.get("Paths", "log", fallback="logs/executor.log")
        }
    return {}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run the pipeline")
    parser.add_argument("-c", "--config", help="Path to config.ini file", type=str)
    parser.add_argument("-iF", "--fasta", help="Input fasta file")
    parser.add_argument("-iH", "--hmm", help="Input hmm file")
    parser.add_argument("-iP", "--iprscan-path", type=str, help="Path to interproscan.sh")
    parser.add_argument("-of", "--output-format-iprscan", type=str, help="Output format for interproscan")
    parser.add_argument("-dB", "--databases", help="Optional: Comma-separated databases to limit InterProScan search, if the single purpose is using GO terms inside the pipeline, pass the string panther,gene3d,smart,pfam,superfamily")
    parser.add_argument("-r", "--resource-dir", help="Resource directory")
    parser.add_argument("-o", "--output-dir", help="Output directory")
    parser.add_argument("-t", "--threads", type=int, help="Number of threads")
    parser.add_argument("-p", "--python", help="Path to the Python executable")
    parser.add_argument("-n", "--nucleotide", action="store_true", help="Flag for nucleotide sequences")
    parser.add_argument("-e", "--eco-codes", nargs="*", help="Space-separated ECO codes")
    parser.add_argument("-l", "--log", type=str, help="Log path")

    args = parser.parse_args()

    # Load config file first
    config = load_config(args.config)

    # Override with any command line arguments that were specified
    cmd_args = {k: v for k, v in vars(args).items() if v is not None}
    config.update(cmd_args)

    # Validate required parameters
    required = ["fasta", "hmm", "iprscan_path", "resource_dir", "output_dir"]
    missing = [param for param in required if param not in config or not config[param]]
    if missing:
        parser.error(f"Missing required parameters: {', '.join(missing)}")

    return argparse.Namespace(**config)

def configure_logging(log_path):
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", filename=log_path)
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
    iprscan_sh_path = args.iprscan_path
    output_format_iprscan = args.output_format_iprscan
    databases = args.databases
    resource_dir = args.resource_dir
    output_dir = args.output_dir
    eco_codes = args.eco_codes
    threads = args.threads
    nucleotide = args.nucleotide
    python_executable = args.python
    logger = configure_logging(args.log)
    log = args.log

    # run_hmmsearch.py
    per_dom_json = os.path.join(output_dir, "hmmsearch_per_domain.json")
    if os.path.exists(per_dom_json):
        logger.info("Output for hmmsearch step %s already exists. Skipping.", per_dom_json)
    else:
        run_hmmsearch_call = (f"{python_executable} run_hmmsearch.py "
                             f"-iF {input_fasta} -iH {input_hmm} "
                             f"-o {output_dir} -l {log}"
                             f"{' -n' if nucleotide else ''}")
        run_command(run_hmmsearch_call, logger)

    # run_iprscan.py
    run_iprscan_done = os.path.join(output_dir, "run_iprscan.done")
    if os.path.exists(run_iprscan_done):
        logger.info("RUN_IPRSCAN.PY --- Skipping, output already exists")
    else:
        run_iprscan_tasks = []
        for subdir in os.listdir(output_dir):
            subdir_path = os.path.join(output_dir, subdir)
            sequence_fasta = os.path.join(subdir_path, "sequence.fasta")
            if os.path.isdir(subdir_path) and not subdir.startswith("PF") and os.path.isfile(sequence_fasta):
                output_base_file = os.path.join(subdir_path, "iprscan")
                cmd = f"{python_executable} run_iprscan.py -iP {iprscan_sh_path} -iF {sequence_fasta} -oB {output_base_file} -oF {output_format_iprscan}"
                if databases:  # Only add -dB if databases is not empty
                    cmd += f" -dB {databases}"
                cmd += f" -l {log}"
                run_iprscan_tasks.append(cmd)
        Parallel(n_jobs=threads)(delayed(run_command)(task, logger) for task in run_iprscan_tasks)
        with open(run_iprscan_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("RUN_IPRSCAN.PY --- Executed")

    # prepare_fasta_per_domain.py
    prepare_fasta_done = os.path.join(output_dir, "prepare_fasta_per_domain.done")
    if os.path.exists(prepare_fasta_done):
        logger.info("PREPARE_FASTA_PER_DOMAIN.PY --- Skipping, output already exists")
    else:
        with open(per_dom_json, "r", encoding="utf-8") as f:
            hits_per_domain = json.load(f)

        prepare_fasta_tasks = [
            f"{python_executable} prepare_fasta_per_domain.py -iJ {per_dom_json} -iD {dom_accession} -r {resource_dir} -o {output_dir} -l {log}"
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
                    run_hmmalign_tasks.append(f"{python_executable} run_hmmalign.py -iDI {domain_info} -l {log}")
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
                    transfer_annotations_tasks.append(f"{python_executable} transfer_annotations.py -iA {dom_align} -r {resource_dir} -o {output_dir} --eco-codes {' '.join(eco_codes)} -l {log}")
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
        merge_sequences_call = f"{python_executable} merge_sequences.py -o {output_dir} -l {log}"
        run_command(merge_sequences_call, logger)
        with open(merge_sequences_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("MERGE_SEQUENCES --- Executed")

    # make_view_jsons.py
    make_view_jsons_done = os.path.join(output_dir, "make_view_jsons.done")
    if os.path.exists(make_view_jsons_done):
        logger.info("MAKE_VIEW_JSONS.PY --- Skipping, output already exists")
    else:
        make_view_jsons_call = f"{python_executable} make_view_jsons.py -o {output_dir} -l {log}"
        run_command(make_view_jsons_call, logger)
        with open(make_view_jsons_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("MAKE_VIEW_JSONS.PY --- Executed")

    logger.info("Pipeline finished successfully")


if __name__ == "__main__":
    main()
