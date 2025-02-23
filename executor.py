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
import psutil
import subprocess
import argparse
import sys
import logging
from configparser import ConfigParser
from joblib import Parallel, delayed
from utils import get_logger

def load_config(config_file=None):
    """Load configuration from INI file"""
    config = ConfigParser()
    if config_file and os.path.exists(config_file):
        config.read(config_file)
        return {
            "fasta": config.get("Inputs", "fasta",
            fallback=None),
            "hmm": config.get("Inputs", "hmm",
            fallback=None),
            "iprscan_path": config.get("Paths", "iprscan_path",
            fallback=None),
            "resource_dir": config.get("Paths", "resource_dir",
            fallback=None),
            "output_dir": config.get("Paths", "output_dir",
            fallback=None),
            "python": config.get("Paths", "python",
            fallback=sys.executable),
            "log": config.get("Paths", "log",
            fallback=None),
            "output_format_iprscan": config.get("Parameters", "output_format_iprscan",
            fallback="TSV, XML, GFF3"),
            "cpu_cores_iprscan": config.getint("Parameters", "cpu_cores_iprscan",
            fallback=8),
            "number_jobs_iprscan": config.getint("Parameters", "number_jobs_iprscan",
            fallback=1),
            "seq_batch_size_iprscan": config.getint("Parameters", "seq_batch_size_iprscan",
            fallback=2000),
            "analyses_iprscan": config.get("Parameters", "analyses_iprscan",
            fallback=""),
            "enable_precalc_iprscan": config.getboolean("Parameters", "enable_precalc_iprscan",
            fallback=False),
            "disable_res_iprscan": config.getboolean("Parameters", "disable_res_iprscan",
            fallback=False),
            "threads": config.getint("Parameters", "threads",
            fallback=2),
            "total_memory": config.getint("Parameters", "total_memory",
            fallback=None),
            "nucleotide": config.getboolean("Parameters", "nucleotide",
            fallback=False),
            "eco_codes": config.get("Parameters", "eco_codes",
            fallback="").split(),
        }
    return {}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run the pipeline")
    parser.add_argument("-c", "--config", help="Path to config.ini file", type=str, required=False, default=None)
    # These are required, but will be checked for later, after merging params with config file
    parser.add_argument("-f", "--fasta", help="Input fasta file",
                        required=False, default=None)
    parser.add_argument("-iH", "--hmm", help="Input hmm file",
                        required=False, default=None)
    parser.add_argument("-iPr", "--iprscan-path", type=str, help="Path to InterProScan.sh",
                        required=False, default=None)
    parser.add_argument("-r", "--resource-dir", help="Resource directory",
                        required=False, default=None)
    parser.add_argument("-o", "--output-dir", help="Output directory",
                        required=False, default=None)
    # Optional parameters
    parser.add_argument("-iOf", "--output-format-iprscan", type=str,
                        help="Output format for InterProScan",
                        required=False, default="TSV, XML, GFF3")
    parser.add_argument("-iCc", "--cpu-cores-iprscan", type=int,
                        help="Number of CPU cores to use per job",
                        required=False, default=8)
    parser.add_argument("-iNj", "--number-jobs-iprscan",
                        type=int, help="Number of parallel jobs to run InterProScan",
                        required=False, default=1)
    parser.add_argument("-iBs", "--seq-batch-size-iprscan",
                        type=int, help="Number of sequences per batch in InterProScan",
                        required=False, default=2000)
    parser.add_argument("-iA", "--analyses-iprscan",
                        help="Optional: Comma-separated analyses to limit InterProScan run. \
                        If you just want what's used in the pipeline (GO terms), pass the string: \
                        panther,gene3d,smart,pfam,superfamily",
                        required=False, default="")
    parser.add_argument("-iPc", "--enable-precalc-iprscan",
                        action="store_true",
                        help="Flag to enable pre-calculated match lookup in InterProScan. \
                        Keep this off unless running with proteins in UniProtKB.",
                        required=False)
    parser.add_argument("-iDr", "--disable-res-iprscan",
                        action="store_true",
                        help="Flag to disable residue-level annotations in InterProScan outputs, \
                        may increase performance, but possible miss some cross-match annotations.",
                        required=False)
    parser.add_argument("-t", "--threads", type=int, help="Number of threads",
                        required=False, default=2)
    parser.add_argument("-m", "--total_memory", type=int,
                        help="Total memory available in GB",
                        required=False, default=None)
    parser.add_argument("-n", "--nucleotide", action="store_true",
                        help="Flag for nucleotide sequences",
                        required=False)
    parser.add_argument("-e", "--eco-codes", nargs="*",
                        help="Space-separated ECO codes",
                        required=False, default="")
    parser.add_argument("-p", "--python",
                        help="Path to the Python executable",
                        required=False, default=sys.executable)
    parser.add_argument("-l", "--log", type=str, help="Log path",
                        required=False, default=None)
    args = parser.parse_args()

    # Load config file first
    config = load_config(args.config)

    # Override with any command line arguments that were specified
    defaults = {}
    for action in parser._actions:
        if action.default is not argparse.SUPPRESS:
            defaults[action.dest] = action.default

    cmd_args = {}
    for k, v in vars(args).items():
        if v is not None and v != defaults.get(k, None):
            cmd_args[k] = v

    config.update(cmd_args)

    # Validate required parameters
    required = ["fasta", "hmm", "iprscan_path", "resource_dir", "output_dir"]
    missing = [param for param in required if param not in config or not config[param]]
    if missing:
        parser.error(f"Missing required parameters: {', '.join(missing)}")

    return argparse.Namespace(**config)

def run_command(command: list, logger: logging.Logger):
    """Run command, capture output and log any errors. Exit if command fails.

    Args:
        command: List comprised of the command and its arguments
        logger: Logger instance
        is_parallel: If True, run command in parallel mode
    """
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            logger.error(
                "EXECUTOR --- RUN_CMD --- Command failed with return code %d: %s\n"
                "STDOUT:\n%s\n"
                "STDERR:\n%s",
                result.returncode,
                command,
                result.stdout,
                result.stderr
            )
            sys.exit(1)
    except subprocess.CalledProcessError as e:
        logger.error("Command failed with return code %d: %s", e.returncode, e.cmd)
        logger.error("STDOUT:\n%s", e.stdout)
        logger.error("STDERR:\n%s", e.stderr)
        sys.exit(1)

def get_seqs_and_count(json_file: str) -> tuple[list[str], int]:
    """Get list of all sequences and total count from all_sequences.json file.

    Args:
        json_file: Path to all_sequences.json

    Returns:
        tuple[list[str], int]: List of all sequence IDs and total count
    """
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)

    # Flatten all batches into a single list
    all_sequences = []
    for batch in data.values():
        all_sequences.extend(batch)

    return all_sequences, len(all_sequences)

def create_sequence_batches(sequences: list[str], batch_size: int) -> list[list[str]]:
    """Create batches of sequences of specified size.

    Args:
        sequences: List of sequence IDs
        batch_size: Number of sequences per batch

    Returns:
        list[list[str]]: List of sequence batches
    """
    return [
        sequences[i:i + batch_size]
        for i in range(0, len(sequences), batch_size)
    ]

def validate_iprscan_resources(
    cpu_cores_iprscan: int,
    seq_batch_size: int,
    logger: logging.Logger,
    total_memory: int,
) -> bool:
    """Validate if resources allow InterProScan to run successfully.
    Note that increasing the number of sequences per batch or cores
    too much may lead to running out of memory. Refer to
    https://interproscan-docs.readthedocs.io/en/latest/ImprovingPerformance.html

    Args:
        cpu_cores_iprscan: CPU cores allocated per InterProScan job
        seq_batch_size: Number of sequences per batch
        logger: Logger instance for output
        total_memory: Available RAM in GB. Auto-detected if None.

    Returns:
        bool: True if resources are sufficient, False otherwise
    """
    # Recommended Constants from InterProScan Docs - Conservative
    memory_per_core = 0.5  # 8GB/16 cores = 0.5GB/core
    system_reserve_gb = 2  # Reserve for OS/other processes
    min_cores = 3  # Minimum cores needed (1 for main process + 2 for worker)
    # Maximum recommended seq batch size, increase at your own risk (+ memory req.)
    max_rec_seq_batch_size = 8000

    # Auto-detect memory if not provided
    if total_memory is None:
        total_memory = psutil.virtual_memory().available / (1024**3)  # Convert bytes to GB
    usable_memory = max(0, total_memory - system_reserve_gb)

    # Validate CPU cores
    if cpu_cores_iprscan < min_cores:
        print("CPU cores", cpu_cores_iprscan)
        print("Min cores", min_cores)
        logger.error(
            "EXECUTOR --- VAL_IPRSCAN_RESOURCES --- Insufficient CPU cores allocated (%d). \
            InterProScan needs at least %d cores to run comfortably",
            cpu_cores_iprscan, min_cores
        )
        return False

    # Validate available memory
    required_memory = cpu_cores_iprscan * memory_per_core
    if required_memory > usable_memory:
        logger.error(
            "EXECUTOR --- VAL_IPRSCAN_RESOURCES --- \
            Insufficient memory for requested cores. Need %.1fGB for %d cores, "
            "but only %.1fGB available after %d system reserve",
            required_memory, cpu_cores_iprscan, usable_memory, system_reserve_gb
        )
        return False

    # Validate batch size
    if seq_batch_size > max_rec_seq_batch_size:
        logger.warning(
            "EXECUTOR --- VAL_IPRSCAN_RESOURCES --- \
            Batch size (%d) exceeds recommended maximum (%d). This may cause memory issues, \
            but if enough cores and memory are available, this may not occur.",
            seq_batch_size, max_rec_seq_batch_size
        )


    logger.info(
        "EXECUTOR --- VAL_IPRSCAN_RESOURCES --- Resource validation passed: %d cores - %.1fGB usable memory (%.1fGB memory required), %d sequences per batch",
        cpu_cores_iprscan, usable_memory, required_memory, seq_batch_size
    )
    return True

def main():
    """Main function for running the pipeline."""
    args = parse_arguments()
    input_fasta = args.fasta
    input_hmm = args.hmm
    iprscan_sh_path = args.iprscan_path
    output_format_iprscan = args.output_format_iprscan
    cpu_cores_iprscan = args.cpu_cores_iprscan
    number_jobs_iprscan = args.number_jobs_iprscan
    seq_batch_size_iprscan = args.seq_batch_size_iprscan
    analyses_iprscan = args.analyses_iprscan
    enable_precalc_iprscan = args.enable_precalc_iprscan
    disable_res_iprscan = args.disable_res_iprscan
    resource_dir = args.resource_dir
    output_dir = args.output_dir
    eco_codes = args.eco_codes
    threads = args.threads
    total_memory = args.total_memory
    nucleotide = args.nucleotide
    python_executable = args.python
    logger, timestamped_log = get_logger(args.log)
    all_sequences_json = os.path.join(output_dir, "all_sequences.json")

    # run_hmmsearch.py
    per_dom_json = os.path.join(output_dir, "hmmsearch_per_domain.json")
    if os.path.exists(per_dom_json):
        logger.info("EXECUTOR --- RUN_HMMSEARCH.PY --- \
        Output for hmmsearch step %s already exists. Skipping.", per_dom_json)
    else:
        run_hmmsearch_call = [
            python_executable,
            "run_hmmsearch.py",
            "-iF", input_fasta,
            "-iH", input_hmm,
            "-o", output_dir,
            "-l", timestamped_log,
        ]
        if nucleotide:
            run_hmmsearch_call.append("-n")
        run_command(run_hmmsearch_call, logger)
        logger.info("EXECUTOR --- RUN_HMMSEARCH.PY --- Executed.")

    # seq_and_batch_prep.py
    if os.path.exists(all_sequences_json):
        logger.info("EXECUTOR --- SEQ_AND_BATCH_PREP.PY --- Output already exists %s. Skipping.", all_sequences_json)
    else:
        seq_batch_prep_call = [
            python_executable,
            "seq_and_batch_prep.py",
            "-iF", input_fasta,
            "-o", output_dir,
            "-b", str(seq_batch_size_iprscan),
            "-l", timestamped_log,
        ]
        if nucleotide:
            seq_batch_prep_call.append("-n")
        run_command(seq_batch_prep_call, logger)
        logger.info("EXECUTOR --- SEQ_AND_BATCH_PREP.PY --- Executed.")

    # Count sequences from hmmsearch_sequences.json
    list_of_sequences, seq_count = get_seqs_and_count(all_sequences_json)
    sequence_batches = create_sequence_batches(list_of_sequences, seq_batch_size_iprscan)

    if seq_count < seq_batch_size_iprscan:
        can_run = validate_iprscan_resources(cpu_cores_iprscan, seq_count, logger, total_memory)
    else:
        can_run = validate_iprscan_resources(cpu_cores_iprscan, seq_batch_size_iprscan, logger, total_memory)

    if not can_run:
        logger.warning("EXECUTOR --- VAL_IPRSCAN_RESOURCES --- No resources available for InterProScan. Exiting pipeline.")
        sys.exit(1)

    # run_iprscan.py
    run_iprscan_done = os.path.join(output_dir, "run_iprscan.done")
    if os.path.exists(run_iprscan_done):
        logger.info("EXECUTOR --- RUN_IPRSCAN.PY --- Skipping, output already exists")
    else:
        run_iprscan_tasks = []
        for batch_idx, sequence_batch in enumerate(sequence_batches, 1):
            batch_seq_str = ",".join(sequence_batch)
            batch_seq_str = batch_seq_str.replace("|", "-")  # Avoid shell parsing issues with pipes

            # General arguments
            cmd = [
                python_executable,
                "run_iprscan.py",
                "-iPr", iprscan_sh_path,
                "-iOf", output_format_iprscan,
                "-iCc", str(cpu_cores_iprscan),
                "-l", timestamped_log,
            ]

            if analyses_iprscan:
                cmd.extend(["-iA", analyses_iprscan])
            if enable_precalc_iprscan:
                cmd.append("-iDpc")
            if disable_res_iprscan:
                cmd.append("-iDr")

            # Mode-specific arguments for 'batch'
            cmd.extend([
                "batch",
                "-sB", batch_seq_str,
                "-sBi", str(batch_idx),
                "-sPd", output_dir,
            ])
            run_iprscan_tasks.append(cmd)

        Parallel(n_jobs=number_jobs_iprscan)(
            delayed(run_command)(task, logger)
            for task in run_iprscan_tasks
        )
        with open(run_iprscan_done, "w", encoding="utf-8") as f:
            # Write the names of all sequences in the batch, open each batch *done file to know which sequences should not be in a new batch.
            f.write("")
        logger.info("EXECUTOR --- RUN_IPRSCAN.PY --- Executed.")

    # prepare_fasta_per_domain.py
    prepare_fasta_done = os.path.join(output_dir, "prepare_fasta_per_domain.done")
    if os.path.exists(prepare_fasta_done):
        logger.info("EXECUTOR --- PREPARE_FASTA_PER_DOMAIN.PY --- Skipping, output already exists")
    else:
        with open(per_dom_json, "r", encoding="utf-8") as f:
            hits_per_domain = json.load(f)

        prepare_fasta_tasks = [
            [
                python_executable,
                "prepare_fasta_per_domain.py",
                "-iJ", per_dom_json,
                "-iD", dom_accession,
                "-r", resource_dir,
                "-o", output_dir,
                "-l", timestamped_log
            ]
            for dom_accession in hits_per_domain
        ]

        Parallel(n_jobs=threads)(
            delayed(run_command)(task, logger)
            for task in prepare_fasta_tasks
        )
        with open(prepare_fasta_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("EXECUTOR --- PREPARE_FASTA_PER_DOMAIN.PY --- Executed.")

    # run_hmmalign.py
    run_hmmalign_done = os.path.join(output_dir, "run_hmmalign.done")
    if os.path.exists(run_hmmalign_done):
        logger.info("EXECUTOR --- RUN_HMMALIGN.PY --- Skipping, output already exists")
    else:
        run_hmmalign_tasks = []
        for subdir in os.listdir(output_dir):
            subdir_path = os.path.join(output_dir, subdir)
            if os.path.isdir(subdir_path) and subdir.startswith("PF"):
                domain_info = os.path.join(subdir_path, "domain_info.json")
                if os.path.isfile(domain_info):
                    run_hmmalign_tasks.append([
                        python_executable,
                        "run_hmmalign.py",
                        "-iDI", domain_info,
                        "-d", subdir,
                        "-l", timestamped_log
                    ])
        Parallel(n_jobs=threads)(
            delayed(run_command)(task, logger)
            for task in run_hmmalign_tasks
        )
        with open(run_hmmalign_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("EXECUTOR --- RUN_HMMALIGN.PY --- Executed.")

    # transfer_annotations.py
    transfer_annotations_done = os.path.join(output_dir, "transfer_annotations.done")
    if os.path.exists(transfer_annotations_done):
        logger.info("EXECUTOR --- TRANSFER_ANNOTATIONS.PY --- Skipping, output already exists")
    else:
        transfer_annotations_tasks = []
        for subdir in os.listdir(output_dir):
            subdir_path = os.path.join(output_dir, subdir)
            if os.path.isdir(subdir_path) and subdir.startswith("PF"):
                dom_aligns = [dom_align for dom_align in glob.glob(os.path.join(subdir_path, "PF*_hmmalign.sth")) if os.path.isfile(dom_align)]
                for dom_align in dom_aligns:
                    transfer_annotations_tasks.append([
                        python_executable,
                        "transfer_annotations.py",
                        "-iA", dom_align,
                        "-r", resource_dir,
                        "-d", subdir,
                        "-o", output_dir,
                        "--eco-codes", *eco_codes,
                        "-l", timestamped_log
                    ])
        logger.debug("EXECUTOR --- TRANSFER_ANNOTATIONS.PY --- Number of Transfer annotations tasks: %s", len(transfer_annotations_tasks))
        Parallel(n_jobs=threads)(
            delayed(run_command)(task, logger)
            for task in transfer_annotations_tasks
        )
        with open(transfer_annotations_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("EXECUTOR --- TRANSFER_ANNOTATIONS.PY --- Executed.")

    # merge_reports_in_sequences.py
    merge_reports_in_sequences = os.path.join(output_dir, "merge_sequences.done")
    if os.path.exists(merge_reports_in_sequences):
        logger.info("EXECUTOR --- MERGE_REPORT_SEQUENCES --- Skipping, output already exists")
    else:
        merge_reports_in_sequences_tasks = []
        for subdir in os.listdir(output_dir):
            subdir_path = os.path.join(output_dir, subdir)
            if os.path.isdir(subdir_path) and not subdir.startswith("PF"):
                merge_reports_in_sequences_tasks.append([
                    python_executable,
                    "merge_reports_in_sequences.py",
                    "-s", subdir,
                    "-sd", subdir_path,
                    "-l", timestamped_log
                ])
        Parallel(n_jobs=threads)(
            delayed(run_command)(task, logger)
            for task in merge_reports_in_sequences_tasks
        )
        with open(merge_reports_in_sequences, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("EXECUTOR --- MERGE_REPORT_SEQUENCES --- Executed.")

    # make_view_jsons.py
    make_view_jsons_done = os.path.join(output_dir, "make_view_jsons.done")
    if os.path.exists(make_view_jsons_done):
        logger.info("EXECUTOR --- MAKE_VIEW_JSONS.PY --- Skipping, output already exists")
    else:
        make_view_jsons_tasks = []
        for subdir in os.listdir(output_dir):
            subdir_path = os.path.join(output_dir, subdir)
            if os.path.isdir(subdir_path) and not subdir.startswith("PF"):
                make_view_jsons_tasks.append([
                    python_executable,
                    "make_view_jsons.py",
                    "-sD", subdir_path,
                    "-s", subdir,
                    "-l", timestamped_log
                ])
        Parallel(n_jobs=threads)(
            delayed(run_command)(task, logger)
            for task in make_view_jsons_tasks
        )
        with open(make_view_jsons_done, "w", encoding="utf-8") as f:
            f.write("")
        logger.info("EXECUTOR --- MAKE_VIEW_JSONS.PY --- Executed.")

    logger.info("EXECUTOR --- Pipeline finished successfully")

if __name__ == "__main__":
    main()
