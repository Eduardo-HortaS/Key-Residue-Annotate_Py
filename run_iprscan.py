"""
run_iprscan.py

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

Runs InterProScan for a sequence in a fasta input.
The user may opt (-db) for the default usage, a run including the panther, gene3d, smart,
pfam, and superfamily databases, or to run with a custom set of available databases (-db=<comma-separated databases string>).

This script uses 6 arguments: required paths to interproscan.sh, input sequence fasta and
output base file name (absolute or relative),
besides optional output format comma-and-space separated string (e.g. TSV, JSON, XML, GFF3)
and/or databases to use comma-separated string (e.g. pfam,gene3d,smart) and log filepath.
The default output format is TSV.

Command Format:

<iprscan_path> -i <input_fasta> -b <output_base_file> -f <output_format> -appl <databases> -goterms -iprlookup -dp

Functions:
    1 - parse_arguments - Parse command-line arguments for running InterProScan.
    2 - create_batch_fasta - Creates a FASTA file containing all sequences in a batch.
    3 - run_interproscan - Run InterProScan with the given arguments.
    4 - split_iprscan_output - Splits InterProScan output files by sequence.

Reference for TSV output columns:
https://interproscan-docs.readthedocs.io/en/latest/OutputFormats.html#tab-separated-values-format-tsv

"""

import os
import sys
import argparse
import logging
import subprocess
from typing import Callable
from utils import get_logger, get_multi_logger

def parse_arguments():
    """Parse command-line arguments for running InterProScan."""
    parser = argparse.ArgumentParser(
        description="Runs InterProScan for a sequences in \
        either single file or batch mode.")
    subparsers = parser.add_subparsers(dest="mode", required=True)

    # Required arguments
    parser.add_argument("-iPr", "--iprscan-path",
                        help="Path to interproscan.sh",
                        required=True, type=str)

    # Single File Mode
    single_file_parser = subparsers.add_parser("single")
    single_file_parser.add_argument("-f", "--fasta",
                        help="Path to input sequence FASTA (Single File Mode)")
    single_file_parser.add_argument("-s", "--sequence",
                        help="Sequence identifier for scoped logging")

    # Batch Mode
    batch_file_parser = subparsers.add_parser("batch")
    batch_file_parser.add_argument("-sB", "--sequence-batch",
                        help="Comma-separated list of sequence IDs to process")
    batch_file_parser.add_argument("-sBi", "--sequence-batch-index",
                        help="Index of this batch (for naming)")
    batch_file_parser.add_argument("-sPd", "--sequence-parent-dir",
                        help="Parent directory of sequence subdirectories")

    # Optional arguments
    parser.add_argument("-iA", "--analyses",
                        help="Optional: Comma-separated analyses to limit InterProScan run. \
                        If you just want what's used in the pipeline (GO terms), \
                        pass the string: panther,gene3d,smart,pfam,superfamily",
                        required=False, type=str, default="")
    # Reminder: action="store_true" is a keyword argument that defaults to False,
    # hence no default value is needed.
    parser.add_argument("-iDpc", "--enable-precalc", action="store_true",
                        help="Flag to enable pre-calculated match lookup in InterProScan. \
                        Keep this off unless running with proteins in UniProtKB.",
                        required=False)
    parser.add_argument("-iDr", "--disable-res", action="store_true",
                        help="Flag to disable residue-level annotations in \
                        InterProScan outputs, may increase performance, \
                        but possible miss some cross-match annotations.",
                        required=False)
    parser.add_argument("-iOf", "--output-format",
                        help="Output format/s, from (TSV, JSON, XML, GFF3). \
                        This pipeline only uses TSV, but the default is \
                        TSV, XML and GFF3 for possible downstream analyses",
                        required=False, type=str, default="TSV,XML,GFF3")
    parser.add_argument("-iCc", "--cpu-cores",
                        help="Number of CPU cores to use per job",
                        required=False, type=int, default=8)
    parser.add_argument("-l", "--log",
                        help="Log path",
                        required=False, type=str, default="logs/run_iprscan.log")
    args = parser.parse_args()

    if args.mode == "single":
        if not args.fasta:
            parser.error("Single mode requires --fasta")
        if not args.sequence:
            parser.error("Single mode requires --sequence")

    if args.mode == "batch":
        if not args.sequence_batch:
            parser.error("Batch mode requires --sequence-batch")
        if not args.sequence_batch_index:
            parser.error("Batch mode requires --sequence")
        if not args.sequence_parent_dir:
            parser.error("Batch mode requires --sequence-parent-dir")

    # Validate batch mode arguments
    if args.sequence_batch and not args.sequence_parent_dir:
        parser.error("--sequence-batch requires --sequence-parent-dir")
    if args.sequence_parent_dir and not args.sequence_batch:
        parser.error("--sequence-parent-dir requires --sequence-batch")

    args.output_format = ', '.join(fmt.strip() for fmt in args.output_format.split(','))

    return args

def create_batch_fasta(
    sequence_parent_dir: str, sequence_batch_dash: list[str], batch_idx: int,
    logger: logging.Logger, multi_logger: Callable) -> str:
    """Creates a FASTA file containing all sequences in the batch.

    Args:
        sequence_parent_dir: Parent directory containing sequence subdirectories
        sequence_batch_dash: List of sequence IDs in this batch, using dashes instead of pipes
        batch_idx: Index of this batch (for naming)
        logger: Logger function

    Returns:
        str: Path to created batch FASTA file
    """
    # Create batches directory
    batches_dir = os.path.join(sequence_parent_dir, "batches")
    os.makedirs(batches_dir, exist_ok=True)

    # Create batch filename
    first_seq = sequence_batch_dash[0]
    last_seq = sequence_batch_dash[-1]
    batch_name = f"batch_{batch_idx}_{first_seq}_to_{last_seq}.fasta"
    batch_path = os.path.join(batches_dir, batch_name)

    try:
        with open(batch_path, "w", encoding="utf-8") as batch_file:
            for seq_id in sequence_batch_dash:
                seq_fasta = os.path.join(sequence_parent_dir, seq_id, "sequence.fasta")
                if not os.path.exists(seq_fasta):
                    multi_logger("warning", "RUN_IPRSCAN --- CREATE_BATCH_FASTA --- Sequence file not found: %s", seq_fasta)
                    continue

                with open(seq_fasta, "r", encoding="utf-8") as f:
                    batch_file.write(f.read())

        logger.info(f"RUN_IPRSCAN --- CREATE_BATCH_FASTA --- Created batch FASTA at {batch_path}")
        return batch_path

    except Exception as e:
        if os.path.exists(batch_path):
            os.unlink(batch_path)
        multi_logger("error", "RUN_IPRSCAN --- CREATE_BATCH_FASTA --- Failed to create batch FASTA: %s", e)
        raise

def split_iprscan_output(
    output_base: str, sequence_batch_pipe: list[str], formats: list[str],
    logger: logging.Logger, multi_logger: Callable) -> None:
    """Splits InterProScan output files by sequence.

    Args:
        output_base: Base path of InterProScan output files
        sequence_ids: List of sequence IDs in the batch
        formats: List of output formats to process
        logger: Logger function
    """
    logger.info("RUN_IPRSCAN --- SPLIT_IPRSCAN --- Using sequence_ids: %s", sequence_batch_pipe)
    for fmt in formats:
        fmt = fmt.strip().lower()
        input_file = f"{output_base}.{fmt}"
        logger.info("RUN_IPRSCAN --- SPLIT_IPRSCAN --- Splitting InterProScan output file: %s", input_file)
        if not os.path.exists(input_file):
            multi_logger("error", "RUN_IPRSCAN --- SPLIT_IPRSCAN --- \
            InterProScan output file not found: %s", input_file)
            continue

        if fmt == "tsv":
            logger.info("RUN_IPRSCAN --- SPLIT_IPRSCAN --- Splitting TSV output")
            # TSV format - split by sequence ID in first column
            with open(input_file, "r", encoding="utf-8") as f:
                # Group lines by sequence ID
                sequence_outputs = {}
                for line in f:
                    if line.strip():
                        seq_id = line.split("\t")[0]
                        logger.info(f"RUN_IPRSCAN --- SPLIT_IPRSCAN --- Sequence ID {seq_id}")
                        if seq_id in sequence_batch_pipe:
                            sequence_outputs.setdefault(seq_id, []).append(line)

                # Write individual sequence files
                for seq_id, lines in sequence_outputs.items():
                    sanitized_seq_id = seq_id.replace("|", "-")
                    seq_dir = os.path.join(os.path.dirname(os.path.dirname(output_base)), sanitized_seq_id)
                    os.makedirs(seq_dir, exist_ok=True)
                    out_file = os.path.join(seq_dir, f"iprscan.{fmt}")
                    with open(out_file, "w", encoding="utf-8") as f:
                        f.writelines(lines)
            logger.info("RUN_IPRSCAN --- SPLIT_IPRSCAN --- \
            Split TSV output for %d sequences", len(sequence_outputs))

        elif fmt in ["xml", "json", "gff3"]:
            logger.info("RUN_IPRSCAN --- SPLIT_IPRSCAN --- \
            Splitting %s output is not yet supported, may be in the future \
            if sufficient need is made evident", fmt.upper())

        else:
            multi_logger("warning", "RUN_IPRSCAN --- SPLIT_IPRSCAN --- \
            Unknown or unsupported InterProScan output format: %s", fmt.upper())

    logger.info("RUN_IPRSCAN --- SPLIT_IPRSCAN --- Finished splitting InterProScan output files")

def run_interproscan(
    iprscan_path: str, enable_precalc: bool, disable_res: bool,
    input_fasta: str, output_basefile: str, output_format: str,
    analyses: str, cpu_cores: int,  multi_logger: Callable) -> None:
    """
    Runs InterProScan with the given arguments.

    Args:
        iprscan_path: Path to interproscan.sh
        enable_precalc: Flag to enable pre-calculated match lookup in InterProScan
        disable_res: Flag to disable residue-level annotations in InterProScan outputs
        input_fasta: Path to input sequence FASTA
        output_basefile: Base path for writing InterProScan output files
        output_format: Output format/s, from (TSV, JSON, XML, GFF3)
        analyses: Comma-separated analyses to limit InterProScan run
        cpu_cores: Number of CPU cores to use per job
        multi_logger: Callable to send warning, error and critical messages to both a main and sequence|batch logger
    """
    formats = [fmt.strip().lower() for fmt in output_format.split(',')]
    additional_formats = [fmt for fmt in formats if fmt in ["xml", "json", "gff3"]]

    if 'tsv' not in formats:
        multi_logger("warning", "RUN_IPRSCAN --- RUN --- Adding TSV output format \
        since it's required for GO term processing in transferring annotations")
        formats.append('tsv')

    if additional_formats:
        multi_logger("info", "RUN_IPRSCAN --- RUN --- Additional output formats detected: (%s). \
        These will be available as-is for downstream analyses but are not processed by this pipeline",
        ', '.join(additional_formats).upper())

    output_format = ', '.join(formats)
    # NOTE: Need to be careful with CPU cores and memory allocated,
    # InterProScan can be very intensive in these regards.
    cmd = f"{iprscan_path} -i {input_fasta} -b {output_basefile} -f {output_format} \
    -goterms -iprlookup -dra --cpu {cpu_cores}"

    if analyses:
        cmd += f" -appl {analyses}"
    if not enable_precalc:
        cmd += " -dp"
    if disable_res:
        cmd += " -dra"
    else:
        cmd += " -etra" # Includes sites in TSV output

    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=True)

    if result.returncode != 0:
        multi_logger("error", "RUN_IPRSCAN --- RUN --- Error running InterProScan: %s",
        result.stderr.decode('utf-8'))
    else:
        multi_logger("info", "RUN_IPRSCAN --- RUN --- InterProScan completed successfully. \
        Output saved to %s.%s", output_basefile, output_format.lower())


def main():
    """Main function, initializes this script"""
    args = parse_arguments()
    main_logger, _ = get_logger(args.log, scope="main")
    sequence_logger = None
    if args.mode == "single":
        sequence_logger, _ = get_logger(args.log, scope="sequence", identifier=args.sequence)
    elif args.mode == "batch":
        sequence_logger, _ = get_logger(
            args.log,
            scope="seq_batch",
            identifier=f"batch_{args.sequence_batch_index}" if hasattr(args, 'sequence_batch_index') else None
        )
    if sequence_logger:
        log_to_both = get_multi_logger([main_logger, sequence_logger])
    else:
        main_logger.error("RUN_IPRSCAN --- MAIN --- Sequence logger not initialized, \
        necessary for more specific logging messages. \
        Absence indicates issue with missing batch arguments. Exiting.")
        sys.exit(1)

    input_fasta = None
    sequence_batch_pipe = None

    if args.mode == "single":
        # Single File Mode
        input_fasta = args.fasta
    elif args.mode == "batch":
        # Batch Mode
        sequence_batch_dash = [sb.strip() for sb in args.sequence_batch.split(',')]
        sequence_batch_pipe = [sb.replace("-", "|") for sb in sequence_batch_dash]
        input_fasta = create_batch_fasta(
            sequence_parent_dir=args.sequence_parent_dir,
            sequence_batch_dash=sequence_batch_dash,
            batch_idx=args.sequence_batch_index,
            logger=sequence_logger,
            multi_logger=log_to_both
        )

    if input_fasta and args.mode == "single":
        output_base = os.path.join(os.path.dirname(input_fasta), "iprscan")
    elif input_fasta and args.mode == "batch":
        output_base = os.path.join(os.path.dirname(input_fasta), f"iprscan_batch_{args.sequence_batch_index}")
    else:
        log_to_both("error", "RUN_IPRSCAN --- MAIN --- No input fasta file found, \
        cannot run InterProScan without sequences. Exiting.")
        sys.exit(1)

    # Run InterProScan
    run_interproscan(
        iprscan_path=args.iprscan_path,
        enable_precalc=args.enable_precalc,
        disable_res=args.disable_res,
        input_fasta=input_fasta,
        output_basefile=output_base,
        output_format=args.output_format,
        analyses=args.analyses,
        cpu_cores=args.cpu_cores,
        multi_logger=log_to_both
    )

    # Split outputs if in batch mode
    if args.mode == "batch" and sequence_batch_pipe:
        formats = [fmt.strip() for fmt in args.output_format.split(',')]
        split_iprscan_output(
            output_base=output_base,
            sequence_batch_pipe=sequence_batch_pipe,
            formats=formats,
            logger=sequence_logger,
            multi_logger=log_to_both)

if __name__ == '__main__':
    main()
