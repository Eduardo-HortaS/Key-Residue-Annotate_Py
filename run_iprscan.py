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
    3 - run_interproscan - Run InterProScan with the given arguments.

Reference for TSV output columns:
https://interproscan-docs.readthedocs.io/en/latest/OutputFormats.html#tab-separated-values-format-tsv

"""

import argparse
import logging
import subprocess
from typing import Callable
from utils import get_logger, get_multi_logger

def parse_arguments():
    """Parse command-line arguments for running InterProScan."""
    parser = argparse.ArgumentParser(description="Runs InterProScan for a sequence in a fasta input.")
    parser.add_argument("-iP", "--iprscan-path", help="Path to interproscan.sh", required=True, type=str)
    parser.add_argument("-iF", "--input-fasta", help="Path to input sequence fasta", required=True, type=str)
    parser.add_argument("-s", "--sequence", help="Sequence to process", required=True, type=str)
    parser.add_argument("-oB", "--output-base-file", help="Path to output base file, ending right before the extension", required=True, type=str)
    parser.add_argument("-oF", "--output-format", help="Output format/s, from (TSV, JSON, XML, GFF3). This pipeline only uses TSV, but the default is TSV, XML and GFF3 for possible downstream analyses", required=False, type=str, default="TSV,XML,GFF3")
    parser.add_argument("-dB", "--databases", help="Optional: Comma-separated databases to limit the search, if the single purpose is using GO terms inside the pipeline, pass the string panther,gene3d,smart,pfam,superfamily", required=False, type=str, default="")
    parser.add_argument("-l", "--log", help="Log path", required=False, type=str, default="logs/run_iprscan.log")
    args = parser.parse_args()
    args.output_format = ', '.join(fmt.strip() for fmt in args.output_format.split(','))
    return args

def run_interproscan(iprscan_path: str, input_fasta: str, output_basefile: str, output_format: str, databases: str, multi_logger: Callable) -> None:
    """
    Runs InterProScan with the given arguments.
    """
    # NOTE: Kinda slow, not sure about CPU usage.
    command = f"{iprscan_path} -i {input_fasta} -b {output_basefile} -f {output_format} -goterms -iprlookup -dp --cpu 1"

    if databases:
        command += f" -appl {databases}"

    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    if result.returncode != 0:
        multi_logger("error", "RUN_IPRSCAN --- Error running InterProScan: %s", result.stderr.decode('utf-8'))
    else:
        multi_logger("info", "RUN_IPRSCAN --- InterProScan completed successfully. Output saved to %s.%s", output_basefile, output_format.lower())


def main(sequence_logger: logging.Logger):
    """Main function, initializes this script"""
    args = parse_arguments()
    iprscan_sh_path = args.iprscan_path
    input_fasta_filepath = args.input_fasta
    output_basefile = args.output_base_file
    output_format = args.output_format
    databases = args.databases
    main_logger = logging.getLogger("main")
    log_to_both = get_multi_logger([main_logger, sequence_logger])
    log_to_both("info", "RUN_IPRSCAN --- Running InterProScan with input fasta: %s",  input_fasta_filepath)

    run_interproscan(iprscan_sh_path, input_fasta_filepath, output_basefile, output_format, databases, log_to_both)

if __name__ == '__main__':
    outer_args = parse_arguments()
    outer_logger, _ = get_logger(outer_args.log, scope="sequence", identifier=outer_args.sequence)
    main(outer_logger)
