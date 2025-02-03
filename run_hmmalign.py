"""
run_hmmalign.py

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

This script contains a function for running hmmalign and needs 1 argument:
a domain info JSON file with paths to the HMM file, seed alignment, domain fasta and output file.

1 - run_hmmalign - Runs hmmalign for the domain in the domain_info JSON.
"""

import argparse
import logging
import json
import subprocess
from utils import get_logger, get_multi_logger
from typing import Callable
# from modules.decorators import measure_time_and_memory
# from memory_profiler import profile

def parse_arguments():
    """Parse command-line arguments for running hmmalign
    with a JSON. Uses a FASTA containing all a domain's hits
    and stores results per domain in its folder.
    Also includes an output directory and a log file path.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description=
    "Runs hmmalign using a domain's hits sequence database against its HMM \
    aiming to generate a <domain_accession>_hmmalign.sto file inside the domain's folder.")
    parser.add_argument("-iDI", "--dom-info", help="Path to domain info JSON with paths", required=True, type=str)
    parser.add_argument("-d", "--domain-accession", help="Domain accession for scoped logging", required=True, type=str)
    parser.add_argument("-l", "--log", help="Log path", \
        required=False, type=str, default="logs/run_hmmalign.log")
    return parser.parse_args()

#@measure_time_and_memory
#@profile
def run_hmmalign(dom_info_json: str, multi_logger: Callable) -> None:
    """
    Runs hmmalign for the domain in the domain_info JSON.
    """
    with open(dom_info_json, 'r', encoding='utf-8') as dom_info_file:
        dom_info_json = json.load(dom_info_file)

    hmm_file_path = dom_info_json['hmm_file']
    seed_alignment_path = dom_info_json['seed_alignment']
    pfam_id_hmmaligned = dom_info_json['pfam_id_hmmaligned']
    dom_fasta = dom_info_json['dom_fasta']

    command = f"hmmalign --outformat Pfam --mapali {seed_alignment_path} {hmm_file_path} {dom_fasta} > {pfam_id_hmmaligned}"

    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    if result.returncode != 0:
        multi_logger("error", "RUN_HMMALIGN --- Error running hmmalign: %s", result.stderr.decode('utf-8'))
    else:
        multi_logger("info", "RUN_HMMALIGN --- Generated: %s", pfam_id_hmmaligned)

def main(domain_logger: logging.Logger, main_logger: logging.Logger):
    """Main function, initializes this script"""
    args = parse_arguments()
    domain_info_json = args.dom_info

    # Can also get main logger if needed
    main_logger = logging.getLogger("main")
    log_to_both = get_multi_logger([main_logger, domain_logger])
    log_to_both("info", "RUN_HMMALIGN --- Running hmmalign for domain info JSON: %s", domain_info_json)

    run_hmmalign(domain_info_json, log_to_both)

if __name__ == '__main__':
    outer_args = parse_arguments()
    main_logger, _ = get_logger(outer_args.log, scope="main")
    domain_logger, _ = get_logger(outer_args.log, scope="domain", identifier=outer_args.domain_accession)
    main(domain_logger, main_logger)
