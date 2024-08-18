"""
prepare_fasta_per_domain.py

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

This script contains 2 functions, one sees suitable domains,
the other prepares a fasta per each domain in the resulting hits per domain JSON.
It needs 4 arguments:
a hits per domain JSON file, the domain accession to prep and the paths to the resource directory and output directory.

1 - can_run_hmmalign - Checks if intermediary files are present for the given domain accession. If so, call prep_domain_fasta.

2 - prep_domain_fasta - Accesses the JSON in search of the given accession and makes a multifasta with all hits contained in it.

Obs.: It'll make a subdir for each domain in the output directory. Also, it'll put the substring
"target/" between target_seq_name and ali range to facilitate parsing in the next step.
That is, that substring denotes a target sequence versus the seed sequences.

"""

import os
import json
import argparse
import logging
from typing import Any
# from modules.decorators import measure_time_and_memory

def parse_arguments():
    """Parse command-line arguments for running hmmsearch
    using a sequence database against target HMMs,
    aiming to generate a hmmsearch_per_domain_output.json.
    Also includes an output directory and a log file path.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description=
    'Generates a temporary multifasta for running hmmalign using a hits per domain JSON.')
    parser.add_argument("-iJ", "--per-dom-json", help="Path to hits per domain json", required=True, type=str)
    parser.add_argument("-iD", "--dom-accession", help="The domain Pfam accession you're prepping for", required=True, type=str)
    parser.add_argument("-r", "--resource-dir", help="Resource dir path", required=True, type=str)
    parser.add_argument("-o", "--output-dir", help="Output dir path", required=True, type=str)
    parser.add_argument("-l", "--log", help="Log path", \
        required=False, type=str, default="logs/prepare_fasta_per_domain.log")
    return parser.parse_args()

def configure_logging(log_path: str) -> logging.Logger:
    """Set up logging for the script."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(filename=log_path, level=logging.DEBUG, filemode='a', \
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    return logging.getLogger()

def can_run_hmmalign(dom_accession: str, resource_dir: str, output_dir: str) -> dict[str, Any]:
    """
    For a given domain accession, checks if the necessary resource files are present.
    Store this and paths to these files and to output in a dictionary,
    returned for outside determination of what domains to prep a fasta for.
    """
    if '.' in dom_accession:
        dom_accession = dom_accession.split('.')[0]
    # print(f"Dom access in can run {dom_accession}")
    hmm_file_path = os.path.join(resource_dir, dom_accession, "domain.hmm")
    seed_alignment_path = os.path.join(resource_dir, dom_accession, "alignment.seed")
    pfam_id_hmmaligned = os.path.join(output_dir, dom_accession, dom_accession + "_hmmalign.sth")

    # Check for file existence and prepare the dictionary
    domain_run_info = {
            "can_align": os.path.isfile(hmm_file_path) and os.path.isfile(seed_alignment_path),
            "hmm_file": hmm_file_path if os.path.isfile(hmm_file_path) else None,
            "seed_alignment": seed_alignment_path if os.path.isfile(seed_alignment_path) else None,
            "pfam_id_hmmaligned": pfam_id_hmmaligned
    }
    return domain_run_info

def prep_domain_fasta(per_dom_json: str, dom_accession: str, output_dir: str, logger: logging.Logger) -> (str | None):
    """
    Loads a hits per domain JSON and searches for a target domain by its accession to
    generate a FASTA containing its hits across all sequences.
    Each domain's files are stored in a domain subdirectory within the output directory.
    """
    fasta_data = ""
    try:
        with open(per_dom_json, 'r', encoding='utf-8') as f:
            hits = json.load(f)
    except IOError as e:
        logger.error(f"Error opening or reading the file: {e}")
        return

    for accession, sequences in hits.items():
        if accession == dom_accession:
            for sequence_hits in sequences.values():
                for hit in sequence_hits:
                    subseq = hit.get('subseq', '')
                    ali_range = hit.get('ali_range', '')
                    target_seq_name = hit.get('target_seq_name', '')
                    header = f">{target_seq_name}target/{ali_range}"
                    fasta_data += f"{header}\n{subseq}\n"

    domain_dir = os.path.join(output_dir, dom_accession)
    os.makedirs(domain_dir, exist_ok=True)

    fasta_filename = f"{dom_accession}_hits.fasta"
    fasta_path = os.path.join(domain_dir, fasta_filename)

    with open(fasta_path, 'w', encoding='utf-8') as fasta_file:
        fasta_file.write(fasta_data)
    return fasta_path

def main(logger: logging.Logger):
    """Main function, initializes this script"""
    args = parse_arguments()
    per_dom_json = args.per_dom_json
    dom_accession = args.dom_accession
    resource_dir = args.resource_dir
    output_dir = args.output_dir

    logger.info(f"Running prepare_fasta_per_domain with arguments: {args}")

    domain_info = can_run_hmmalign(dom_accession, resource_dir, output_dir)
    if domain_info['can_align']:
        dom_fasta_path = prep_domain_fasta(per_dom_json, dom_accession, output_dir, logger)
        domain_info['dom_fasta'] = dom_fasta_path
        output_json_path = os.path.join(output_dir, dom_accession, 'domain_info.json')
        os.makedirs(os.path.dirname(output_json_path), exist_ok=True)
        try:
            with open(output_json_path, 'w', encoding='utf-8') as f:
                json.dump(domain_info, f, indent=4)
            logger.info(f"Information for {dom_accession} was written to {output_json_path}")
        except IOError as e:
            logger.error(f"Error writing the file: {e}")
    else:
        logger.warning(f"Couldn't find the necessary files for domain {dom_accession}")

if __name__ == '__main__':
    outer_args = parse_arguments()
    outer_logger = configure_logging(outer_args.log)
    main(outer_logger)
