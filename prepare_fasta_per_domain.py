"""
prepare_fasta_per_domain.py

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

This script contains 2 functions, one sees if a domain in the hits per domain JSON is suitable,
the other prepares a fasta for it if so. It needs 4 arguments:
a hmmsearch hits per domain JSON, the domain accession to prep
and the paths to the resource directory and output directory.

1 - can_run_hmmalign - Checks if intermediary files are present for the given domain accession in resource dir.
If so, call prep_domain_fasta.

2 - prep_domain_fasta - Accesses the JSON in search of the given accession and makes a multifasta with all hits contained in it.

Obs.: It'll make a subdir for each valid domain in the output directory. Also, it'll put the substring
"target/" between target_seq_name and ali range to facilitate parsing in the transfer_annotations step:
signalling that substring denotes a target sequence versus the seed sequences.

"""

import os
import json
import argparse
import logging
from typing import Any, Callable
from utils import get_logger, get_multi_logger
# from modules.decorators import measure_time_and_memory

def parse_arguments():
    """Parse command-line arguments for finding out if a given domain from those
    in the hmmsearch hits per domain JSON is suitable - presents required intermediary files.
    Also includes paths to the resource directory and output directory, and a log file path.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description=
    'Generates a temporary multifasta for running hmmalign using a hits per domain JSON.')
    parser.add_argument("-iJ", "--per-dom-json", help="Path to hits per domain json", required=True, type=str)
    parser.add_argument("-iD", "--domain-accession", help="The domain Pfam accession you're prepping for", required=True, type=str)
    parser.add_argument("-r", "--resource-dir", help="Resource dir path", required=True, type=str)
    parser.add_argument("-o", "--output-dir", help="Output dir path", required=True, type=str)
    parser.add_argument("-l", "--log", help="Log path", \
        required=False, type=str, default="logs/prepare_fasta_per_domain.log")
    return parser.parse_args()

def can_run_hmmalign(dom_accession: str, resource_dir: str, output_dir: str) -> dict[str, Any]:
    """
    For a given domain accession, checks if the necessary resource files are present.
    Required: HMM file and seed alignment
    Optional but need at least one: conservations or annotations file
    """
    if '.' in dom_accession:
        dom_accession = dom_accession.split('.')[0]

    hmm_file_path = os.path.join(resource_dir, dom_accession, "domain.hmm")
    seed_alignment_path = os.path.join(resource_dir, dom_accession, "alignment.seed")
    conservations_file_path = os.path.join(resource_dir, dom_accession, "conservations.json")
    annotations_file_path = os.path.join(resource_dir, dom_accession, "annotations.json")
    pfam_id_hmmaligned = os.path.join(output_dir, dom_accession, dom_accession + "_hmmalign.sth")

    # Required files check
    required_files_exist = (os.path.isfile(hmm_file_path) and
                          os.path.isfile(seed_alignment_path))

    # At least one optional file must exist
    optional_file_exists = (os.path.isfile(conservations_file_path) or
                          os.path.isfile(annotations_file_path))

    domain_run_info = {
        "can_align": required_files_exist and optional_file_exists,
        "hmm_file": hmm_file_path if os.path.isfile(hmm_file_path) else None,
        "seed_alignment": seed_alignment_path if os.path.isfile(seed_alignment_path) else None,
        "conservations": conservations_file_path if os.path.isfile(conservations_file_path) else None,
        "annotations": annotations_file_path if os.path.isfile(annotations_file_path) else None,
        "pfam_id_hmmaligned": pfam_id_hmmaligned
    }
    return domain_run_info

def prep_domain_fasta(per_dom_json: str, dom_accession: str, output_dir: str, domain_logger: logging.Logger, multi_logger: Callable) -> (str | None):
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
        multi_logger("error", "PREPARE_FASTA_PER_DOMAIN --- Error opening or reading file %s: %s", per_dom_json, e)
        return None

    domain_logger.info("PREPARE_FASTA_PER_DOMAIN --- Preparing fasta for domain %s", dom_accession)

    for accession, sequences in hits.items():
        if accession == dom_accession:
            for sequence_hits in sequences.values():
                for hit in sequence_hits:
                    subseq = hit.get('subseq', '')
                    ali_range = hit.get('ali_range', '')
                    target_seq_name = hit.get('target_seq_name', '')
                    header = f">{target_seq_name}target/{ali_range}"
                    fasta_data += f"{header}\n{subseq}\n"

    if not fasta_data:
        multi_logger("warning", "PREPARE_FASTA_PER_DOMAIN --- No hits found for domain %s", dom_accession)
        return None

    domain_dir = os.path.join(output_dir, dom_accession)
    os.makedirs(domain_dir, exist_ok=True)

    fasta_filename = f"{dom_accession}_hits.fasta"
    fasta_path = os.path.join(domain_dir, fasta_filename)

    try:
        with open(fasta_path, 'w', encoding='utf-8') as fasta_file:
            fasta_file.write(fasta_data)
        domain_logger.info("PREPARE_FASTA_PER_DOMAIN --- Generated fasta for domain %s at %s", dom_accession, fasta_path)
        return fasta_path
    except IOError as e:
        multi_logger("error", "PREPARE_FASTA_PER_DOMAIN --- Error writing FASTA file %s: %s", fasta_path, e)
        return None

def main():
    """Main function, initializes this script"""
    args = parse_arguments()
    per_dom_json = args.per_dom_json
    dom_accession = args.domain_accession
    resource_dir = args.resource_dir
    output_dir = args.output_dir
    log_path = args.log

    main_logger, _ = get_logger(log_path, scope="main")
    domain_logger, _ = get_logger(log_path, scope="domain", identifier=dom_accession)
    log_to_both = get_multi_logger([main_logger, domain_logger])

    domain_logger.info("PREPARE_FASTA_PER_DOMAIN --- Running prepare_fasta_per_domain with arguments: %s", args)

    domain_info = can_run_hmmalign(dom_accession, resource_dir, output_dir)
    if domain_info['can_align']:
        dom_fasta_path = prep_domain_fasta(per_dom_json, dom_accession, output_dir, domain_logger, log_to_both)
        if dom_fasta_path:
            domain_info['dom_fasta'] = dom_fasta_path
            output_json_path = os.path.join(output_dir, dom_accession, 'domain_info.json')
            os.makedirs(os.path.dirname(output_json_path), exist_ok=True)
            try:
                with open(output_json_path, 'w', encoding='utf-8') as f:
                    json.dump(domain_info, f, indent=4)
                domain_logger.info("PREPARE_FASTA_PER_DOMAIN --- Information for  %s was written to %s", dom_accession, output_json_path)
            except IOError as e:
                log_to_both("error", "PREPARE_FASTA_PER_DOMAIN --- Error writing domain info to %s: %s", output_json_path, e)
    else:
        log_to_both("warning", "PREPARE_FASTA_PER_DOMAIN --- Missing required files for domain %s", dom_accession)

if __name__ == '__main__':
    main()
