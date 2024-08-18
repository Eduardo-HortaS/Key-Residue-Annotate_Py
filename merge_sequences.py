"""
merge_sequences.py

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

This script contains 2 functions:

1 - merge_sequences.py - merges every sequence's PF*_report.json into a single JSON, with structure:
    report[sequence][domain] = {<pair's data>} and stores in the output directory. Also, makes a reported_domains.txt
    for debugging purposes, this will be removed later. Rquired argument: output dir path.

2 - write_tsv_representation - Converts the aggregated report into a TSV representation and writes it to a file.
    Requires the path to the aggregated report JSON and the path to write the output TSV file.

"""

import os
import argparse
import logging
import json


def parse_arguments():
    """
    Parse command-line arguments for merging every sequence's
    PF*_report.json in a single JSON, with structure:
    report[sequence][domain] = {<pair's data>}
    and stores in the output directory.
    Only needed argument is the output dir path.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """

    parser = argparse.ArgumentParser(description=
    "Merges every sequence's PF*_report.json into a single JSON, with structure: \
    report[sequence][domain] = {<pair's data>} in the output directory.")
    parser.add_argument("-o", "--output-dir", help="Path to output dir", required=True, type=str)
    parser.add_argument("-l", "--log", help="Log path", \
        required=False, type=str, default="logs/merge_sequences.log")
    return parser.parse_args()

def configure_logging(log_path: str) -> logging.Logger:
    """Set up logging for the script."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(filename=log_path, level=logging.DEBUG, filemode='a', \
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    return logging.getLogger()

def merge_sequences(output_dir: str, logger: logging.Logger) -> None:
    """Merges every sequence's PF*_report.json into a single JSON, with structure:
    report[sequence][domain] = {<pair's data>} in the output directory.
    For debugging purposes, also stores a list of unique domains in what_domains.txt.
    Returns the path to aggregated_report.json."""
    aggregated_report_path = os.path.join(output_dir, "aggregated_report.json")
    reported_domains_path = os.path.join(output_dir, "reported_domains.txt")

    # Check if aggregated_report.json already exists
    if os.path.exists(aggregated_report_path):
        logger.error(f"Error: {aggregated_report_path} already exists. Please do not re-run this script, intended use is one-time only, after transfer_annotations.py")
        return

    aggregated_report = {}

    # Clear the contents of what_domains.txt before appending new data
    open(reported_domains_path, "w", encoding="utf-8").close()

    unique_domains = set()

    for root, _, files in os.walk(output_dir):
        for file in files:
            if file.endswith("_report.json"):
                domain = file.split("_")[0]
                if domain.startswith("PF") and domain not in unique_domains:
                    unique_domains.add(domain)
                    with open(reported_domains_path, 'a', encoding='utf-8') as domains_path:
                        domains_path.write(domain + "\n")
                with open(os.path.join(root, file), 'r', encoding='utf-8') as report_file:
                    report = json.load(report_file)
                    sequence = os.path.basename(root)
                    if sequence not in aggregated_report:
                        aggregated_report[sequence] = {}
                    if domain not in aggregated_report[sequence]:
                        aggregated_report[sequence][domain] = report

    with open(aggregated_report_path, "w", encoding="utf-8") as aggregated_report_file:
        json.dump(aggregated_report, aggregated_report_file, indent=4)
    logger.info(f"Wrote unique domains to {reported_domains_path}")
    logger.info(f"Wrote an aggregated report to {aggregated_report_path}")
    return aggregated_report_path

def write_tsv_representation(aggregated_report_path: str, output_tsv_path: str, logger: logging.Logger) -> None:
    """
    Converts the aggregated report into a TSV representation and writes it to a file.
    Requires the path to the aggregated report JSON and the path to write the output TSV file.
    """

    # Initialize the TSV header
    tsv_header = [
        "Sequence ID", "Domain", "Position", "Anno ID", "Count",
        "Evidence - Key", "Evidence - Rep", "Evidence - Count",
        "Type", "Description"
    ]
    additional_keys_set = set()

    # First pass to collect all unique additional keys
    with open(aggregated_report_path, 'r', encoding='utf-8') as aggregated_report_file:
        aggregated_report = json.load(aggregated_report_file)
        for sequence, domains in aggregated_report.items():
            for domain, positions in domains.items():
                for position, annotations in positions.items():
                    for anno_id, details in annotations.items():
                        additional_keys = details.get("additional_keys", {})
                        for key in additional_keys.keys():
                            additional_keys_set.add(key)

    # Add additional keys to the TSV header
    for key in additional_keys_set:
        tsv_header.extend([
            f'"{key}" - Rep', f'"{key}" - Count', f'"{key}" - Value'
        ])

    # Initialize the TSV representation with the header
    tsv_rep = "\t".join(tsv_header) + "\n"

    # Second pass to generate the TSV rows
    with open(aggregated_report_path, 'r', encoding='utf-8') as aggregated_report_file:
        aggregated_report = json.load(aggregated_report_file)
        for sequence, domains in aggregated_report.items():
            for domain, positions in domains.items():
                for position, annotations in positions.items():
                    for anno_id, details in annotations.items():
                        essentials = details.get("essentials", {})
                        evidence = details.get("evidence", {})
                        additional_keys = details.get("additional_keys", {})

                        type_ = essentials.get("type", "")
                        description = essentials.get("description", "")
                        count = essentials.get("count", "")

                        # Handle multiple evidence entries
                        evidence_keys = []
                        evidence_reps = []
                        evidence_counts = []

                        for evidence_key, evidence_value in evidence.items():
                            evidence_keys.append(evidence_key)
                            evidence_reps.append(evidence_value.get("rep_entry_name", ""))
                            evidence_counts.append(evidence_value.get("count", ""))

                        # Join multiple evidence entries with a delimiter
                        evidence_key_str = "; ".join(evidence_keys)
                        evidence_rep_str = "; ".join(evidence_reps)
                        evidence_count_str = "; ".join(map(str, evidence_counts))

                        # Prepare row with basic details
                        row = [
                            sequence, domain, position, anno_id, count,
                            evidence_key_str, evidence_rep_str, evidence_count_str,
                            type_, description
                        ]

                        # Add additional key values, such as ligand_id, ligand_label...
                        for additional_key in additional_keys_set:
                            key_details = additional_keys.get(additional_key, {})
                            rep_entry_name = key_details.get("rep_entry_name", "")
                            key_count = key_details.get("count", "")
                            value = key_details.get("value", "")
                            row.extend([rep_entry_name, key_count, value])

                        # Append the row to the TSV representation
                        tsv_rep += "\t".join(map(str, row)) + "\n"

    with open(output_tsv_path, 'w', encoding='utf-8') as output_tsv_file:
        output_tsv_file.write(tsv_rep)
    logger.info(f"Wrote TSV representation to {output_tsv_path}")

def main(logger: logging.Logger):
    """Main function, initializes this script"""
    args = parse_arguments()
    output_dir = args.output_dir
    logger.info(f"Running merge_sequences with arguments: {args}")

    aggregated_report_path = merge_sequences(output_dir, logger)
    output_tsv = os.path.join(output_dir, "initial_visualization.tsv")
    write_tsv_representation(aggregated_report_path, output_tsv_path=output_tsv, logger=logger)

if __name__ == '__main__':
    outer_args = parse_arguments()
    outer_logger = configure_logging(outer_args.log)
    main(outer_logger)