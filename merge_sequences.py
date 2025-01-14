"""
merge_sequences.py

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

def merge_sequences(output_dir: str, logger: logging.Logger) -> str:
    """Merges every sequence's PF*_report.json into a single aggregated_report JSON, with structure:
    report[sequence][domain] = {<pair's data>} in the output directory.
    Only processes JSONs in sequence dirs, skips domain dirs.
    For debugging purposes, also stores a list of unique domains in reported_domains.txt.
    Returns the path to aggregated_report.json."""
    aggregated_report_path = os.path.join(output_dir, "aggregated_report.json")
    reported_domains_path = os.path.join(output_dir, "reported_domains.txt")

    # Check if aggregated_report.json already exists
    if os.path.exists(aggregated_report_path):
        logger.error(f"Error: {aggregated_report_path} already exists. Please do not re-run this script, intended use is one-time only, after transfer_annotations.py")
        return

    aggregated_report = {}
    unique_domains = set()

    # Clear the contents of reported_domains.txt before appending new data
    open(reported_domains_path, "w", encoding="utf-8").close()


    for root, _, files in os.walk(output_dir):
        base_dir = os.path.basename(root)
        if base_dir.startswith("PF") or root == output_dir:
            continue
        for file in files:
            if file.endswith("_report.json"):
                domain = file.split("_")[0]
                if domain.startswith("PF") and domain not in unique_domains:
                    unique_domains.add(domain)
                    with open(reported_domains_path, 'a', encoding='utf-8') as domains_path:
                        domains_path.write(domain + "\n")
                report_path = os.path.join(root, file)
                try:
                    with open(report_path, 'r', encoding='utf-8') as report_file:
                        sequence_report = json.load(report_file)
                        sequence_name = sequence_report["sequence_id"]
                        domain_data = sequence_report["domain"]
                        aggregated_report.setdefault(sequence_name, {})
                        aggregated_report[sequence_name].update(domain_data)
                except json.JSONDecodeError:
                    logger.error(f"Failed to parse JSON from {report_path}")
                except Exception as e:
                    logger.error(f"Error processing {report_path}: {str(e)}")

    try:
        with open(aggregated_report_path, "w", encoding="utf-8") as aggregated_report_file:
            json.dump(aggregated_report, aggregated_report_file, indent=4)
        logger.info(f"Successfully wrote {len(unique_domains)} unique domains to {reported_domains_path}")
        logger.info(f"Successfully wrote aggregated report for {len(aggregated_report)} sequences to {aggregated_report_path}")
    except Exception as e:
        logger.error(f"Failed to write aggregated report to {aggregated_report_path} - Error: {str(e)}")

    return aggregated_report_path

def extract_structured_details(details: dict, key_name: str) -> tuple[list, list, list, list]:
    """
    Extracts common structured data (accession, name, count) for a given key.
    Returns (values, accessions, names, counts).
    """
    values = []
    accessions = []
    names = []
    counts = []

    for value, details in details.get(key_name, {}).items():
        values.append(value)
        accessions.append(details.get('rep_primary_accession', ''))
        names.append(details.get('rep_mnemo_name', ''))
        counts.append(str(details.get('count', '')))

    return values, accessions, names, counts

def extract_additional_keys_details(additional_keys: dict) -> tuple[list, list, list, list, list]:
    """Extracts details from additional_keys structure."""
    all_keys = []
    all_values = []
    all_keys_acc = []
    all_keys_rep = []
    all_keys_count = []

    for key, value_dict in additional_keys.items():
        for value, details in value_dict.items():
            all_keys.append(key)
            all_values.append(value)
            all_keys_acc.append(details.get('rep_primary_accession', ''))
            all_keys_rep.append(details.get('rep_mnemo_name', ''))
            all_keys_count.append(str(details.get('count', '')))

    return all_keys, all_values, all_keys_acc, all_keys_rep, all_keys_count


def write_tsv_representation(aggregated_report_path: str, output_tsv_path: str) -> None:
    """Converts aggregated report to TSV."""
    tsv_header = [
        # Core identifiers
        "Sequence ID", "Domain", "Position", "Anno ID", "Hit Status",
        "Type", "Description", "Count",
        # Evidence details
        "Evidence Value", "Evidence Accession", "Evidence Name", "Evidence Count",
        # Paired position details
        "Paired Position", "Paired Accession", "Paired Name", "Paired Count",
        # Additional keys details
        "Additional Key Type", "Additional Value",
        "Additional Accession", "Additional Name", "Additional Count",
        # Conservation data
        "Conservation Score", "Conservation Hit",
        # GO terms
        "GO Annot Names", "GO Terms", "GO Jaccard Index",
        # Indices data
        "Annotation Matches", "Annotation Misses",
        "Conservation Matches", "Conservation Misses"
    ]

    tsv_lines = []
    with open(aggregated_report_path, 'r', encoding="utf-8") as report_file:
        report = json.load(report_file)

        for sequence, domains in report.items():
            for domain, sequence_data in domains.items():
                # Get indices data
                anno_matches = '; '.join(sequence_data.get('annotations', {}).get('indices', {}).get('matches', []))
                anno_misses = '; '.join(sequence_data.get('annotations', {}).get('indices', {}).get('misses', []))
                cons_matches = '; '.join(sequence_data.get('conservations', {}).get('indices', {}).get('matches', []))
                cons_misses = '; '.join(sequence_data.get('conservations', {}).get('indices', {}).get('misses', []))

                for pos, annots in sequence_data.get('annotations', {}).get('positions', {}).items():
                    for anno_id, details in annots.items():
                        # Basic data
                        essentials = details.get('essentials', {})
                        hit_status = "Hit" if details.get('hit', False) else "Miss"

                        # Evidence data
                        ev_values, ev_acc, ev_names, ev_counts = extract_structured_details(details, 'evidence')

                        # Paired position data
                        paired_pos, paired_acc, paired_names, paired_counts = extract_structured_details(details, 'paired_position')

                        # Additional keys data
                        add_types, add_values, add_acc, add_names, add_counts = extract_additional_keys_details(
                            details.get('additional_keys', {})
                        )

                        # Conservation data
                        cons_data = sequence_data.get('conservations', {}).get('positions', {}).get(pos, {})

                        # GO data
                        go_names = []
                        go_terms = []
                        jaccard_indices = []
                        for go_mnemo_name, go_details in details.get('GO', {}).items():
                            go_names.append(go_mnemo_name)
                            go_terms.extend(go_details.get('terms', {}).keys())
                            jaccard_indices.append(str(go_details.get('jaccard_index', '')))

                        row = [
                            sequence, domain, pos, anno_id, hit_status,
                            essentials.get('type', ''), essentials.get('description', ''),
                            str(essentials.get('count', '')),
                            '; '.join(ev_values), '; '.join(ev_acc), '; '.join(ev_names), '; '.join(ev_counts),
                            '; '.join(paired_pos), '; '.join(paired_acc), '; '.join(paired_names), '; '.join(paired_counts),
                            '; '.join(add_types), '; '.join(add_values),
                            '; '.join(add_acc), '; '.join(add_names), '; '.join(add_counts),
                            str(cons_data.get('conservation', '')), str(cons_data.get('hit', '')),
                            '; '.join(go_names), '; '.join(go_terms), '; '.join(jaccard_indices),
                            anno_matches, anno_misses, cons_matches, cons_misses
                        ]
                        tsv_lines.append('\t'.join(row))

    # Write all lines at once after processing
    with open(output_tsv_path, 'w', encoding='utf-8') as tsv_file:
        tsv_file.write('\t'.join(tsv_header) + '\n')
        tsv_file.write('\n'.join(tsv_lines))

def main(logger: logging.Logger):
    """Main function, initializes this script"""
    args = parse_arguments()
    output_dir = args.output_dir
    logger.info(f"Running merge_sequences with arguments: {args}")

    aggregated_report_path = merge_sequences(output_dir, logger)
    output_tsv = os.path.join(output_dir, "initial_visualization.tsv")
    write_tsv_representation(aggregated_report_path, output_tsv_path=output_tsv)

if __name__ == '__main__':
    outer_args = parse_arguments()
    outer_logger = configure_logging(outer_args.log)
    main(outer_logger)