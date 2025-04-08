"""
merge_reports_in_sequences.py

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

This script merges multiple domain-specific PF*_report.json files for a sequence
into a single aggregated_report.json file. It is designed to run after transfer_annotations.py.

The main function:
merge_sequences - combines a sequence's domain reports into a single JSON with structure:
report[sequence][domain] = {<pair's data>} and stores as aggregated_report.json in the sequence directory.

Required command-line arguments:
- sequence: Sequence identifier for scoped logging
- sequence-dir: Path to the sequence directory containing PF*_report.json files
"""

import os
import argparse
import logging
import json
from typing import Callable
from utils import get_logger, get_multi_logger

def parse_arguments():
    """
    Parse command-line arguments for merging a sequence's
    PF*_report.json files into a single JSON, with structure:
    report[sequence][domain] = {<pair's data>}
    and stores in the sequence directory.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """

    parser = argparse.ArgumentParser(description=
    "Merges a sequence's PF*_report.json files into a single JSON, with structure: \
    report[sequence][domain] = {<pair's data>} in the sequence directory.")
    parser.add_argument("-s", "--sequence", help="Sequence identifier for scoped logging", required=True, type=str)
    parser.add_argument("-sd", "--sequence-dir", help="Sequence directory within output dir", required=True, type=str)
    parser.add_argument("-l", "--log", help="Log path", required=False, type=str, default="logs/merge_sequences.log")
    return parser.parse_args()

def merge_sequences(sequence_dir: str, multi_logger: Callable, logger: logging.Logger) -> str:
    """Merges a sequence's PF*_report.json files into a single aggregated_report JSON, with structure:
    report[sequence][domain] = {<pair's data>} in the sequence directory.
    Returns the path to aggregated_report.json."""

    aggregated_report_path = os.path.join(sequence_dir, "aggregated_report.json")

    # Check if aggregated_report.json already exists
    if os.path.exists(aggregated_report_path):
        multi_logger("error", "Error: %s already exists. Please do not re-run this script, intended use is one-time only, after transfer_annotations.py", aggregated_report_path)
        return None

    aggregated_report = {}

    for file in os.listdir(sequence_dir):
        if file.endswith("_report.json"):
            report_path = os.path.join(sequence_dir, file)
            try:
                with open(report_path, 'r', encoding='utf-8') as report_file:
                    sequence_report = json.load(report_file)
                    sequence_name = sequence_report["sequence_id"]
                    domain_data = sequence_report["domain"]
                    aggregated_report.setdefault(sequence_name, {})
                    aggregated_report[sequence_name].update(domain_data)
            except json.JSONDecodeError:
                multi_logger("error", "Failed to parse JSON from %s", report_path)
            except Exception as e:
                multi_logger("error", "Error processing %s: %s", report_path, str(e))

    try:
        with open(aggregated_report_path, "w", encoding="utf-8") as aggregated_report_file:
            json.dump(aggregated_report, aggregated_report_file, indent=4)
        logger.info("Successfully wrote aggregated report for sequence to %s", aggregated_report_path)
    except Exception as e:
        multi_logger("error", "Failed to write aggregated report to %s - Error: %s", aggregated_report_path, str(e))

    return aggregated_report_path

def main():
    """Main function, initializes this script"""
    args = parse_arguments()
    sequence_dir = args.sequence_dir
    main_logger, _ = get_logger(args.log, scope="main")
    sequence_logger, _ = get_logger(args.log, scope="sequence", identifier=args.sequence)
    log_to_both = get_multi_logger([main_logger, sequence_logger])

    log_to_both("info", "MERGE_SEQUENCES --- Running merge_sequences with arguments: %s", args)
    merge_sequences(sequence_dir, log_to_both, sequence_logger)

if __name__ == '__main__':
    main()
