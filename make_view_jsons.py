"""
make_view_jsons.py

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

This script contains functions to create JSON structures properly formatted
for consumption by Nightingale components in the report app.

Functions:
1. process_sequence_report: Reads per-sequence report JSONs and processes annotations
2. transform_to_ranges: Converts position-based annotations to range-based format
3. aggregate_range_annotations: Combines annotation data for continuous ranges
4. write_range_views: Outputs data in a Nightingale-compatible format

Input: Per-sequence JSON reports from transfer_annotations.py
Output: Range-based JSON views for Nightingale visualization
"""


import os
import sys
import json
import logging
import argparse
from collections import defaultdict
from typing import Tuple, Dict, Callable
from utils import convert_lists_to_original_types, convert_sets_and_tuples_to_lists, convert_defaultdict_to_dict, get_logger, get_multi_logger

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Transform position-based annotations into range-based views")
    parser.add_argument("-sD", "--sequence-dir", help="Path to sequence directory with report.jsons", required=True)
    parser.add_argument("-s", "--sequence", help="Sequence to process", required=True)
    parser.add_argument("-l", "--log", help="Log file path", default="logs/make_views_jsons.log")
    return parser.parse_args()

def track_position_data(target: Dict, key: str, pos_str: str, value: Dict, track_key: str) -> None:
    """Helper to track per-position data like counts and hits."""
    if track_key not in target[key]:
        target[key][track_key] = {}
    target[key][track_key][pos_str] = value[track_key]

def merge_nested_data(source: Dict, target: Dict, pos_str: str) -> None:
    """Recursively merge nested dictionaries at range_id, tracking per-position data."""
    if "hit" in source:
        if "hit" not in target:
            target["hit"] = {}
        target["hit"][pos_str] = source["hit"]

    for key, value in source.items():
        if isinstance(value, dict):
            if key not in target:
                target[key] = defaultdict(dict)

            # Track position-specific data
            if "count" in value:
                track_position_data(target, key, pos_str, value, "count")

            if key == "GO":
                if pos_str not in target["GO"]:
                    target["GO"][pos_str] = {}
                target["GO"][pos_str].update(value)
                continue

            # Recurse for nested structures
            merge_nested_data(value, target[key], pos_str)
        else:
            if key not in ["count", "hit", "GO"]:  # Skip specially handled keys
                target[key] = value

def aggregate_range_positions(interval_dict: Dict, range_id: str, range_tuple: Tuple[int, int],
                          data_type: str) -> Tuple[Dict, Dict]:
    """Aggregate position data and metadata for a range."""
    start, end = range_tuple
    aggregated = defaultdict(lambda: defaultdict(dict))
    metadata = {}

    # Get metadata
    if data_type == "annotations":
        metadata = {k: v for k, v in interval_dict["annotations"].items()
                   if k != "positions"}
    elif data_type == "conservations":
        metadata = {k: v for k, v in interval_dict["conservations"].items()
                   if k != "positions"}

    # Process positions
    if data_type == "annotations":
        for pos in range(start, end + 1):
            pos_str = str(pos)
            if pos_str not in interval_dict["annotations"]["positions"]:
                continue
            pos_data = interval_dict["annotations"]["positions"][pos_str]
            if range_id not in pos_data:
                continue
            merge_nested_data(pos_data[range_id], aggregated, pos_str)

    elif data_type == "conservations":
        for pos in range(start, end + 1):
            pos_str = str(pos)
            if pos_str in interval_dict["conservations"]["positions"]:
                cons_pos_data = interval_dict["conservations"]["positions"][pos_str]
                # Track per-position conservation, hit and residue data
                if "conservation" not in aggregated:
                    aggregated["conservation"] = {}
                if "hit" not in aggregated:
                    aggregated["hit"] = {}
                if "residue" not in aggregated:
                    aggregated["residue"] = {}

                aggregated["conservation"][pos_str] = cons_pos_data["conservation"]
                aggregated["hit"][pos_str] = cons_pos_data["hit"]
                aggregated["residue"][pos_str] = cons_pos_data["residue"]

    return dict(aggregated), metadata

def transform_to_ranges(interval_dict: Dict, multi_logger: Callable) -> Dict:
    """Transform position-based data to range-based format."""
    range_based = defaultdict(lambda: defaultdict(dict))

    if not interval_dict:
        return {}

    for data_type, ranges_field in [
        ("annotations", "annotation_ranges"),
        ("conservations", "conservation_ranges")
    ]:
        ranges_dict = interval_dict.get(ranges_field, {})

        for range_id, ranges_data in ranges_dict.items():
            ranges_list = ranges_data.get("ranges", [])

            if not isinstance(ranges_list, list):
                multi_logger("warning", "MAKE_VIEW - TRANSFORM_2_RANGES - %s invalid format for %s", ranges_field, range_id)
                continue

            try:
                # Process each range in the ranges list
                for start, end in ranges_list:
                    range_key = f"{start}-{end}"

                    # Get both position data and metadata
                    position_data, metadata = aggregate_range_positions(
                        interval_dict, range_id, (start, end), data_type
                    )

                    # Store position data under range key
                    range_based[range_id][range_key] = position_data

                    # Add metadata at the same level as range_id
                    if metadata:
                        range_based[data_type] = metadata

            except (TypeError, ValueError) as e:
                multi_logger("warning", "MAKE_VIEW - TRANSFORM_2_RANGES - Range processing error for %s: %s", range_id, e)
                continue

    return dict(range_based)

def process_sequence_report(report_path: str, sequence_id: str, logger: logging.Logger, multi_logger: Callable) -> Dict:
    """Process the aggregated report file for a specific sequence.

    Args:
        report_path: Path to aggregated_report.json
        sequence_id: ID of sequence to process from aggregated report
        logger: Logger instance
        multi_logger: Multi-logger callable
    """
    logger.info("MAKE_VIEW - PROC_SEQ_REP - Starting to process report: %s", report_path)

    try:
        with open(report_path, 'r', encoding='utf-8') as f:
            # Convert lists to original types for processing
            data = convert_lists_to_original_types(json.load(f))
    except (FileNotFoundError, IOError) as e:
        multi_logger("error", "MAKE_VIEW - PROC_SEQ_REP - Failed to open report file %s: %s", report_path, e)
        raise

    if not data:
        multi_logger("warning", "MAKE_VIEW - PROC_SEQ_REP - Empty report for sequence %s - no domains were found", sequence_id)
        return None

    if sequence_id not in data:
        multi_logger("error", "MAKE_VIEW - PROC_SEQ_REP - Sequence %s not found in report", sequence_id)
        raise ValueError(f"Sequence {sequence_id} not found in aggregated report")

    transformed = {
        "sequence_id": sequence_id,
        "domain": {}
    }

    sequence_data = data[sequence_id]
    for domain_id, domain_data in sequence_data.items():
        transformed["domain"][domain_id] = {}

        for interval_key, interval_data in domain_data.get("hit_intervals", {}).items():
            transformed["domain"][domain_id][interval_key] = transform_to_ranges(interval_data, multi_logger)
            logger.info(
                "MAKE_VIEW - PROC_SEQ_REP - Processed %s - %s - %s",
                sequence_id, domain_id, interval_key
            )
            transformed_data = transformed["domain"][domain_id][interval_key]
            transformed_data_to_list = convert_sets_and_tuples_to_lists(transformed_data)
            # Convert defaultdict to dict for logging
            transformed_data_dict = convert_defaultdict_to_dict(transformed_data_to_list)
            # Log the transformed data as a JSON-formatted string
            logger.info(
                "MAKE_VIEW - PROC_SEQ_REP - Transformed data: %s",
                json.dumps(transformed_data_dict, indent=2)
            )

    return transformed

def write_range_views(transformed_data: Dict, sequence_dir: str, clean_sequence_id: str, logger: logging.Logger, multi_logger: Callable) -> None:
    """Write range-based views to files."""

    try:
        logger.info("MAKE VIEW - Writing range views for sequence %s", clean_sequence_id)
        print(transformed_data)
        for domain_id in transformed_data["domain"]:
            domain_data = {
                "sequence_id": transformed_data["sequence_id"],
                "domain": {
                    domain_id: transformed_data["domain"][domain_id]
                }
            }
            output_path = os.path.join(sequence_dir, f"{domain_id}_ranges.json")
            logger.debug("MAKE VIEW - Writing data to: %s", output_path)
            cleaned_domain_data = convert_sets_and_tuples_to_lists(domain_data)
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(cleaned_domain_data, f, indent=2)
            logger.info("MAKE VIEW - Successfully wrote range view for %s - %s", clean_sequence_id, domain_id)

    except (PermissionError, OSError) as e:
        multi_logger("error", "MAKE VIEW - Failed to write range view for %s: %s", clean_sequence_id, e)
        raise

def main():
    """Main execution function."""
    args = parse_arguments()
    main_logger, _ = get_logger(args.log, scope="main")
    clean_sequence_id = args.sequence
    sequence_id = clean_sequence_id.replace("-", "|")
    sequence_logger, _ = get_logger(args.log, scope="sequence", identifier=clean_sequence_id)
    sequence_dir = args.sequence_dir
    sequence_logger.info("MAKE VIEW - Will process aggregated_report.json in %s", sequence_dir)

    log_to_both = get_multi_logger([main_logger, sequence_logger])

    aggregated_report = os.path.join(sequence_dir, "aggregated_report.json")
    if not os.path.exists(aggregated_report):
        log_to_both("error", "MAKE VIEW - No aggregated_report.json found at %s", aggregated_report)
        sys.exit(1)

    try:
        transformed_data = process_sequence_report(
            aggregated_report,
            sequence_id,
            sequence_logger,
            log_to_both
        )

        if transformed_data is None:
            log_to_both("info", "MAKE_VIEW - No domains were found for sequence %s - skipping view generation", clean_sequence_id)
            sys.exit(0)  # Exit cleanly, this is not an error case

        write_range_views(transformed_data, sequence_dir, clean_sequence_id, sequence_logger, log_to_both)
    except (IOError, json.JSONDecodeError, KeyError) as e:
        log_to_both("error", "Error processing aggregated report: %s", str(e))
        sys.exit(1)

if __name__ == '__main__':
    main()