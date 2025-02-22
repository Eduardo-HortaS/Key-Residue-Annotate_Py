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
import json
import logging
import argparse
from collections import defaultdict
from typing import Tuple, Dict, Callable
from utils import convert_lists_to_original_types, convert_sets_and_tuples_to_lists, convert_defaultdict_to_dict, get_logger, get_multi_logger

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Transform position-based annotations into range-based views")
    parser.add_argument("-o", "--output-dir", help="Path to output directory where sequence subdirs are", required=True)
    parser.add_argument("-s", "--sequence", help="Sequence to process", required=True)
    parser.add_argument("-l", "--log", help="Log file path", default="logs/make_views_jsons.log")
    return parser.parse_args()

def track_position_data(target: Dict, key: str, pos_str: str, value: Dict, track_key: str) -> None:
    """Helper to track per-position data like counts and hits."""
    if track_key not in target[key]:
        target[key][track_key] = {}
    target[key][track_key][pos_str] = value[track_key]

def merge_nested_data(source: Dict, target: Dict, pos_str: str, range_id: str) -> None:
    """Recursively merge nested dictionaries at range_id, tracking counts and hits per position."""
    if "hit" in source:
        if "hit" not in target:
            target["hit"] = source["hit"]

    for key, value in source.items():
        if isinstance(value, dict):
            if key not in target:
                target[key] = defaultdict(dict)

            # Track position-specific data
            if "count" in value:
                track_position_data(target, key, pos_str, value, "count")
            if "hit" in value:
                track_position_data(target, key, pos_str, value, "hit")

            if key == "GO":
                if pos_str not in target["GO"]:
                    target["GO"][pos_str] = {}
                target["GO"][pos_str].update(value)
                continue

            # Recurse for nested structures
            merge_nested_data(value, target[key], pos_str, range_id)
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
            merge_nested_data(pos_data[range_id], aggregated, pos_str, range_id)

    elif data_type == "conservations":
        conservation_data = defaultdict(list)
        for pos in range(start, end + 1):
            pos_str = str(pos)
            if pos_str in interval_dict["conservations"]["positions"]:
                cons_pos_data = interval_dict["conservations"]["positions"][pos_str]
                conservation_data["conservation"].append(cons_pos_data["conservation"])
                conservation_data["hit"].append(cons_pos_data["hit"])

        if conservation_data:
            # Aggregate conservation data for range
            aggregated["conservation"] = sum(conservation_data["conservation"]) / len(conservation_data["conservation"])
            aggregated["hit"] = all(conservation_data["hit"])  # True if all positions hit

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

def process_sequence_report(report_path: str, logger: logging.Logger, multi_logger: Callable) -> Dict:
    """Process a single sequence report file."""
    if logger:
        logger.info("MAKE_VIEW - PROC_SEQ_REP - Starting to process report: %s", report_path)
    try:
        with open(report_path, 'r', encoding='utf-8') as f:
            # Convert lists to original types for processing
            data = convert_lists_to_original_types(json.load(f))
    except (FileNotFoundError, IOError) as e:
        multi_logger("error", "MAKE_VIEW - PROC_SEQ_REP - Failed to open report file %s: %s", report_path, e)
        raise

    sequence_id = data["sequence_id"]
    transformed = {"sequence_id": sequence_id, "domain": {}}

    for domain_id, domain_data in data["domain"].items():
        transformed["domain"][domain_id] = {}

        for interval_key, interval_data in domain_data.get("hit_intervals", {}).items():
            transformed["domain"][domain_id][interval_key] = transform_to_ranges(interval_data, multi_logger)
            if logger:
                logger.info("MAKE_VIEW - PROC_SEQ_REP - Processed %s - %s - %s", sequence_id, domain_id, interval_key)
                transformed_data = transformed["domain"][domain_id][interval_key]
                transformed_data_to_list = convert_sets_and_tuples_to_lists(transformed_data)
                # Convert defaultdict to dict for logging
                transformed_data_dict = convert_defaultdict_to_dict(transformed_data_to_list)
                # Log the transformed data as a JSON-formatted string
                logger.info("MAKE_VIEW - PROC_SEQ_REP - Transformed data: %s", json.dumps(transformed_data_dict, indent=2))


    return transformed

def write_range_views(transformed_data: Dict, output_dir: str, logger: logging.Logger, multi_logger: Callable) -> None:
    """Write range-based views to files."""
    clean_sequence_id = transformed_data["sequence_id"].replace("|", "-")
    try:
        logger.info("MAKE VIEW - Writing range views for sequence %s", clean_sequence_id)
        sequence_dir = os.path.join(output_dir, clean_sequence_id)
        os.makedirs(sequence_dir, exist_ok=True)
        logger.debug("MAKE VIEW - Created or verified directory: %s", sequence_dir)

        for domain_id in transformed_data["domain"]:
            output_path = os.path.join(sequence_dir, f"{domain_id}_ranges.json")
            logger.debug("MAKE VIEW - Writing data to: %s", output_path)
            cleaned_transformed_data = convert_sets_and_tuples_to_lists(transformed_data)
            with open(output_path, 'w', encoding='utf-8') as f:
                json.dump(cleaned_transformed_data, f, indent=2)
            logger.info("MAKE VIEW - Successfully wrote range view for %s - %s", clean_sequence_id, domain_id)

    except (PermissionError, OSError) as e:
        multi_logger("error", "MAKE VIEW - Failed to write range view for %s: %s", clean_sequence_id, e)
        raise

def main():
    """Main execution function."""
    args = parse_arguments()
    main_logger, _ = get_logger(args.log, scope="main")
    sequence_logger, _ = get_logger(args.log, scope="sequence", identifier=args.sequence)
    output_dir = args.output_dir
    sequence_logger.info("MAKE VIEW - Starting to process sequence reports in %s", output_dir)

    log_to_both = get_multi_logger([main_logger, sequence_logger])

    for file in os.scandir(output_dir):
        if file.name.endswith('_report.json'):
            try:
                transformed_data = process_sequence_report(file.path, sequence_logger, log_to_both)
                write_range_views(transformed_data, output_dir, sequence_logger, log_to_both)
            except (IOError, json.JSONDecodeError) as e:
                log_to_both("error", "Error processing %s: %s", file.path, str(e))
                continue

    # Create done file
    with open(os.path.join(output_dir, "make_views_jsons.done"), 'w', encoding='utf-8') as f:
        f.write('')

if __name__ == '__main__':
    main()