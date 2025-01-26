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
4. write_range_views: Outputs processed data in Nightingale-compatible format

Input: Per-sequence JSON reports from transfer_annotations.py
Output: Range-based JSON views for Nightingale visualization
"""


import os
import json
import logging
import argparse
from collections import defaultdict
from typing import Tuple, Dict
from utils import convert_lists_to_original_types, convert_sets_and_tuples_to_lists

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Transform position-based annotations into range-based views")
    parser.add_argument("-o", "--output-dir", help="Path to output directory", required=True)
    parser.add_argument("-l", "--log", help="Log file path", default="logs/make_views_jsons.log")
    return parser.parse_args()

def configure_logging(log_path: str) -> logging.Logger:
    """Configure logging."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger()

def aggregate_range_annotations(interval_dict: Dict, anno_id: str, range_tuple: Tuple[int, int]) -> Dict:
    """Aggregate annotations for a specific range."""
    start, end = range_tuple
    aggregated = defaultdict(lambda: defaultdict(dict))

    for pos in range(start, end + 1):
        pos_str = str(pos)
        if pos_str not in interval_dict["annotations"]["positions"]:
            continue

        pos_data = interval_dict["annotations"]["positions"][pos_str]
        if anno_id not in pos_data:
            continue

        anno_data = pos_data[anno_id]

        # Aggregate essentials
        if "essentials" not in aggregated:
            aggregated["essentials"] = anno_data["essentials"]
        else:
            aggregated["essentials"]["count"] += anno_data["essentials"]["count"]

        # Aggregate evidence
        for ev_type, ev_data in anno_data.get("evidence", {}).items():
            if ev_type not in aggregated["evidence"]:
                aggregated["evidence"][ev_type] = ev_data.copy()
            else:
                aggregated["evidence"][ev_type]["count"] += ev_data["count"]

        # Copy hit status
        aggregated["hit"] = anno_data.get("hit", False)

    return dict(aggregated)

def transform_to_ranges(interval_dict: Dict) -> Dict:
    """Transform position-based annotations to range-based format."""
    range_based = defaultdict(dict)

    for anno_id, ranges_data in interval_dict.get("annotation_ranges", {}).items():
        for start, end in ranges_data.get("ranges", []):
            range_key = f"{start}-{end}"
            range_based[anno_id][range_key] = aggregate_range_annotations(
                interval_dict, anno_id, (start, end)
            )

    return dict(range_based)

def process_sequence_report(report_path: str, logger: logging.Logger = None) -> Dict:
    """Process a single sequence report file."""
    with open(report_path, 'r', encoding='utf-8') as f:
        # Convert lists to original types for processing
        data = convert_lists_to_original_types(json.load(f))

    sequence_id = data["sequence_id"]
    transformed = {"sequence_id": sequence_id, "domain": {}}

    for domain_id, domain_data in data["domain"].items():
        transformed["domain"][domain_id] = {}

        for interval_key, interval_dict in domain_data.get("hit_intervals", {}).items():
            transformed["domain"][domain_id][interval_key] = transform_to_ranges(interval_dict)
            if logger:
                logger.info(f"Processed {sequence_id} - {domain_id} - {interval_key}")
                logger.info(f"Transformed data: {transformed['domain'][domain_id][interval_key]}")

    return transformed

def write_range_views(transformed_data: Dict, output_dir: str, logger: logging.Logger) -> None:
    """Write range-based views to files."""
    sequence_id = transformed_data["sequence_id"]
    sequence_dir = os.path.join(output_dir, sequence_id)
    os.makedirs(sequence_dir, exist_ok=True)

    for domain_id in transformed_data["domain"]:
        output_path = os.path.join(sequence_dir, f"{domain_id}_ranges.json")

        # Convert sets and tuples to lists for JSON serialization
        cleaned_transformed_data = convert_sets_and_tuples_to_lists(transformed_data)

        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(cleaned_transformed_data, f, indent=2)

        logger.info(f"Wrote range view for {sequence_id} - {domain_id}")

def main():
    """Main execution function."""
    args = parse_arguments()
    logger = configure_logging(args.log)

    try:
        # Process all sequence directories
        for entry in os.scandir(args.output_dir):
            if entry.is_dir() and not entry.name.startswith('PF'):
                for file in os.scandir(entry.path):
                    if file.name.endswith('_report.json'):
                        logger.info("Processing %s", file.path)
                        try:
                            transformed_data = process_sequence_report(file.path, logger)
                            write_range_views(transformed_data, args.output_dir, logger)
                        except (IOError, json.JSONDecodeError) as e:
                            logger.error("Error processing %s: %s", file.path, str(e))
                            continue

        # Create done file
        with open(os.path.join(args.output_dir, "make_views_jsons.done"), 'w', encoding='utf-8') as f:
            f.write('')

    except Exception as e:
        logger.error("Error in main execution: %s", str(e))
        raise

if __name__ == '__n__':
    main()