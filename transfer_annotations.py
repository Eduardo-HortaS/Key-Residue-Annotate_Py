"""
transfer_annotations.py

Copyright 2025 Eduardo Horta Santos <GitHub: Eduardo-HortaS>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any pairedr version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.

This script contains all functions dealing directly with transferring annotations.

Intended usage:

Runs for a single domain with 3 required arguments:
paths to the domain's hmmalign result and resource and output dirs.
Optionally takes a list of desired ECO codes to filter annotations.
Execution begins by calling find_and_map_annots with 2 arguments:
a list of hmmalign result lines and the loaded annotations dict.

The function call order is as follows:
parse_arguments -> configure_logging -> get_pfam_id_from_hmmalign_result -> get_annotation_filepath -> read_files ->
find_and_map_annots -> map_and_filter_annot_pos -> validate_annotations -> process_annotations -> make_anno_total_dict
(Anno_Total is None?) YES -> DO NOTHING --- position effectively skipped
(Anno_Total is None?) NO -> PROCEED CHECK PAIREABLE ANNOTATIONS
(paireable annotation?) NO -> add_to_transfer_dict
(paireable annotation?) YES -> map_and_filter_annot_pos -> validate_paired_annotations ->  make_anno_total_dict -> process_annotation -> add_to_transfer_dict
"""

import os
import json
import re
import argparse
import logging
import traceback
import copy
import pandas as pd
from typing import Any, Dict, Optional, Generator, Tuple
# from modules.decorators import measure_time_and_memory
# from memory_profiler import profile


def parse_arguments():
    """Parse command-line arguments for transferring annotations
    per domain using its hmmalign alignment.
    Generates PF*_report.json for each domain found, in every sequence.
    Also needs paths to output and resource directories.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description=
    'Generates a temporary multifasta for running hmmalign using a hits per domain JSON.')
    parser.add_argument("-iA", "--dom-align", required=True, type=str, help="Path to domain's hmmalign alignment")
    parser.add_argument("-r", "--resource-dir", required=True, type=str, help="Resource dir path")
    parser.add_argument("-o", "--output-dir", required=True, type=str, help="Output dir path")
    parser.add_argument("-e", "--eco-codes", required=False, default=[], nargs="*", help="Space-separated ECO codes to filter annotations")
    parser.add_argument("-l", "--log", required=False, default="logs/transfer_annotations.log", type=str, help="Log path")

    args = parser.parse_args()

    # Clean eco codes - remove any quotes, brackets and commas
    if args.eco_codes:
        args.eco_codes = [code.strip('",[]') for code in args.eco_codes]

    return args

def configure_logging(log_path: str) -> logging.Logger:
    """Set up logging for the script."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(filename=log_path, level=logging.DEBUG, filemode='a', \
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    return logging.getLogger()

def get_pfam_id_from_hmmalign_result(hmmalign_result: str) -> str:
    """
    Extracts the PFAM ID from the hmmalign result filename.
    """
    return os.path.splitext(os.path.basename(hmmalign_result))[0].split('_')[-2]

def get_annotation_filepath(resource_dir: str, pfam_id: str) -> str:
    """
    Constructs and returns the filepaths for annotations JSON based on the PFAM ID.
    """
    base_dir = os.path.join(resource_dir, pfam_id)
    annotations_filepath = os.path.join(base_dir, "annotations.json")
    conservations_filepath = os.path.join(base_dir, "conservations.json")
    return annotations_filepath, conservations_filepath

def read_files(hmmalign_result: str, annotations_filepath: str) -> tuple[list[str], dict]:
    """
    Reads and returns the content of the hmmalign result and annotations files,
    respectively, as lists of lines and a loaded JSON object.
    """
    with open(hmmalign_result, 'r', encoding="utf-8") as hmmaligned_file:
        hmmalign_lines = [line.rstrip('\n') for line in hmmaligned_file]
        # # # # hmmalign_lines = hmmaligned_file.readlines()
    with open(annotations_filepath, 'r', encoding="utf-8") as annotations_file:
        annotations = json.load(annotations_file)
    return hmmalign_lines, annotations

### New Helper Functions
def iterate_aligned_sequences(
    seq_a: str,
    seq_b: str,
    start_a: int,
    start_b: int,
    end_a: int,
    end_b: int
) -> Generator[Tuple[int, Optional[int], Optional[int], str, str], None, None]:
    """
    Yields index, current positions (pos_a, pos_b), and characters (char_a, char_b).
    Positions are incremented only when the respective character is alphabetic.
    """
    pos_a = None
    pos_b = None
    for index, (char_a, char_b) in enumerate(zip(seq_a, seq_b)):
        if char_a.isalpha():
            pos_a = start_a if pos_a is None else pos_a + 1
        if char_b.isalpha():
            pos_b = start_b if pos_b is None else pos_b + 1

        yield index, pos_a, pos_b, char_a, char_b

        if (pos_a is not None and pos_a == end_a) or (pos_b is not None and pos_b == end_b):
            break

#@measure_time_and_memory
#@profile
def find_and_map_annots(
    logger: logging.Logger,
    hmmalign_lines: list,
    annotations: dict,
    good_eco_codes: list) -> dict:
    """
    From the hmmalign alignment lines, finds target lines and stores data for them in target_info,
    organized by target name and hit interval. Then, finds seed sequences in hmmalign alignment
    for which we have annotations and calls map_and_filter_annot_pos to map positions
    by iteratin on both target and seed sequences. Makes processed_annotations set
    to keep track of visited positions across calls, where relevant, and sends it to map_and_filter_annot_pos.

    Args:
        logger (logging.Logger): Logger object for logging.
        hmmalign_lines (list): List of lines from the hmmalign alignment file.
        annotations (dict): Loaded annotations JSON object.
        good_eco_codes (list): List of ECO codes to filter annotations.

    Returns:
        dict: Transfer dictionary containing annotations for each sequence-domain pair.
    """
    transfer_dict = {}
    target_info = {}
    processed_annotations = set()

    for line in hmmalign_lines:
        if "target/" in line and not line.startswith("#"):
            logger.debug(f"---> DEBUG --- FIND_AND_MAP --- Target Line: {line}")
            parts = line.split("target/")
            target_name = parts[0].strip()
            target_relevant = parts[1].strip()
            target_parts = target_relevant.split()
            if len(target_parts) >= 2:
                target_hit_interval = target_parts[0].split("/")[1]
                target_hit_start = int(target_hit_interval.split("-")[0])
                target_hit_end = int(target_hit_interval.split("-")[1])
                target_hit_sequence = target_parts[1]
                if target_name not in target_info:
                    target_info[target_name] = {}
                if target_hit_interval not in target_info[target_name]:
                    target_info[target_name][target_hit_interval] = []
                target_info[target_name][target_hit_interval].append((target_hit_start, target_hit_end, target_hit_sequence))

    if not target_info:
        logger.error("---> ERROR --- FIND_AND_MAP --- No target sequences found in hmmalign lines!")
        return transfer_dict

    for entry_mnemo_name, entry_annotations in annotations.items():
        for line in hmmalign_lines:
            if line.startswith("#") or 'target/' in line:
                continue
            splitted = line.split('/')
            entry_mnemo_name_in_line = splitted[0]
            if entry_mnemo_name == entry_mnemo_name_in_line:
                offset_start = int(splitted[1].split('-')[0])
                offset_end = int(splitted[1].split('-')[1].split()[0])
                annot_sequence = splitted[1].split()[1]
                for target_name, target_per_interval in target_info.items():
                    for target_hit_interval, target_info_list in target_per_interval.items():
                        for target_hit_start, target_hit_end, target_hit_sequence in target_info_list:
                            entry_annotations_copy = copy.deepcopy(entry_annotations)
                            # # DEBUG, Delete Later
                            # if entry_mnemo_name == "Q5CUW0_CRYPI":
                            map_and_filter_annot_pos(
                                logger=logger,
                                good_eco_codes=good_eco_codes,
                                target_sequence=target_hit_sequence,
                                target_name=target_name,
                                target_hit_start=target_hit_start,
                                target_hit_end=target_hit_end,
                                offset_start=offset_start,
                                offset_end=offset_end,
                                annot_sequence=annot_sequence,
                                entry_mnemo_name=entry_mnemo_name,
                                entry_annotations=entry_annotations_copy,
                                transfer_dict=transfer_dict,
                                processed_annotations=processed_annotations
                            )
                            logger.debug(f"---> DEBUG --- FIND_AND_MAP --- Mapped, filtered and possibly added to transfer_dict: target {target_name} and annotated {entry_mnemo_name} at target hit interval {target_hit_interval}")
    return transfer_dict

def read_conservations_and_annotations(conservations_filepath: str, annotations_filepath: str) -> tuple[dict, dict]:
    """
    Reads the conservations and annotations JSON files.
    Returns appropriate empty structures if files are empty.

    Returns:
        tuple[dict, dict]: (conservations, annotations) where:
            - conservations: {"sequence_id/start-end": {"position": conservation_score, ...}}
            - annotations: {"entry_name": {"position": [{"type": str, ...}]}}
    """
    with open(conservations_filepath, 'r', encoding="utf-8") as conservations_file:
        conservations = json.load(conservations_file)
        if not conservations:
            conservations = {"sequence_id/range": {}}

    with open(annotations_filepath, 'r', encoding="utf-8") as annotations_file:
        annotations = json.load(annotations_file)
        if not annotations:
            annotations = {"entry_name": {}}

    return conservations, annotations

def parse_go_annotations(go_column: str) -> list:
    """
    Parses the GO Annotations column from the InterProScan TSV output.

    Args:
        go_column (str): The GO Annotations column value (e.g., "GO:0005515|GO:0006302").

    Returns:
        list: A list of GO terms (e.g., ["GO:0005515", "GO:0006302"]).
    """
    if not go_column or go_column == "-" or not go_column.strip():
        return []
    return [term.split("(")[0].strip() for term in go_column.split("|") if term]

def gather_go_terms_for_target(target_name: str, go_terms_dir: str) -> set:
    """
    Gathers GO terms for a single target from iprscan.tsv if present.
    Note: GO terms are kept as a string of GO terms separated by "|".
    The GO terms base dir is expected to be the same as the output dir,
    since this function requires the previous scripts to have been run.
    """
    go_term_filepath = os.path.join(go_terms_dir, target_name, "iprscan.tsv")
    target_gos = set()
    if os.path.exists(go_term_filepath):
        with open(go_term_filepath, 'r', encoding="utf-8") as go_term_file:
            iprscan_df = pd.read_csv(go_term_file, sep='\t', header=None)
            iprscan_df.columns = [
                "Protein Accession", "Sequence MD5 digest", "Sequence Length", "Analysis",
                "Signature Accession", "Signature Description", "Start Location", "Stop Location",
                "Score", "Status", "Date", "InterPro Accession", "InterPro Description",
                "GO Annotations", "Pathways Annotations"
            ]
            target_go_annots = iprscan_df["GO Annotations"].dropna().tolist()
            for go_entry in target_go_annots:
                target_gos.update(parse_go_annotations(go_entry))
    return target_gos


def get_alignment_sequences(hmmalign_lines: list, target_id: str, conservation_id: str):
    """
    From hmmalign_lines, extracts the target_seq and conservation_seq matching the given IDs.
    Returns (target_seq, conservation_seq).
    """
    target_seq = None
    conservation_seq = None

    for line in hmmalign_lines:
        if line.startswith(target_id) or line.startswith(conservation_id):
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                sequence = parts[1]
                if seq_id == conservation_id:
                    conservation_seq = sequence
                elif seq_id == target_id:
                    target_seq = sequence
            if target_seq and conservation_seq:
                break

    if target_seq is None:
        raise ValueError(f"Target ID '{target_id}' not found in hmmalign_lines.")
    if conservation_seq is None:
        raise ValueError(f"Conservation ID '{conservation_id}' not found in hmmalign_lines.")

    return target_seq, conservation_seq


def populate_conservation(
    transfer_dict: dict,
    pfam_id: str,
    target_name: str,
    target_seq: str,
    conservation_seq: str,
    conservation_key: str,
    conservations: dict,
    conserved_positions: list,
    target_hit_start: int,
    target_hit_end: int,
    conservation_start: int,
    conservation_end: int,
    logger: Optional[logging.Logger] = None
) -> None:
    """
    Iterates over aligned sequences and populates transfer_dict with conservation info.
    """
    processed_conserved_pos = set()
    for index, counter_cons_pos, counter_target_pos, char_cons, char_target in iterate_aligned_sequences(
        seq_a=conservation_seq,
        seq_b=target_seq,
        start_a=conservation_start,
        start_b=target_hit_start,
        end_a=conservation_end,
        end_b=target_hit_end
    ):
        if counter_target_pos is None or counter_cons_pos is None:
            continue

        counter_cons_pos_str = str(counter_cons_pos)
        if counter_cons_pos_str in conserved_positions:
            processed_conserved_pos.add(counter_cons_pos_str)
            counter_target_pos_str = str(counter_target_pos)
            index_str = str(index)
            is_match = char_cons == char_target and char_cons not in '.-' and char_target not in '.-'
            cons_type = "matches" if is_match else "misses"

            score_data = conservations[conservation_key][counter_cons_pos_str]
            counter_target_pos_info = {
                "conservation": score_data,
                "hit": is_match
            }

            transfer_dict[pfam_id]['sequence_id'][target_name]["conservations"]["positions"][counter_target_pos_str] = counter_target_pos_info
            transfer_dict[pfam_id]['sequence_id'][target_name]["conservations"]["indices"][cons_type].add(counter_target_pos_str)

            # Add position conversion info
            transfer_dict[pfam_id]['sequence_id'][target_name]['position_conversion']['target_to_aln'][counter_target_pos_str] = index_str
            aln_dict = transfer_dict[pfam_id]['sequence_id'][target_name]['position_conversion']['aln_to_target']
            if index_str not in aln_dict:
                aln_dict[index_str] = set()
            aln_dict[index_str].add(counter_target_pos_str)

        missing = set(conserved_positions) - processed_conserved_pos
        if missing and logger:
            logger.warning(f"Missing conserved positions: {missing}")
        if counter_target_pos == target_hit_end or counter_cons_pos == conservation_end:
            break


def populate_go_data_for_annotations(
    transfer_dict: dict,
    pfam_id: str,
    target_name: str,
    annotations: dict,
    go_terms_annot_key: str,
    target_gos: set
) -> None:
    """
    Fills GO data into the transfer_dict for each annotation found, using any overlapping GO terms.
    """
    go_info_cache = {}
    for pos, position_data in transfer_dict[pfam_id]['sequence_id'][target_name]['annotations']['positions'].items():
        for anno_id, anno_data in position_data.items():
            go_data_for_annotation = {}
            for _, evidence_data in anno_data.get("evidence", {}).items():
                annot_name = evidence_data.get("rep_mnemo_name")
                if annot_name and annot_name in annotations and go_terms_annot_key in annotations[annot_name]:
                    cache_key = (target_name, annot_name)
                    if cache_key not in go_info_cache:
                        matched_go_info = {}
                        annotation_go_set = set()

                        for subontology, go_terms_with_meanings in annotations[annot_name][go_terms_annot_key].items():
                            for go_term, meaning in go_terms_with_meanings.items():
                                annotation_go_set.add(go_term)
                                if go_term in target_gos:
                                    matched_go_info[go_term] = meaning

                        target_go_set = set(target_gos)
                        intersection = len(target_go_set.intersection(annotation_go_set))
                        union = len(target_go_set.union(annotation_go_set))
                        jaccard_index = intersection / union if union > 0 else 0

                        go_info_cache[cache_key] = {
                            annot_name: {
                                "terms": matched_go_info,
                                "jaccard_index": jaccard_index
                            }
                        }
                    go_data_for_annotation.setdefault("GO", {}).update(go_info_cache[cache_key])
            if go_data_for_annotation:
                anno_data.update(go_data_for_annotation)


def cleanup_improve_transfer_dict(
    transfer_dict: dict,
    pfam_id: str,
    hmmalign_lines: list,
    conservations_filepath: str,
    annotations_filepath: str,
    output_dir: str,
) -> dict:
    """
    Cleans up and improves the transfer dictionary by:
    - Setting transfer_dict[pfam_id] from transfer_dict["DOMAIN"]
    - Adding conservation data
    - Adding GO term data
    """

    go_terms_annot_key = "0"
    # print("Ayooooooo - Transfer dict")
    # debug_dict = convert_sets_to_lists(transfer_dict)
    # print(json.dumps(debug_dict, indent=4))
    transfer_dict[pfam_id] = transfer_dict.pop("DOMAIN")

    # Read external JSON files
    conservations, annotations = read_conservations_and_annotations(conservations_filepath, annotations_filepath)

    # For each target in the dictionary, gather GO terms, get aligned sequences and fill conservation and GO data
    for target_name in transfer_dict[pfam_id]['sequence_id']:
        target_gos = gather_go_terms_for_target(target_name, output_dir)
        target_hit_start = int(transfer_dict[pfam_id]['sequence_id'][target_name]['hit_start'])
        target_hit_end = int(transfer_dict[pfam_id]['sequence_id'][target_name]['hit_end'])
        target_id = f"{target_name}target//{target_hit_start}-{target_hit_end}"

        # There's only one key in conservations, a sequence from the
        # original alignment used as a reference for the
        # conserved positions between the original and query/new alignments.
        conservation_key = list(conservations.keys())[0]
        conserved_positions = list(conservations[conservation_key].keys())

        # Parse the numeric range from the conservation key
        conservation_range = conservation_key.split("/")[1]
        conservation_start, conservation_end = map(int, conservation_range.split("-"))

        # Extract sequences from hmmalign
        target_seq, conservation_seq = get_alignment_sequences(hmmalign_lines, target_id, conservation_key)

        # Populate the conservation data if sequences were found
        if target_seq and conservation_seq:
            populate_conservation(
                transfer_dict=transfer_dict,
                pfam_id=pfam_id,
                target_name=target_name,
                target_seq=target_seq,
                conservation_seq=conservation_seq,
                conservation_key=conservation_key,
                conservations=conservations,
                conserved_positions=conserved_positions,
                target_hit_start=target_hit_start,
                target_hit_end=target_hit_end,
                conservation_start=conservation_start,
                conservation_end=conservation_end,
            )

        # Populate GO data for each annotation
        populate_go_data_for_annotations(
            transfer_dict=transfer_dict,
            pfam_id=pfam_id,
            target_name=target_name,
            annotations=annotations,
            go_terms_annot_key=go_terms_annot_key,
            target_gos=target_gos
        )

    return transfer_dict

def convert_sets_to_lists(data):
    if isinstance(data, dict):
        return {k: convert_sets_to_lists(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [convert_sets_to_lists(elem) for elem in data]
    elif isinstance(data, set):
        return sorted(list(data))
    else:
        return data


def write_reports(
    logger: logging.Logger,
    transfer_dict: dict,
    output_dir: str,
) -> None:
    """
    1. Writes the entire transfer_dict to a single JSON file under:
       output_dir/pfam_id/pfam_id_report.json
    2. Writes each target_name's data under:
       output_dir/target_name/pfam_id_report.json, preserving pfam_id structure
    """
    transfer_dict = convert_sets_to_lists(transfer_dict)
    pfam_id = next(iter(transfer_dict.keys()))

    # Write the entire transfer_dict to the pfam_id subdir
    pfam_dir = os.path.join(output_dir, pfam_id)
    os.makedirs(pfam_dir, exist_ok=True)
    entire_report_path = os.path.join(pfam_dir, pfam_id + "_report.json")

    structured_report = {
        "domain_id": pfam_id,
        "sequences": transfer_dict[pfam_id]['sequence_id']
    }

    with open(entire_report_path, 'w', encoding="utf-8") as report_file:
        json.dump(structured_report, report_file, indent=4)
    logger.debug(f"---> DEBUG --- WRITE_REPORT --- Wrote entire Transfer Report: {entire_report_path}")

    # Write individual files for each target_name, preserving pfam_id structure
    for target_name, target_data in transfer_dict[pfam_id]['sequence_id'].items():
        if target_data:
            sequence_dict = {
                "sequence_id": target_name,
                "domain": {
                    pfam_id: target_data
                }
            }
            logger.debug(f"---> DEBUG --- WRITE_REPORT --- Data to write for {target_name}!")
        else:
            sequence_dict = {
                "sequence_id": target_name,
                "domain": {
                    pfam_id: {"annotations": "None"}
                }
            }
            logger.debug(f"---> DEBUG --- WRITE_REPORT --- No data to write for {target_name}!")

        safe_target_name = target_name.replace("|", "-")
        target_report_filepath = os.path.join(output_dir, safe_target_name, pfam_id + "_report.json")
        os.makedirs(os.path.join(output_dir, safe_target_name), exist_ok=True)

        with open(target_report_filepath, 'w', encoding="utf-8") as report_file:
            json.dump(sequence_dict, report_file, indent=4)

        logger.debug(f"---> DEBUG --- WRITE_REPORT --- Wrote Transfer Report for {target_name}-{pfam_id}: {target_report_filepath}")


# def write_reports(
#     logger: logging.Logger,
#     transfer_dict: dict,
#     output_dir: str,
# ) -> None:
#     """
#     1. Writes the entire transfer_dict to a single JSON file under:
#        output_dir/pfam_id/pfam_id_report.json
#     2. Optionally writes each target_name's data under:
#        output_dir/target_name/pfam_id_report.json
#     """
#     transfer_dict = convert_sets_to_lists(transfer_dict)
#     pfam_id = next(iter(transfer_dict.keys()))

#     # Write the entire transfer_dict to the pfam_id subdir
#     pfam_dir = os.path.join(output_dir, pfam_id)
#     os.makedirs(pfam_dir, exist_ok=True)
#     entire_report_path = os.path.join(pfam_dir, pfam_id + "_report.json")
#     with open(entire_report_path, 'w', encoding="utf-8") as report_file:
#         json.dump(transfer_dict, report_file, indent=4)
#     logger.debug(f"---> DEBUG --- WRITE_REPORT --- Wrote entire Transfer Report: {entire_report_path}")

#     # Then we may write individual files for each target_name under sequence_id
#     for target_name, target_data in transfer_dict[pfam_id]['sequence_id'].items():
#         if target_data:
#             json_data = json.dumps(target_data, indent=4)
#             logger.debug(f"---> DEBUG --- WRITE_REPORT --- Data to write for {target_name}!")
#         else:
#             json_data = json.dumps({"Annotations": "None"}, indent=4)
#             logger.debug(f"---> DEBUG --- WRITE_REPORT --- No data to write for {target_name}!")

#         safe_target_name = target_name.replace("|", "-")
#         target_report_filepath = os.path.join(output_dir, safe_target_name, pfam_id + "_report.json")
#         os.makedirs(os.path.join(output_dir, safe_target_name), exist_ok=True)

#         with open(target_report_filepath, 'w', encoding="utf-8") as report_file:
#             report_file.write(json_data)

#         logger.debug(f"---> DEBUG --- WRITE_REPORT --- Wrote Transfer Report for {target_name}-{pfam_id}: {target_report_filepath}")

#@measure_time_and_memory
##@profile
def map_and_filter_annot_pos(
    logger: logging.Logger,
    good_eco_codes: list,
    target_sequence: str,
    target_name: str,
    target_hit_start: int,
    target_hit_end: int,
    offset_start: int,
    offset_end: int,
    annot_sequence: str,
    entry_mnemo_name: str,
    entry_annotations: Dict[str, Any],
    transfer_dict: dict,
    processed_annotations: set = None,
    paired_annot_pos_str: str = None,
    caller_target_pos_str: str = None,
    paired_annotation_dict: Optional[Dict[str, Any]] = None) -> Optional[tuple[bool, dict]]:
    """
    Asks to map positions between annotated and target sequences.
    If target paired annot position and caller target position are provided,
    this is an associated call for a paired position sent to validate_paired_annotations.
    Otherwise, it is a call for a single position sent to validate_annotations.
    Makes counter_annot_pos and counter_target_pos to track the current position in both sequences.
    """
    counter_annot_pos = None
    counter_target_pos = None

    if paired_annot_pos_str and caller_target_pos_str and paired_annotation_dict:
        paired_tuple_res_hit_result_dict = validate_paired_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes,
            target_sequence=target_sequence,
            target_name=target_name,
            target_hit_start=target_hit_start,
            target_hit_end=target_hit_end,
            offset_start=offset_start,
            offset_end=offset_end,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name,
            paired_annotation_dict=paired_annotation_dict,
            entry_annotations=entry_annotations,
            counter_target_pos=counter_target_pos,
            counter_annot_pos=counter_annot_pos,
            paired_annot_pos_str=paired_annot_pos_str,
            caller_target_pos_str=caller_target_pos_str,
        )
        return paired_tuple_res_hit_result_dict
    try:
        validate_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes,
            target_sequence=target_sequence,
            target_name=target_name,
            target_hit_start=target_hit_start,
            target_hit_end=target_hit_end,
            offset_start=offset_start,
            offset_end=offset_end,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name,
            entry_annotations=entry_annotations,
            transfer_dict=transfer_dict,
            processed_annotations=processed_annotations,
            counter_target_pos=counter_target_pos,
            counter_annot_pos=counter_annot_pos
        )
        return None
    except Exception as e:
        logger.error(f"---> ERROR --- MAP_AND_FILTER --- Error in validate_annotations: {e}")
        traceback.print_exc()
        raise

def add_to_transfer_dict(
    hit: bool,
    logger: logging.Logger,
    transfer_dict: dict,
    target_name: str,
    target_sequence_continuous: str,
    target_hit_start: int,
    target_hit_end: int,
    anno_id: str,
    anno_total: dict,
    entry_mnemo_name: str,
    entry_primary_accession: str,
    paired_position_res_hit: Optional[bool] = None,
    paired_anno_id: Optional[str] = "",
    paired_anno_total: Optional[dict] = None) -> None:
    """
    Adds anno_total data to the transfer dictionary.
    The dictionary is structured by target name, hit | miss branches,
    target position, annotation ID (type + description).
    Further, you'll have essentials (type, description, count), evidence
    and various possible additional annotation details, mainly for BINDING annotations.
    Keeps tabs of how many times each annotation ID is found, as well as if paired data is provided
    (paired_position_res_hit, paired_anno_id, paired_anno_total), processes that data analogously after the early pair data.
    """
    try:
        additional_keys = {k: v for k, v in anno_total.items()
                           if k not in ['type', 'description', 'count', 'evidence',
                                         'target_position', 'annot_position', 'paired_target_position', 'annot_amino_acid', 'target_amino_acid', 'index_position']}
    except AttributeError as ae:
        logger.error("---> ERROR --- ADD_TO_TRANSFER_DICT --- AttributeError: Anno_total: %s", anno_total)
        logger.error(f"---> ERROR --- ADD_TO_TRANSFER_DICT --- AttributeError: {ae} - target_name: {target_name}, entry_mnemo_name: {entry_mnemo_name} \n", traceback.format_exc())
        raise

    paired_target_position_str = anno_total.get('paired_target_position', None)

    # Initialize transfer_dict structure if needed
    transfer_dict.setdefault("DOMAIN", {})
    transfer_dict["DOMAIN"].setdefault("sequence_id", {})

    if target_name not in transfer_dict["DOMAIN"]["sequence_id"]:
        transfer_dict["DOMAIN"]["sequence_id"][target_name] = {
            "sequence": target_sequence_continuous,
            "length": len(target_sequence_continuous),
            "hit_start": target_hit_start,
            "hit_end": target_hit_end,
            "annotations": {'positions': {}, 'indices': {'matches': set(), 'misses': set()}},
            "conservations": {'positions': {}, 'indices': {'matches': set(), 'misses': set()}},
            "position_conversion": {
                "target_to_aln": {},
                "aln_to_target": {}
            }
        }

    # Process main annotation
    _add_single_annotation(
        hit=hit,
        transfer_dict=transfer_dict,
        target_name=target_name,
        anno_id=anno_id,
        anno_total=anno_total,
        entry_mnemo_name=entry_mnemo_name,
        entry_primary_accession=entry_primary_accession,
        additional_keys=additional_keys
    )

    # Process paired annotation if provided
    if all(v is not None for v in [paired_position_res_hit, paired_anno_id, paired_anno_total, paired_target_position_str]):
        try:
            paired_additional_keys = {k: v for k, v in paired_anno_total.items()
                                    if k not in ['type', 'description', 'count', 'evidence',
                                                 'target_position', 'annot_position', 'paired_target_position', 'annot_amino_acid', 'target_amino_acid', 'index_position']}
        except AttributeError as ae:
            logger.error("---> ERROR --- ADD_TO_TRANSFER_DICT --- AttributeError: paired pair anno total: %s", paired_anno_total)
            logger.error(f"---> ERROR --- ADD_TO_TRANSFER_DICT --- AttributeError: {ae} - target_name: {target_name}, entry_mnemo_name: {entry_mnemo_name} \n", traceback.format_exc())
            raise

        _add_single_annotation(
            hit=paired_position_res_hit,
            transfer_dict=transfer_dict,
            target_name=target_name,
            anno_id=paired_anno_id,
            anno_total=paired_anno_total,
            entry_mnemo_name=entry_mnemo_name,
            entry_primary_accession=entry_primary_accession,
            additional_keys=paired_additional_keys
        )

        ### DELETE AFTER TESTS
        # print("Paired Target Position String:")
        # print(json.dumps(paired_target_position_str, indent=4))

        # print("Additional Keys:")
        # print(json.dumps(additional_keys, indent=4))

        # print("Transfer Dictionary for Target Name:")
        # print(json.dumps(transfer_dict[target_name], indent=4))


def _add_single_annotation(
    hit: bool,
    transfer_dict: dict,
    target_name: str,
    anno_id: str,
    anno_total: dict,
    entry_mnemo_name: str,
    entry_primary_accession: str,
    additional_keys: dict) -> None:
    """Helper function to add a single annotation to transfer_dict"""

    # Extract positions in target and annotated numbering - for early pair data
    target_position_str = anno_total.get('target_position', None)
    annot_position_str = anno_total.get('annot_position', None)
    index_position_str = anno_total.get('index_position', None)
    paired_target_position_str = anno_total.get('paired_target_position', None)

    # Information to use in basic populating
    evidence_value = anno_total.get('evidence', None)
    type_value = anno_total.get('type', None)
    description_value = anno_total.get('description', None)
    count_value = anno_total.get('count', None)
    annot_amino_value = anno_total.get('annot_amino_acid', None)
    target_amino_value = anno_total.get('target_amino_acid', None)

    if target_position_str not in transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions']:
        transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str] = {}
        # NEW, Add target position to one of the indices according to hit type
        index_type = 'matches' if hit else 'misses'
        transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['indices'][index_type].add(target_position_str)
        # NEW, Add target position and index to position_conversion
        transfer_dict['DOMAIN']['sequence_id'][target_name]['position_conversion']['target_to_aln'][target_position_str] = index_position_str
        if index_position_str not in transfer_dict['DOMAIN']['sequence_id'][target_name]['position_conversion']['aln_to_target']:
            transfer_dict['DOMAIN']['sequence_id'][target_name]['position_conversion']['aln_to_target'][index_position_str] = set()
        transfer_dict['DOMAIN']['sequence_id'][target_name]['position_conversion']['aln_to_target'][index_position_str].add(target_position_str)

    if anno_id not in transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str]:
        essentials = {
            'type': type_value,
            'description': description_value,
            'count': count_value,
            'annot_amino_acid': annot_amino_value,
            'target_amino_acid': target_amino_value
        }
        transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id] = {}
        transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id].setdefault('essentials', essentials)
    else:
        transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['essentials']['count'] += 1

    # Was the comparison between the two sequences a hit or a miss? Keep at {} level
    transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id].setdefault('hit', hit)

    if evidence_value:
        if 'evidence' not in transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]:
            transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id].setdefault('evidence', {})
        if evidence_value not in transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['evidence']:
            transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['evidence'][evidence_value] = {
                "rep_primary_accession": entry_primary_accession,
                "rep_mnemo_name": entry_mnemo_name,
                "count": 1
            }
        else:
            transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['evidence'][evidence_value]["count"] += 1

    if paired_target_position_str:
        if 'paired_position' not in transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]:
            transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id].setdefault('paired_position', {})
        if paired_target_position_str not in transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['paired_position']:
            transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['paired_position'][paired_target_position_str] = {
                "rep_primary_accession": entry_primary_accession,
                "rep_mnemo_name": entry_mnemo_name,
                "count": 1
            }
        else:
            transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['paired_position'][paired_target_position_str]["count"] += 1

    # New, helps analyze annotations that, in their description,
    # relate positions in annotated sequence uniprot numbering.
    # Note that those with sequence ranges are a bit tricky yet, but I'll leave it as is for now
    if type_value in ['CROSSLNK', 'DISULFID', 'MUTAGEN'] and description_value and re.search(r'\b[A-Z]-\d+', description_value):
        additional_keys['annot_position'] = annot_position_str

    if additional_keys and type_value in ['BINDING', 'ACT_SITE', 'CROSSLNK', 'DISULFID', 'MUTAGEN']:
        if 'additional_keys' not in transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]:
            transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id].setdefault('additional_keys', {})
        for key, value in additional_keys.items():
            if key not in transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['additional_keys']:
                transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['additional_keys'][key] = {}
            if value not in transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['additional_keys'][key]:
                transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['additional_keys'][key][value] = {
                    "rep_primary_accession": entry_primary_accession,
                    "rep_mnemo_name": entry_mnemo_name,
                    "count": 1
                }
            else:
                transfer_dict['DOMAIN']['sequence_id'][target_name]['annotations']['positions'][target_position_str][anno_id]['additional_keys'][key][value]["count"] += 1


#@measure_time_and_memory
#@profile
def make_anno_total_dict(
    good_eco_codes: list,
    entry_mnemo_name: str,
    annotation_dict: Dict[str, Any],
    counter_target_pos_str: str,
    counter_annot_pos_str: str,
    index: int,
    target_amino: str,
    logger: Optional[logging.Logger] = None, # Delete in production
    entry_annotations: Optional[Dict[str, Any]] = None,
    caller_target_pos_str: Optional[str] = None) -> Dict[str, Any]:
    """
    For any given annotation in appropriate format (list with 1 dict(annotation_dict) inside), paired or not, produces an anno_total dict and
    associated data to return as a dict to be extracted from in accordance to the calling function's needs.
    Used by process_annotation and validate_paired_annotations.
    """
    anno_type = ""
    anno_id = ""
    anno_total = None
    annotation = annotation_dict

    # DEBUGGING, still keeping logger and entry_annotations for now, will DELETE in PRODUCTION!
    logger.debug(f"---> DEBUG --- MAKE_ANNO --- Annotation Dict: {annotation_dict}")
    logger.debug(f"---> DEBUG --- MAKE_ANNO --- Entry Annotations List: {entry_annotations.get(counter_annot_pos_str)}")
    logger.debug(f"---> DEBUG --- MAKE_ANNO --- Good ECO codes : {good_eco_codes}")

    paired_annot_pos_str = annotation.get('paired_position', None)
    anno_type = annotation['type']
    anno_desc = annotation.get('description', None)
    anno_id = f"{anno_type} | {anno_desc}"
    anno_count = 1
    anno_evidence = {entry_mnemo_name: annotation.get('evidence', None)}
    annot_amino_acid = annotation.get('aminoacid', None)

    if anno_evidence[entry_mnemo_name] and good_eco_codes:
        if isinstance(anno_evidence[entry_mnemo_name], list):
            evidence_string = anno_evidence[entry_mnemo_name][0]
        else:
            evidence_string = anno_evidence[entry_mnemo_name]
        evidence_items = evidence_string.split(', ')
        valid_evidence_items = []
        for item in evidence_items:
            code, _, _ = item.strip().partition('|')  # Split by the first occurrence of '|'
            if code.strip() in good_eco_codes:
                valid_evidence_items.append(item.strip())  # Preserve the whole item if the code is valid
        anno_evidence[entry_mnemo_name] = ', '.join(valid_evidence_items)

    if anno_evidence[entry_mnemo_name]:
        anno_total = {
            "type": anno_type,
            "description": anno_desc,
            "count": anno_count,
            "evidence": anno_evidence[entry_mnemo_name],
            "target_position": counter_target_pos_str,
            "annot_position": counter_annot_pos_str,
            "index_position": str(index),
            "annot_amino_acid": annot_amino_acid,
            "target_amino_acid": target_amino
        }
        if paired_annot_pos_str and caller_target_pos_str:
            anno_total['paired_target_position'] = caller_target_pos_str
        # Changed to include ACT_SITE, because of ligand_ccd_id in BioLiP sourced annotations
        if anno_type in ['BINDING', 'ACT_SITE']:
            keys_to_exclude = {'type', 'description', 'count', 'evidence',
                               'paired_position', 'aminoacid', 'entry'}
            for key, value in annotation.items():
                if key not in keys_to_exclude:
                    anno_total[key] = value

    return {
        "annotation": annotation,
        "anno_type": anno_type,
        "anno_id": anno_id,
        "anno_total": anno_total,
        "paired_annot_pos_str": paired_annot_pos_str
        }

#@measure_time_and_memory
#@profile
def process_annotation(
    res_hit: bool,
    logger: logging.Logger,
    good_eco_codes: list,
    entry_mnemo_name: str,
    target_name: str,
    target_hit_start: int,
    target_hit_end: int,
    annotation_dict: Dict[str, Any],
    entry_annotations: Dict[str, Any],
    transfer_dict: dict,
    target_sequence: str,
    offset_start: int,
    offset_end: int,
    annot_sequence: str,
    counter_target_pos_str: str,
    counter_annot_pos_str: str,
    index: int,
    target_amino: str,
    processed_annotations: set,
    annotation_key: tuple
) -> None:
    """
    Processes annotations to add them to the transfer dictionary by calling make_anno_total dict
    with anno_total containing at least type, description, count, and evidence, plus associated data.
    Additionally, it adds paired position if the other member also hits (successfull map_and_filter_annot_pos for it),
    and adds additional keys for BINDING annotations, if present.
    """
    result_dict = make_anno_total_dict(
        good_eco_codes=good_eco_codes,
        entry_mnemo_name=entry_mnemo_name,
        annotation_dict=annotation_dict,
        counter_target_pos_str=counter_target_pos_str,
        counter_annot_pos_str=counter_annot_pos_str,
        index=index,
        target_amino=target_amino,
        logger=logger, # Delete in production
        entry_annotations=entry_annotations
    )
    annotation = result_dict['annotation']
    anno_type = result_dict['anno_type']
    anno_id = result_dict['anno_id']
    anno_total = result_dict['anno_total']
    paired_annot_pos_str = result_dict['paired_annot_pos_str']

    if anno_total:
        entry_primary_accession = annotation.get('entry', None)
        target_sequence_continuous = ''.join([char.upper() for char in target_sequence if char.isalpha()])
        if anno_type in ['DISULFID', 'CROSSLNK', 'SITE', 'BINDING'] and paired_annot_pos_str is not None:
            paired_annotations = entry_annotations.get(paired_annot_pos_str, [])
            paired_annotation_dict = next(
                (annotation for annotation in paired_annotations if annotation['type'] == anno_type), None
            )
            if paired_annotation_dict:
                paired_position_res_hit, paired_result_dict = map_and_filter_annot_pos(
                    logger=logger,
                    good_eco_codes=good_eco_codes,
                    target_sequence=target_sequence,
                    target_name=target_name,
                    target_hit_start=target_hit_start,
                    target_hit_end=target_hit_end,
                    offset_start=offset_start,
                    offset_end=offset_end,
                    annot_sequence=annot_sequence,
                    entry_mnemo_name=entry_mnemo_name,
                    entry_annotations=entry_annotations,
                    transfer_dict=transfer_dict,
                    paired_annot_pos_str=paired_annot_pos_str,
                    caller_target_pos_str=counter_target_pos_str,
                    paired_annotation_dict=paired_annotation_dict
                )
                paired_anno_total = paired_result_dict['anno_total']
                paired_anno_id = paired_result_dict['anno_id']
                paired_target_position_str = paired_anno_total['target_position']
            else:
                paired_position_res_hit = False
                paired_anno_total = None
                paired_anno_id = None
                paired_target_position_str = None
            anno_total['paired_target_position'] = paired_target_position_str
            logger.debug(f"---> DEBUG --- PROCESS_ANNOT --- Paired Position Valid for target {target_name} and annotated {entry_mnemo_name}")
            try:
                add_to_transfer_dict(
                    hit=res_hit,
                    logger=logger,
                    transfer_dict=transfer_dict,
                    target_name=target_name,
                    target_sequence_continuous=target_sequence_continuous,
                    target_hit_start=target_hit_start,
                    target_hit_end=target_hit_end,
                    anno_id=anno_id,
                    anno_total=anno_total,
                    entry_mnemo_name=entry_mnemo_name,
                    entry_primary_accession=entry_primary_accession,
                    paired_position_res_hit=paired_position_res_hit,
                    paired_anno_id=paired_anno_id,
                    paired_anno_total=paired_anno_total,
                )
                processed_annotations.add(annotation_key)
                if paired_target_position_str:
                    paired_annotation_key = (
                        entry_mnemo_name,
                        target_name,
                        paired_annot_pos_str,
                        paired_target_position_str,
                        anno_type
                    )
                    processed_annotations.add(paired_annotation_key)
                logger.debug(
                    f"---> DEBUG --- PROCESS_ANNOT --- SUCCESS PA Added to transfer_dict for target {target_name} and "
                    f"annotated {entry_mnemo_name} at target {counter_target_pos_str} and "
                    f"annotated {counter_annot_pos_str} --- ORIGIN: Type {anno_type} and "
                    f"Evidence {anno_total['evidence']}"
                )
            except Exception as e:
                logger.error(f"---> ERROR --- PROCESS_ANNOT --- Error in add_to_transfer_dict: {e}")
                traceback.print_exc()
                raise
        else:
            add_to_transfer_dict(
                hit=res_hit,
                logger=logger,
                transfer_dict=transfer_dict,
                target_name=target_name,
                target_sequence_continuous=target_sequence_continuous,
                target_hit_start=target_hit_start,
                target_hit_end=target_hit_end,
                anno_id=anno_id,
                anno_total=anno_total,
                entry_mnemo_name=entry_mnemo_name,
                entry_primary_accession=entry_primary_accession
            )
            processed_annotations.add(annotation_key)
            logger.debug(
                f"---> DEBUG --- PROCESS_ANNOT --- SUCCESS PA Added to transfer_dict for target {target_name} and "
                f"annotated {entry_mnemo_name} at target {counter_target_pos_str} and "
                f"annotated {counter_annot_pos_str} --- ORIGIN: Type {anno_type} and "
                f"Evidence {anno_total['evidence']}"
            )
    else:
        logger.debug(f"---> DEBUG --- PROCESS_ANNOT --- NO ANNOTATION FOR SINGLE --- Anno Total was None for target {target_name} and annotated {entry_mnemo_name} at target {counter_target_pos_str} and annotated {counter_annot_pos_str} --- ORIGIN: No evidence of desired type")
        # PENDING - Conservation: Add code to deal with conserved positions with no annotations
        processed_annotations.add(annotation_key)

#@measure_time_and_memory
#@profile
def validate_paired_annotations(
    logger: logging.Logger,
    good_eco_codes: list,
    target_sequence: str,
    target_name: str,
    target_hit_start: int,
    target_hit_end: int,
    offset_start: int,
    offset_end: int,
    annot_sequence: str,
    entry_mnemo_name: str,
    paired_annotation_dict: Dict[str, Any],
    entry_annotations: Dict[str, Any],
    counter_target_pos: Optional[None],
    counter_annot_pos: Optional[None],
    paired_annot_pos_str: str,
    caller_target_pos_str: str,
) -> tuple[bool, dict | None]:
    """
    Validates the 2nd member of a paired position by checking
    if the target paired annot position from caller is a valid match.
    Success and failure both mean adding the annotation to the transfer dictionary
    in the process_annotation parent function, difference lies in where (match | miss).
    """
    paired_position_res_hit = False
    paired_result_dict = {}
    paired_annot_pos_int = int(paired_annot_pos_str)
    last_target_pos = False

    logger.debug(f"---> DEBUG --- VAL_PAIRED --- Running for target {target_name}/{target_hit_start}-{target_hit_end} and annotated {entry_mnemo_name}/{offset_start}-{offset_end}")

    while not last_target_pos:
        for index, counter_annot_pos, counter_target_pos, char_annot, char_target in iterate_aligned_sequences(
            seq_a=annot_sequence,
            seq_b=target_sequence,
            start_a=offset_start,
            start_b=target_hit_start,
            end_a=offset_end,
            end_b=target_hit_end
            ):
                if counter_annot_pos is None or counter_target_pos is None:
                    continue

                if counter_annot_pos == paired_annot_pos_int:
                    counter_target_pos_str = str(counter_target_pos)

                    # DEBUGGING INFO
                    logger.debug("\n --- DEBUG --- VAL_PAIRED --- Helpful Info for paired match | miss \n")
                    start_index = max(0, index - 3)
                    end_index = min(len(annot_sequence), index + 4)
                    annot_window = annot_sequence[start_index:end_index]
                    target_window = target_sequence[start_index:end_index]
                    logger.debug(f"---> DEBUG --- VAL_PAIRED --- Counter Tar Pos {counter_target_pos_str} and amino acid {target_sequence[index]} + Counter annot Pos (annotated) {str(counter_annot_pos)} and amino acid {annot_sequence[index]}")
                    logger.debug(f"---> DEBUG --- VAL_PAIRED --- Annot Window: {annot_window} + Target Window: {target_window}")

                    paired_target_amino = target_sequence[index]
                    paired_position_res_hit = bool(target_sequence[index] == annot_sequence[index])

                    paired_result_dict = make_anno_total_dict(
                        good_eco_codes=good_eco_codes,
                        entry_mnemo_name=entry_mnemo_name,
                        annotation_dict=paired_annotation_dict,
                        counter_target_pos_str=counter_target_pos_str, # kinda paired_target_pos_str
                        counter_annot_pos_str=paired_annot_pos_str,
                        index=index,
                        target_amino=paired_target_amino,
                        logger=logger, # Delete in production
                        entry_annotations=entry_annotations,
                        caller_target_pos_str=caller_target_pos_str,
                    )

                    return paired_position_res_hit, paired_result_dict

                if counter_target_pos == target_hit_end or counter_annot_pos == offset_end:
                    last_target_pos = True
                    break

    return paired_position_res_hit, paired_result_dict

#@measure_time_and_memory
#@profile
def validate_annotations(
    logger: logging.Logger,
    good_eco_codes: list,
    target_sequence: str,
    target_name: str,
    target_hit_start: int,
    target_hit_end: int,
    offset_start: int,
    offset_end: int,
    annot_sequence: str,
    entry_mnemo_name: str,
    entry_annotations: Dict[str, Any],
    transfer_dict: dict,
    processed_annotations: set,
    counter_target_pos: None,
    counter_annot_pos: None) -> None:
    """
    Iterate on both annotated and target sequences to find columns where the annotated sequence has annotations. See if the amino acid at target sequence is the same as the one on the annotated sequence, storing the result in a boolean "res_hit".
    Regardless of hit, send to process_annotation, granted that processed_annotations does not contain the annotation key (entry_mnemo_name, target_name, counter_annot_pos_str, counter_target_pos_str and anno_type).
    """

    annotated_annot_pos_list = [int(key) for key in entry_annotations.keys() if key != '0']
    if not any(pos in range(offset_start, offset_end + 1) for pos in annotated_annot_pos_list):
        return
    last_target_pos = False

    while not last_target_pos:
        for index, counter_annot_pos, counter_target_pos, char_annot, char_target in iterate_aligned_sequences(
            seq_a=annot_sequence,
            seq_b=target_sequence,
            start_a=offset_start,
            start_b=target_hit_start,
            end_a=offset_end,
            end_b=target_hit_end):

            if counter_annot_pos is None or counter_target_pos is None:
                continue
            counter_annot_pos_str = str(counter_annot_pos)

            if counter_annot_pos_str in entry_annotations:
                counter_target_pos_str = str(counter_target_pos)

                # DEBUGGING INFO
                logger.debug(f"\n --- DEBUG --- VAL_ANNOTS --- Helpful Info for Single/Caller Annotation Match \n NAMES: target {target_name} and annot {entry_mnemo_name} \n POSITIONS: target {counter_target_pos_str} and annot {counter_annot_pos_str} \n")
                start_index = max(0, index - 3)
                end_index = min(len(annot_sequence), index + 4)
                annot_window = annot_sequence[start_index:end_index]
                target_window = target_sequence[start_index:end_index]
                logger.debug(f"---> DEBUG --- VAL_ANNOTS --- Counter annot Pos {counter_annot_pos_str} and amino acid {annot_sequence[index]} + Counter Tar Pos {str(counter_target_pos)} and amino acid {target_sequence[index]}")
                logger.debug(f"---> DEBUG --- VAL_ANNOTS --- Annot Window: {annot_window} + Target Window: {target_window}")

                target_amino = target_sequence[index]
                res_hit = bool(target_sequence[index] == annot_sequence[index])
                for annotation_dict in entry_annotations[counter_annot_pos_str]:
                    anno_type = annotation_dict.get('type', None)
                    annotation_key = (
                        entry_mnemo_name,
                        target_name,
                        counter_annot_pos_str,
                        counter_target_pos_str,
                        anno_type
                    )
                    # print(f"Annotation Key: {annotation_key}")
                    # print(f"Index: {index}")
                    if annotation_key in processed_annotations:
                        continue
                    try:
                        process_annotation(
                            res_hit=res_hit,
                            logger=logger,
                            good_eco_codes=good_eco_codes,
                            entry_mnemo_name=entry_mnemo_name,
                            target_name=target_name,
                            target_hit_start=target_hit_start,
                            target_hit_end=target_hit_end,
                            annotation_dict=annotation_dict,
                            entry_annotations=entry_annotations,
                            transfer_dict=transfer_dict,
                            target_sequence=target_sequence,
                            offset_start=offset_start,
                            offset_end=offset_end,
                            annot_sequence=annot_sequence,
                            counter_target_pos_str=counter_target_pos_str,
                            counter_annot_pos_str=counter_annot_pos_str,
                            index=index,
                            target_amino=target_amino,
                            processed_annotations=processed_annotations,
                            annotation_key=annotation_key
                        )
                    except Exception as e:
                        logger.error(f"---> ERROR --- VAL_ANNOTS --- Error in validate_annotations: {e}")
                        traceback.print_exc()
                        raise

            if counter_target_pos == target_hit_end or counter_annot_pos == offset_end:
                last_target_pos = True
                break


def main(logger: logging.Logger):
    """Main function, initializes this script"""
    args = parse_arguments()
    dom_align = args.dom_align
    resource_dir = args.resource_dir
    output_dir = args.output_dir
    good_eco_codes = args.eco_codes
    logger.info("---> MAIN --- Running transfer_annotations.py for %s --- ", dom_align)

    pfam_id = get_pfam_id_from_hmmalign_result(dom_align)
    annotations_filepath, conservations_filepath = get_annotation_filepath(resource_dir, pfam_id)
    hmmalign_lines, annotations = read_files(dom_align, annotations_filepath)

    try:
        logger.debug(f"---> DEBUG --- MAIN --- Good ECO Codes to Filter by {good_eco_codes}")
        transfer_dict = find_and_map_annots(logger, hmmalign_lines, annotations, good_eco_codes)
        # DEBUG, Delete Later
        # debug_transfer_dict = convert_sets_to_lists(transfer_dict)
        logger.info(f"Annotations filepath: {annotations_filepath}")
        # logger.info(f"Transfer dictionary: {json.dumps(transfer_dict, indent=4)}")
        logger.info(f"Writing report to: {os.path.join(output_dir, 'MCRB_ECOLI', 'PF00244_report.json')}")
        if not transfer_dict:
            logger.info("---> MAIN --- Transfer Dict was EMPTY")
        else:
            logger.info("---> MAIN --- Transfer Dict FILLED")
    except (KeyError, IndexError, AttributeError) as e:
        error_info = traceback.format_exc()
        logger.error("---> MAIN --- ERROR transferring annotations for Pfam ID %s: %s\n%s", pfam_id, e, error_info)
        raise
    improved_transfer_dict = cleanup_improve_transfer_dict(transfer_dict, pfam_id, hmmalign_lines, conservations_filepath, annotations_filepath, output_dir)
    write_reports(logger, improved_transfer_dict, output_dir)
    # write_reports(logger, transfer_dict, output_dir)

if __name__ == '__main__':
    outer_args = parse_arguments()
    outer_logger = configure_logging(outer_args.log)
    main(outer_logger)
