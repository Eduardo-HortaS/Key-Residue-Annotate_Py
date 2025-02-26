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
parse_arguments -> get_logger -> get_multi_logger ->
get_pfam_id_from_hmmalign_result -> get_annotation_filepath -> read_files ->

(Annotations file missing?)
YES -> annotations = {"sequence_id": {}} -> setup_for_conservations_only -> extract_target_info_from_hmmalign ->
cleanup_improve_transfer_dict -> read_conservations_and_annotations -> get_alignment_sequences ->
populate_conservation -> gather_go_terms_for_target -> populate_go_data_for_annotations -> write_reports

NO -> find_and_map_annots -> extract_target_info_from_hmmalign -> map_and_filter_annot_pos ->
validate_annotations -> iterate_aligned_sequences -> process_annotation -> make_anno_total_dict

(Anno_Total is None?)
YES -> DO NOTHING, JUST TRACK --- position effectively skipped

NO -> PROCEED CHECK PAIREABLE ANNOTATIONS

(paireable annotation?)
NO -> add_to_transfer_dict -> _add_single_annotation ->
update_position_ranges -> merge_adjacent_ranges

YES -> map_and_filter_annot_pos -> validate_paired_annotations ->
iterate_aligned_sequences -> make_anno_total_dict -> process_annotation -> add_to_transfer_dict ->
_add_single_annotation -> update_position_ranges -> merge_adjacent_ranges

ANNOTATIONS PRESENT FINAL STEPS -> cleanup_improve_transfer_dict -> read_conservations_and_annotations ->
get_alignment_sequences -> populate_conservation -> gather_go_terms_for_target
-> populate_go_data_for_annotations -> write_reports

"""

import os
import json
import re
import argparse
import logging
import traceback
import copy
from typing import Any, Dict, Optional, Generator, Tuple, Callable
import pandas as pd
from utils import get_logger, get_multi_logger
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
    "Generates a temporary multifasta for running hmmalign using a hits per domain JSON.")
    parser.add_argument("-iA", "--dom-align", required=True, type=str, help="Path to domain's hmmalign alignment")
    parser.add_argument("-r", "--resource-dir", required=True, type=str, help="Resource dir path")
    parser.add_argument("-d", "--domain-accession", help="Domain accession for scoped logging", required=True, type=str)
    parser.add_argument("-o", "--output-dir", required=True, type=str, help="Output dir path")
    parser.add_argument("-e", "--eco-codes", required=False, default=[], nargs="*", help="Space-separated ECO codes to filter annotations")
    parser.add_argument("-l", "--log", required=False, default="logs/transfer_annotations.log", type=str, help="Log path")

    args = parser.parse_args()

    # Clean eco codes - remove any quotes, brackets and commas
    if args.eco_codes:
        args.eco_codes = [code.strip('",[]') for code in args.eco_codes]

    return args

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

    try:
        with open(annotations_filepath, 'r', encoding="utf-8") as annotations_file:
            annotations = json.load(annotations_file)
    except (FileNotFoundError, IOError):
        annotations = {"sequence_id": {}}

    return hmmalign_lines, annotations

### New Helper Functions
def iterate_aligned_sequences(
    source_sequence: str,
    target_sequence: str,
    source_start: int,
    target_start: int,
    source_end: int,
    target_end: int
) -> Generator[Tuple[int, Optional[int], Optional[int], str, str], None, None]:
    """
    Yields index, current positions (source_pos, target_pos), and characters (source_char, target_char).
    Positions are incremented only when the respective character is alphabetic.
    For target sequence, positions are None when encountering deletion relative to consensus ('-') or insertion gap markers ('.').
    """
    source_pos = None
    target_pos = None
    last_valid_target_pos = None
    last_valid_source_pos = None

    for index, (source_char, target_char) in enumerate(zip(source_sequence, target_sequence)):
        # Increment counters for any letter
        # Handle conservation or annotated sequence
        if source_char.isalpha():
            source_pos = source_start if last_valid_source_pos is None else last_valid_source_pos + 1
            last_valid_source_pos = source_pos
        else: # Gap character in source, set to None
            source_pos = None

        if target_char.isalpha():
            target_pos = target_start if last_valid_target_pos is None else last_valid_target_pos + 1
            last_valid_target_pos = target_pos
        else: # Deletion ("-") or Insertion (".") gap - set None
            target_pos = None

        yield index, source_pos, target_pos, source_char, target_char

        if (source_pos is not None and source_pos == source_end) or (last_valid_target_pos is not None and last_valid_target_pos == target_end):
            break

def extract_target_info_from_hmmalign(logger: logging.Logger, multi_logger: Callable, hmmalign_lines: list) -> dict:
    """Extracts target sequence information from hmmalign lines for both annotations and conservations.
    Args:
        hmmalign_lines (list): List of lines from the hmmalign alignment file.
        logger (logging.Logger): Logger object for logging.
        multi_logger (Callable): Callable for logging to multiple loggers.
    Returns:
        dict: Format {target_name: {target_hit_interval: [(target_hit_start, target_hit_end, target_hit_sequence),...]}}
        where target_hit_interval is "start-end"
    """
    target_info = {}
    for line in hmmalign_lines:
        if "target/" in line and not line.startswith("#"):
            logger.debug("TRANSFER_ANNOTS --- EXTRACT_TARGET_INFO --- Target Line: \n %s", line)
            try:
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
                    target_info[target_name][target_hit_interval].append(
                        (target_hit_start, target_hit_end, target_hit_sequence)
                    )
            except(IndexError, ValueError) as e:
                multi_logger("warning", "TRANSFER_ANNOTS --- EXTRACT_TARGET_INFO --- Invalid format in line: %s. Error: %s", line, e)
                continue
    return target_info

def setup_for_conservations_only(logger: logging.Logger, multi_logger: Callable, hmmalign_lines: list, pfam_id: str) -> dict:
    """Creates transfer dictionary structure for cases where only conservation data is available.
    Args:
        logger (logging.Logger): Logger object for logging.
        multi_logger (Callable): Callable for logging to multiple loggers.
        hmmalign_lines (list): List of lines from the hmmalign alignment file.
        pfam_id (str): Pfam domain accession being processed.
    Returns:
        dict: Transfer dict with necessary keys for conservations: sequence_id, hit_intervals,
        sequence, length, hit_start/end, annotations (stays empty later),
        conservations, position_conversion, annotation_ranges.
    """

    transfer_dict = { pfam_id: {"sequence_id": {}}}
    target_info = extract_target_info_from_hmmalign(logger, multi_logger, hmmalign_lines)

    if not target_info:
        multi_logger("error", "TRANSFER_ANNOTS --- SETUP_CONS_ONLY --- No target sequences found in hmmalign lines!")
        transfer_dict = {}
        return transfer_dict

    for target_name, target_hit_interval in target_info.items():
        if target_name not in transfer_dict[pfam_id]["sequence_id"]:
            transfer_dict[pfam_id]["sequence_id"][target_name] = {"hit_intervals": {}}

        for interval_key, target_info_list in target_hit_interval.items():
            for target_hit_start, target_hit_end, target_sequence in target_info_list:
                target_sequence_continuous = ''.join([char.upper() for char in target_sequence if char.isalpha()])
                interval_dict = {
                    "sequence": target_sequence_continuous,
                    "length": len(target_sequence_continuous),
                    "hit_start": target_hit_start,
                    "hit_end": target_hit_end,
                    "annotations": {},
                    "conservations": {
                        "positions": {},
                        "indices": {"matches": set(), "misses": set()}
                    },
                    "position_conversion": {
                        "target_to_aln": {},
                        "aln_to_target": {}
                    },
                    "annotation_ranges": {},
                    "conservation_ranges": {}
                }
                transfer_dict[pfam_id]["sequence_id"][target_name]["hit_intervals"][interval_key] = interval_dict

    return transfer_dict

#@measure_time_and_memory
#@profile
def find_and_map_annots(
    logger: logging.Logger,
    multi_logger: Callable,
    hmmalign_lines: list,
    annotations: dict,
    good_eco_codes: list) -> dict:
    """
    -> extract_target_info_from_hmmalign() = From the hmmalign alignment lines,
    finds target lines and stores data for them in target_info,
    organized by target name and hit interval. Then, finds seed sequences in hmmalign alignment
    for which we have annotations and calls map_and_filter_annot_pos() to map positions
    by iterating on both target and seed sequences. Makes processed_annotations set
    to keep track of visited positions across calls, where relevant, and sends it to map_and_filter_annot_pos.

    Args:
        logger (logging.Logger): Logger object for logging.
        multi_logger (Callable): Callable for logging to multiple loggers.
        hmmalign_lines (list): List of lines from the hmmalign alignment file.
        annotations (dict): Loaded annotations JSON object.
        good_eco_codes (list): List of ECO codes to filter annotations.

    Returns:
        dict: Transfer dictionary containing annotations for each sequence-domain pair.
    """
    transfer_dict = {}
    processed_annotations = set()

    target_info = extract_target_info_from_hmmalign(logger, multi_logger, hmmalign_lines)

    if not target_info:
        multi_logger("error", "TRANSFER_ANNOTS --- FIND_AND_MAP --- No target sequences found in hmmalign lines!")
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
                            map_and_filter_annot_pos(
                                logger=logger,
                                multi_logger=multi_logger,
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
                            logger.debug(f"TRANSFER_ANNOTS --- FIND_AND_MAP --- Mapped, filtered and possibly added to transfer_dict: target {target_name} and annotated {entry_mnemo_name} at target hit interval {target_hit_interval}")
    return transfer_dict

def read_conservations_and_annotations(conservations_filepath: str, annotations_filepath: str) -> tuple[dict, dict]:
    """Reads conservations and annotations JSON files.

    Args:
        conservations_filepath: Path to conservations JSON
        annotations_filepath: Path to annotations JSON

    Returns:
        tuple: (conservations, annotations) where:
            conservations: {"sequence_id/start-end": {"position": score}}
            annotations: {"sequence_id": {"position": [{"type": str}]}}

    Note:
        Returns empty structures for missing or invalid files
    """
    try:
        with open(conservations_filepath, 'r', encoding="utf-8") as conservations_file:
            conservations = json.load(conservations_file)
            if not conservations:
                conservations = {"sequence_id/range": {}}
    except (FileNotFoundError, IOError):
        conservations = {"sequence_id/range": {}}

    try:
        with open(annotations_filepath, 'r', encoding="utf-8") as annotations_file:
            annotations = json.load(annotations_file)
            if not annotations:
                annotations = {"sequence_id": {}}
    except (FileNotFoundError, IOError):
        annotations = {"sequence_id": {}}

    return conservations, annotations

def parse_go_annotations(go_column: str) -> list:
    """Extracts GO terms from InterProScan TSV column.

    Args:
        go_column: GO terms string (e.g., "GO:0005515|GO:0006302")

    Returns:
        list: Clean GO terms (e.g., ["GO:0005515", "GO:0006302"])
    """
    if not go_column or go_column == "-" or not go_column.strip():
        return []
    return [term.split("(")[0].strip() for term in go_column.split("|") if term]

# NOTE: May be adjusted to improve matching accuracy later, based on empirical testing.
def check_interval_overlap(
    multi_logger: Callable, iprscan_start: int, iprscan_end: int,
    hit_start: int, hit_end: int, margin_percent: float = 0.1) -> bool:
    """
    Checks if iprscan interval substantially overlaps with hit interval.
    Uses percentage-based margins for flexibility in interval matching,
    considering the different algorithms used by iprscan to find domains.

    Args:
        multi_logger: Callable for logging to multiple loggers
        iprscan_start: Start position from InterProScan
        iprscan_end: End position from InterProScan
        hit_start: Start position from HMMER hit
        hit_end: End position from HMMER hit
        margin_percent: Allowed margin as fraction of hit length (default 10%)

    Returns:
        bool: True if intervals overlap, False otherwise
    """
    if not all(isinstance(x, (int, float)) and x >= 0
              for x in [iprscan_start, iprscan_end, hit_start, hit_end, margin_percent]):
        multi_logger("error", "TRANSFER_ANNOTS --- CHECK_INTERVAL --- All positions must be non-negative numbers")
        return False

    if not isinstance(margin_percent, float) or margin_percent <= 0 or margin_percent > 1:
        multi_logger("error", "TRANSFER_ANNOTS --- CHECK_INTERVAL --- margin_percent must be a float between 0 and 1 (exclusive, inclusive)")
        return False

    hit_length = hit_end - hit_start + 1
    margin = int(hit_length * margin_percent)

    # Check if iprscan interval substantially overlaps with hit interval
    return (iprscan_start >= hit_start - margin and
            iprscan_end <= hit_end + margin)

def gather_go_terms_for_target(
    multi_logger: Callable, target_name: str, pfam_id: str,
    go_terms_dir: str, interpro_conv_id: str, hit_start: int, hit_end: int) -> set:
    """
    Gathers GO terms for a single target sequence from iprscan.tsv if
    adequate (meet at least one of 3 conditions) lines are present.
    Note: GO terms are kept as a string with them separated by "|".
    The GO terms base dir is expected to be the same as the output dir,
    since this function requires the previous scripts to have been run.\
    Checks for:
    - Matching InterProScan lines with the PFAM ID or InterPro ID.
    - Matching intervals with the hit start and end.

    Args:
        logger (logging.Logger): Logger object for logging.
        target_name (str): Target sequence name.
        pfam_id (str): PFAM domain accession being processed.
        go_terms_dir (str): Directory containing target sequence subdirs where iprscan.tsv lies.
        interpro_conv_id (str): InterPro ID for the domain.
        hit_start (int): Start position of the hit.
        hit_end (int): End position of the hit.

    Returns:
        set: GO terms found for the target sequence.
    """
    sanitized_target_name = target_name.replace("|", "-")
    go_term_filepath = os.path.join(go_terms_dir, sanitized_target_name, "iprscan.tsv")
    target_gos = set()

    if not os.path.exists(go_term_filepath):
        multi_logger(
            "warning",
            "TRANSFER_ANNOTS --- GO_TERMS_TARGET --- Missing iprscan.tsv file for Pfam ID %s and Target Name %s at %s",
            pfam_id, sanitized_target_name, go_term_filepath
        )
        return target_gos

    found_matching_accession = False
    found_matching_interval = False

    with open(go_term_filepath, 'r', encoding="utf-8") as go_term_file:
        iprscan_df = pd.read_csv(go_term_file, sep='\t', header=None)
        iprscan_df.columns = [
            "Protein Accession", "Sequence MD5 digest", "Sequence Length", "Analysis",
            "Signature Accession", "Signature Description", "Start Location", "Stop Location",
            "Score", "Status", "Date", "InterPro Accession", "InterPro Description",
            "GO Annotations", "Pathways Annotations"
        ]

        for _, row in iprscan_df.iterrows():
            if interpro_conv_id:
                matches_accession = (row["Signature Accession"] == pfam_id or
                                    row["InterPro Accession"] == interpro_conv_id)
            else:
                matches_accession = row["Signature Accession"] == pfam_id

            matches_interval = check_interval_overlap(
                multi_logger, row["Start Location"], row["Stop Location"], hit_start, hit_end
            )

            # Track if we found any matches
            found_matching_accession = found_matching_accession or matches_accession
            found_matching_interval = found_matching_interval or matches_interval

            if matches_accession or matches_interval:
                if not pd.isna(row["GO Annotations"]) and not row["GO Annotations"].strip() == "-":
                    target_gos.update(parse_go_annotations(row["GO Annotations"]))

    if not (found_matching_accession or found_matching_interval):
        multi_logger(
            "warning",
            "TRANSFER_ANNOTS --- GO_TERMS_TARGET --- No usable InterProScan lines in iprscan.tsv for %s-%s",
            pfam_id, sanitized_target_name
        )
        return set()

    return target_gos


def get_alignment_sequences(
    hmmalign_lines: list, target_id: str, conservation_id: str,
    logger: logging.Logger, multi_logger: Callable
    ) -> tuple[str | None, str | None]:
    """From hmmalign_lines, extracts the target_seq and conservation_seq matching the given IDs.

    Args:
        hmmalign_lines: Alignment file content
        target_id: Target sequence identifier
        conservation_id: Conservation sequence identifier
        logger: Domain-specific logger
        multi_logger: Callable for logging to multiple loggers

    Returns:
        tuple: (target_sequence, conservation_sequence) or (None, None) on failure
    """
    target_seq = None
    conservation_seq = None

    try:
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
            multi_logger("error", "TRANSFER_ANNOTS --- GET_ALIGN_SEQS --- Target ID '%s' not found in hmmalign_lines.", target_id)
            return None, None
        if conservation_seq is None:
            multi_logger("error", "TRANSFER_ANNOTS --- GET_ALIGN_SEQS --- Conservation ID '%s' not found in hmmalign_lines.", conservation_id)
            return None, None

        logger.debug("TRANSFER_ANNOTS --- GET_ALIGN_SEQS --- Successfully extracted alignment sequences for target '%s'", target_id)
        return target_seq, conservation_seq

    except Exception as e:
        multi_logger("error", "TRANSFER_ANNOTS --- GET_ALIGN_SEQS --- Error extracting sequences: %s", str(e))
        return None, None


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
    interval_key: str,
    conservation_start: int,
    conservation_end: int,
    logger: Optional[logging.Logger] = None,
) -> None:
    """Updates transfer_dict with conservation scores.

    Processes aligned sequences to map conserved positions and scores
    between target and conservation sequences.

    Args:
        transfer_dict: Transfer dictionary to update
        pfam_id: PFAM domain identifier
        target_name: Target sequence name
        target_seq: Target sequence
        conservation_seq: Conservation reference sequence
        conservation_key: Key (target_id/start-end) for conservation scores in conservations
        conservations: Whole conservation scores dictionary
        conserved_positions: List of conserved positions
        target_hit_start: Start position of the hit in the target sequence
        target_hit_end: End position of the hit in the target sequence
        interval_key: Key for current hit interval in the transfer_dict
        conservation_start: Start position in conservation reference sequence
        conservation_end: End position in conservation reference sequence
        logger: Logger object for registering debug messages
    """
    ranges_dict_key = "conservation_ranges"
    range_id = "conserved_positions"
    processed_conserved_pos = set()
    interval_data = transfer_dict[pfam_id]["sequence_id"][target_name]["hit_intervals"][interval_key]
    conservations_dict = interval_data["conservations"]

    for index, counter_cons_pos, counter_target_pos, char_cons, char_target in iterate_aligned_sequences(
        source_sequence=conservation_seq,
        target_sequence=target_seq,
        source_start=conservation_start,
        target_start=target_hit_start,
        source_end=conservation_end,
        target_end=target_hit_end
    ):
        if counter_target_pos is None or counter_cons_pos is None:
            continue

        if char_cons.islower() or char_target.islower():
            continue

        counter_cons_pos_str = str(counter_cons_pos)
        if counter_cons_pos_str in conserved_positions:
            processed_conserved_pos.add(counter_cons_pos_str)
            counter_target_pos_str = str(counter_target_pos)
            index_str = str(index)

            conserved_data = conservations[conservation_key][counter_cons_pos_str]
            conserved_amino = conserved_data["amino_acid"]
            score_data = conserved_data["conservation"]

            is_match = (char_target == conserved_amino and char_target not in '.-' and conserved_amino not in '.-')
            cons_type = "matches" if is_match else "misses"

            # Add to conservation_ranges-conserved_positions
            update_position_ranges(
                interval_data, counter_target_pos_str,
                ranges_dict_key, range_id
            )

            counter_target_pos_info = {
                "conservation": score_data,
                "residue": conserved_amino,
                "hit": is_match
            }

            conservations_dict["positions"][counter_target_pos_str] = counter_target_pos_info
            conservations_dict["indices"][cons_type].add(counter_target_pos_str)
            # Add position conversion info
            interval_data["position_conversion"]["target_to_aln"][counter_target_pos_str] = index_str
            aln_dict = interval_data["position_conversion"]["aln_to_target"]
            if index_str not in aln_dict:
                aln_dict[index_str] = counter_target_pos_str
        if counter_target_pos == target_hit_end or counter_cons_pos == conservation_end:
            break
    if logger:
        missing = set(conserved_positions) - processed_conserved_pos
        if missing:
            logger.debug("TRANSFER_ANNOTS --- POP_CONS --- Pfam ID: %s - Target Name: %s - Conserv-Rep Name: %s - Missing conserved positions: %s", pfam_id, target_name, conservation_key.split("/")[0], missing)

def populate_go_data_for_annotations(
    transfer_dict: dict,
    pfam_id: str,
    target_name: str,
    annotations: dict,
    go_terms_annot_key: str,
    target_gos: set
) -> None:
    """
    Adds GO data into transfer_dict for each anno_id found, using any overlapping GO terms
    between the target and conserved/annotated sequences to calculate similarity scores by Jaccard Index.

    Args:
        transfer_dict: Transfer dictionary to update
        pfam_id: PFAM domain identifier
        target_name: Target sequence name
        annotations: Annotations dictionary
        go_terms_annot_key: Key for GO terms in annotations
        target_gos: Set of GO terms found for the target sequence
    """
    go_info_cache = {}
    for interval_key in transfer_dict[pfam_id]["sequence_id"][target_name]["hit_intervals"]:
        interval_data = transfer_dict[pfam_id]["sequence_id"][target_name]["hit_intervals"][interval_key]
        for _, position_data in interval_data["annotations"]["positions"].items():
            for _, anno_data in position_data.items():
                go_data_for_annotation = {}
                for _, evidence_data in anno_data.get("evidence", {}).items():
                    annot_name = evidence_data.get("rep_mnemo_name")
                    if annot_name and annot_name in annotations and go_terms_annot_key in annotations[annot_name]:
                        cache_key = (target_name, annot_name)
                        if cache_key not in go_info_cache:
                            matched_go_info = {}
                            annotation_go_set = set()
                            for _, go_terms_with_meanings in annotations[annot_name][go_terms_annot_key].items():
                                for go_term, meaning in go_terms_with_meanings.items():
                                    annotation_go_set.add(go_term)
                                    if go_term in target_gos:
                                        matched_go_info[go_term] = meaning
                            target_go_set = set(target_gos)
                            intersection = len(target_go_set.intersection(annotation_go_set))
                            union = len(target_go_set.union(annotation_go_set))
                            jaccard_index = round(intersection / union, 4) if union > 0 else 0
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
    logger: logging.Logger,
    multi_logger: Callable,
    transfer_dict: dict,
    pfam_id: str,
    hmmalign_lines: list,
    conservations_filepath: str,
    annotations_filepath: str,
    output_dir: str,
    pfam_interpro_map_filepath: str
) -> dict:
    """Main function for enhancing transfer dictionary with conservation and GO data.

    Coordinates data population from conservations.json and annotations.json.
    Uses helper functions for sequence extraction from alignment,
    GO term retrieval from iprscan.tsv,
    and positional and GO term conservation scores calculation (latter 2 done by the populate_*).

    Args:
        logger: Process logging handler
        multi_logger: Callable for logging to multiple loggers
        transfer_dict: Dictionary to enhance
        pfam_id: Domain identifier
        hmmalign_lines: Alignment file content
        *_filepath: Paths to required JSON files
        output_dir: Directory for output files

    Returns:
        dict: Enhanced transfer dictionary with format:
            {"domain": {pfam_id: {...}}}
    """
    go_terms_annot_key = "0"

    if not transfer_dict:
        multi_logger("warning",
        "TRANSFER_ANNOTS --- CLEANUP_IMPROV_TD --- Empty transfer dictionary for Pfam ID %s - skipping data population", pfam_id)
        return {"domain": {"pfam_id": {}}}

    if "DOMAIN" in transfer_dict:
        transfer_dict[pfam_id] = transfer_dict.pop("DOMAIN")

    conservations, annotations = read_conservations_and_annotations(conservations_filepath, annotations_filepath)
    has_valid_conservations = bool(conservations) and conservations != {"sequence_id/range": {}} and any("/" in key for key in conservations.keys())
    has_valid_annotations = bool(annotations) and annotations != {"sequence_id": {}} and any(isinstance(annotations.get(key, {}).get("0", {}), dict) for key in annotations)
    if not has_valid_conservations and not has_valid_annotations:
        multi_logger("warning", "TRANSFER_ANNOTS --- CLEANUP_IMPROV_TD --- Both conservations and annotations data are empty or invalid - skipping data population")
        pfam_data = transfer_dict[pfam_id]
        return {"domain": {pfam_id: pfam_data}}

    mapping = pd.read_csv(pfam_interpro_map_filepath, sep="\t", header=0)
    matching_row = mapping.loc[mapping["Pfam_ID"] == pfam_id, "InterPro_ID"]

    if matching_row.empty:
        multi_logger("warning",
                    "TRANSFER_ANNOTS --- CLEANUP_IMPROV_TD --- No matching InterPro ID found for Pfam ID %s - proceeding regardless",
                    pfam_id)
        interpro_conv_id = ""
        logger.debug("TRANSFER_ANNOTS --- CLEANUP_IMPROV_TS --- Empty matching row: %s", matching_row)
    else:
        interpro_conv_id = matching_row.values[0]
        logger.debug("TRANSFER_ANNOTS --- CLEANUP_IMPROV_TS --- Found matching row: %s", matching_row)

    # For each target in the dictionary, gather GO terms, get aligned sequences and fill conservation and GO data
    for target_name in transfer_dict[pfam_id]["sequence_id"]:
        for interval_key in transfer_dict[pfam_id]["sequence_id"][target_name]["hit_intervals"]:
            target_hit_start = transfer_dict[pfam_id]["sequence_id"][target_name]["hit_intervals"][interval_key]["hit_start"]
            target_hit_end = transfer_dict[pfam_id]["sequence_id"][target_name]["hit_intervals"][interval_key]["hit_end"]
            target_id = f"{target_name}target//{target_hit_start}-{target_hit_end}"

            if has_valid_conservations:
                # Will always have 1 key in conservations,
                # a sequence from the original alignment used as a reference for the
                # conserved positions between the original and query/new alignments.
                conservation_key = list(conservations.keys())[0]
                conserved_positions = list(conservations[conservation_key].keys())

                # Extract sequences from hmmalign alignment
                target_seq, conservation_seq = get_alignment_sequences(
                    hmmalign_lines,
                    target_id,
                    conservation_key,
                    logger,
                    multi_logger
                )
                # Populate conservation data if sequences were found
                if target_seq and conservation_seq:
                    conservation_start = int(conservation_key.split("/")[1].split("-")[0])
                    conservation_end = int(conservation_key.split("/")[1].split("-")[1])
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
                        interval_key=interval_key,
                        conservation_start=conservation_start,
                        conservation_end=conservation_end,
                        logger=logger
                    )
            if has_valid_annotations:
                target_gos = gather_go_terms_for_target(
                    multi_logger, target_name, pfam_id,
                    output_dir, interpro_conv_id, target_hit_start, target_hit_end
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

    pfam_data = transfer_dict[pfam_id]
    return {"domain": {pfam_id: pfam_data}}

def convert_sets_and_tuples_to_lists(data):
    """Recursively converts sets and tuples to lists for JSON serialization.

    Special handling for annotation_ranges structure where:
        - positions (set) -> sorted list
        - ranges (list of tuples) -> list of lists

    Args:
        data: Input data structure containing sets/tuples

    Returns:
        Converted data structure with lists instead of sets/tuples
    """
    if isinstance(data, dict):
        # Special handling for annotation_ranges structure
        if "positions" in data and "ranges" in data:
            return {
                "positions": sorted(list(data["positions"])),
                "ranges": [list(r) for r in data["ranges"]]
            }
        return {k: convert_sets_and_tuples_to_lists(v) for k, v in data.items()}
    elif isinstance(data, tuple):
        return list(convert_sets_and_tuples_to_lists(elem) for elem in data)
    elif isinstance(data, list):
        return [convert_sets_and_tuples_to_lists(elem) for elem in data]
    elif isinstance(data, set):
        return sorted(list(data))
    else:
        return data

def write_reports(
    logger: logging.Logger,
    multi_logger: Callable,
    transfer_dict: dict,
    output_dir: str,
) -> None:
    """Writes transfer results to JSON files in two formats.

    Outputs transfer dictionary data as:
    1. Complete report: output_dir/pfam_id/pfam_id_report.json
       Contains all targets and their annotations for that domain.
    2. Per-target reports: output_dir/target_name/pfam_id_report.json
       Individual target-domain data.

    Args:
        logger: Process logging handler
        transfer_dict: Transfer results to be written
        output_dir: Base output directory

    Note:
        Converts set/tuple data to lists for JSON serialization
        while preserving original dictionary structure. Thi must be
        reversed for downstream consumption, use utils.convert_lists_to_original_types.
    """
    transfer_dict = convert_sets_and_tuples_to_lists(transfer_dict)

    # First check if we have an empty/invalid dictionary
    if not transfer_dict.get("domain"):
        multi_logger("warning", "TRANSFER_ANNOTS --- WRITE_REPORT --- Transfer dictionary is empty - no value assigned to 'domain' - skipping report writing")
        return

    # Get the pfam_id - either literal "pfam_id" or actual ID
    pfam_id = next(iter(transfer_dict["domain"].keys()))
    if transfer_dict["domain"][pfam_id] == {}:
        multi_logger("warning", "TRANSFER_ANNOTS --- WRITE_REPORT --- Transfer dictionary is empty - no value assigned to pfam_id %s - skipping report writing", pfam_id)
        return

    sequence_data = transfer_dict["domain"][pfam_id]["sequence_id"]

    # Write the entire transfer_dict to the pfam_id subdir
    pfam_dir = os.path.join(output_dir, pfam_id)
    os.makedirs(pfam_dir, exist_ok=True)
    entire_report_path = os.path.join(pfam_dir, pfam_id + "_report.json")

    structured_report = {
        "domain_id": pfam_id,
        "sequences": sequence_data
    }

    with open(entire_report_path, 'w', encoding="utf-8") as report_file:
        json.dump(structured_report, report_file, indent=4)
    logger.debug(f"TRANSFER_ANNOTS --- WRITE_REPORT --- Wrote entire Transfer Report: {entire_report_path}")

    # Write individual files for each target_name, preserving pfam_id structure
    for target_name, target_data in sequence_data.items():
        if target_data:
            sequence_dict = {
                "sequence_id": target_name,
                "domain": {
                    pfam_id: target_data
                }
            }
            logger.debug(f"TRANSFER_ANNOTS --- WRITE_REPORT --- Data to write for {target_name}!")
        else:
            sequence_dict = {
                "sequence_id": target_name,
                "domain": {
                    pfam_id: {"annotations": "None"}
                }
            }
            logger.debug(f"TRANSFER_ANNOTS --- WRITE_REPORT --- No data to write for {target_name}, mock !")

        safe_target_name = target_name.replace("|", "-")
        target_report_filepath = os.path.join(output_dir, safe_target_name, pfam_id + "_report.json")
        os.makedirs(os.path.join(output_dir, safe_target_name), exist_ok=True)

        with open(target_report_filepath, 'w', encoding="utf-8") as report_file:
            json.dump(sequence_dict, report_file, indent=4)

        logger.debug(f"TRANSFER_ANNOTS --- WRITE_REPORT --- Wrote Transfer Report for {target_name}-{pfam_id}: {target_report_filepath}")

#@measure_time_and_memory
##@profile
def map_and_filter_annot_pos(
    logger: logging.Logger,
    multi_logger: Callable,
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
    """Routes annotation position mapping between target and source/annotated sequences.

    Core routing function in annotation transfer pipeline. Handles two cases:
    1. Paired positions: Validates second member of paired annotations
    2. Single positions: Routes normal annotation validation

    Args:
        logger: Process logging handler
        multi_logger: Callable for logging to multiple loggers
        good_eco_codes: Valid evidence codes
        target_*: Target sequence data (sequence, name, boundaries)
        offset_*: Alignment boundaries (start, end)
        annot_*: Source sequence data
        entry_*: Source entry data (name, annotations)
        transfer_dict: Output dictionary for results
        processed_*: Set of handled annotations
        paired_*: Optional paired annotation data
        caller_target_pos_str: First position of pair in target

    Returns:
        For paired positions: (match_status, annotation_dict)
        For single positions: None
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
            processed_annotations=processed_annotations
        )
        return paired_tuple_res_hit_result_dict
    try:
        validate_annotations(
            logger=logger,
            multi_logger=multi_logger,
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
        logger.error(f"TRANSFER_ANNOTS ---> ERROR --- MAP_AND_FILTER --- Error in validate_annotations: {e}")
        traceback.print_exc()
        raise

def add_to_transfer_dict(
    hit: bool,
    multi_logger: Callable,
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
    Adds annotation data to transfer dictionary with position tracking.

    Core function for building transfer_dict structure. Handles both single and paired
    annotations, tracking positions, ranges, and relationships. Uses helper functions:
    - _add_single_annotation: Processes individual annotation entries
    - _update_annotation_ranges: Tracks continuous position ranges
    - _merge_adjacent_ranges: Combines overlapping ranges for storage

    Args:
        hit: Whether position matches between target/annotation
        logger: Process logging handler
        multi_logger: Callable for logging to multiple loggers
        transfer_dict: Output dictionary being built
        target_*: Target sequence information
        anno_*: Annotation data (id, total dict)
        entry_*: Source entry identifiers
        paired_*: Optional paired annotation data

    Structure built:
        transfer_dict["DOMAIN"]["sequence_id"][target_name]["hit_intervals"][interval] = {
            sequence data,
            annotations data,
            conservations data,
            position mappings,
            annotations/continuous ranges
        }
    """
    def extract_additional_keys(anno_total: dict, multi_logger: Callable) -> dict:
        try:
            return {
            k: v for k, v in anno_total.items()
            if k not in ["type", "description", "count", "evidence",
            "target_position", "annot_position", "paired_target_position",
            "annot_amino_acid", "target_amino_acid", "index_position",
            "gapped_paired", "insert_column_paired", "paired_target_pos_unreachable"]}
        except AttributeError as ae:
            multi_logger("error", "TRANSFER_ANNOTS --- ADD_TO_TRANSFER_DICT --- AttributeError for Anno_Total %s", anno_total)
            multi_logger("error", "TRANSFER_ANNOTS --- ADD_TO_TRANSFER_DICT --- AttributeError: %s - target_name: %s, entry_mnemo_name: %s\n %s",
                        ae, target_name, entry_mnemo_name, traceback.format_exc())
            raise

    additional_keys = extract_additional_keys(anno_total, multi_logger)
    paired_target_position_str = anno_total.get("paired_target_position", "0")

    # Initialize transfer_dict structure if needed
    transfer_dict.setdefault("DOMAIN", {})
    transfer_dict["DOMAIN"].setdefault("sequence_id", {})

    if target_name not in transfer_dict["DOMAIN"]["sequence_id"]:
        transfer_dict["DOMAIN"]["sequence_id"][target_name] = {
            "hit_intervals": {}
        }
    interval_key = f"{target_hit_start}-{target_hit_end}"

    if interval_key not in transfer_dict["DOMAIN"]["sequence_id"][target_name]["hit_intervals"]:
        transfer_dict["DOMAIN"]["sequence_id"][target_name]["hit_intervals"][interval_key] = {
            "sequence": target_sequence_continuous,
            "length": len(target_sequence_continuous),
            "hit_start": target_hit_start,
            "hit_end": target_hit_end,
            "annotations": {"positions": {}, "indices": {"matches": set(), "misses": set()}},
            "conservations": {"positions": {}, "indices": {"matches": set(), "misses": set()}},
            "position_conversion": {
                "target_to_aln": {},
                "aln_to_target": {}
            },
            "annotation_ranges": {},
            "conservation_ranges": {}
        }

    # Process main annotation
    _add_single_annotation(
        hit=hit,
        transfer_dict=transfer_dict,
        target_name=target_name,
        interval_key=interval_key,
        anno_id=anno_id,
        anno_total=anno_total,
        entry_mnemo_name=entry_mnemo_name,
        entry_primary_accession=entry_primary_accession,
        additional_keys=additional_keys
    )

    # Process paired annotation if provided
    if all(v is not None for v in [paired_anno_id, paired_anno_total]) and paired_target_position_str.isdigit():
        paired_additional_keys = extract_additional_keys(paired_anno_total, multi_logger)

        _add_single_annotation(
            hit=paired_position_res_hit,
            transfer_dict=transfer_dict,
            target_name=target_name,
            interval_key=interval_key,
            anno_id=paired_anno_id,
            anno_total=paired_anno_total,
            entry_mnemo_name=entry_mnemo_name,
            entry_primary_accession=entry_primary_accession,
            additional_keys=paired_additional_keys
        )

def _add_single_annotation(
    hit: bool,
    transfer_dict: dict,
    target_name: str,
    interval_key: str,
    anno_id: str,
    anno_total: dict,
    entry_mnemo_name: str,
    entry_primary_accession: str,
    additional_keys: dict) -> None:
    """Processes single annotation entry for add_to_transfer_dict.

    Handles position tracking, evidence logging, paired relationships,
    and special annotation types (BINDING, ACT_SITE, etc).
    Called for both main and second pair (if applicable) annotations.
    """
    ranges_dict_key = "annotation_ranges"
    # Extract positions in target and annotated numbering - for early pair data
    target_position_str = anno_total.get("target_position", None)
    annot_position_str = anno_total.get("annot_position", None)
    index_position_str = anno_total.get("index_position", None)
    paired_target_position_str = anno_total.get("paired_target_position", None)

    # Information to use in basic populating
    evidence_value = anno_total.get("evidence", None)
    type_value = anno_total.get("type", None)
    description_value = anno_total.get("description", None)
    count_value = anno_total.get("count", None)
    annot_amino_value = anno_total.get("annot_amino_acid", None)
    target_amino_value = anno_total.get("target_amino_acid", None)

    interval_dict = transfer_dict["DOMAIN"]["sequence_id"][target_name]["hit_intervals"][interval_key]

    # Initial setup for target position, including all references to it (annotations-positions,
    # annotations-indices-index_type, position_conversion-target_to_aln, position_conversion-aln_to_target)
    if target_position_str not in interval_dict["annotations"]["positions"]:
        interval_dict["annotations"]["positions"][target_position_str] = {}
        index_type = "matches" if hit else "misses"
        interval_dict["annotations"]["indices"][index_type].add(target_position_str)
        interval_dict["position_conversion"]["target_to_aln"][target_position_str] = index_position_str
        if index_position_str not in interval_dict["position_conversion"]["aln_to_target"]:
            interval_dict["position_conversion"]["aln_to_target"][index_position_str] = target_position_str

    positions_dict = interval_dict["annotations"]["positions"][target_position_str]

    # Populate essentials under anno_id for the position
    if anno_id not in positions_dict:
        essentials = {
            "type": type_value,
            "description": description_value,
            "count": count_value,
            "annot_amino_acid": annot_amino_value,
            "target_amino_acid": target_amino_value
        }
        positions_dict[anno_id] = {}
        positions_dict[anno_id].setdefault("essentials", essentials)
    else:
        positions_dict[anno_id]["essentials"]["count"] += 1

    positions_dict[anno_id].setdefault("hit", hit)

    update_position_ranges(
        interval_dict=interval_dict,
        target_position_str=target_position_str,
        ranges_dict_key=ranges_dict_key,
        range_id=anno_id
        )

    if evidence_value:
        if "evidence" not in positions_dict[anno_id]:
            positions_dict[anno_id].setdefault("evidence", {})
        if evidence_value not in positions_dict[anno_id]["evidence"]:
            positions_dict[anno_id]["evidence"][evidence_value] = {
                "rep_primary_accession": entry_primary_accession,
                "rep_mnemo_name": entry_mnemo_name,
                "count": 1
            }
        else:
            positions_dict[anno_id]["evidence"][evidence_value]["count"] += 1

    if "paired_target_position" in anno_total:
        if "paired_position" not in positions_dict[anno_id]:
            positions_dict[anno_id].setdefault("paired_position", {})
        if paired_target_position_str not in positions_dict[anno_id]["paired_position"]:
            pair_info = {
                "rep_primary_accession": entry_primary_accession,
                "rep_mnemo_name": entry_mnemo_name,
                "count": 1
            }
            # Possibly add gap or insert status to paired position
            if anno_total.get("gapped_paired", False):
                pair_info["gapped_target"] = True
            elif anno_total.get("insert_column_paired", False):
                pair_info["insert_column"] = True
            elif anno_total.get("paired_target_pos_unreachable", False):
                pair_info["unreachable_target"] = True
            positions_dict[anno_id]["paired_position"][paired_target_position_str] = pair_info
        else:
            positions_dict[anno_id]["paired_position"][paired_target_position_str]["count"] += 1

    # Helps analyze annotations that, in their description,
    # relate positions in annotated sequence uniprot numbering.
    # Note that those with sequence ranges are a bit tricky yet, but I"ll leave it as is for now
    if type_value in ["CROSSLNK", "DISULFID", "MUTAGEN"] and description_value and re.search(r"\b[A-Z]-\d+", description_value):
        additional_keys["annot_position"] = annot_position_str

    if additional_keys and type_value in ["BINDING", "ACT_SITE", "CROSSLNK", "DISULFID", "MUTAGEN"]:
        if "additional_keys" not in positions_dict[anno_id]:
            positions_dict[anno_id].setdefault("additional_keys", {})
        for key, value in additional_keys.items():
            if key not in positions_dict[anno_id]["additional_keys"]:
                positions_dict[anno_id]["additional_keys"][key] = {}
            if value not in positions_dict[anno_id]["additional_keys"][key]:
                positions_dict[anno_id]["additional_keys"][key][value] = {
                    "rep_primary_accession": entry_primary_accession,
                    "rep_mnemo_name": entry_mnemo_name,
                    "count": 1
                }
            else:
                positions_dict[anno_id]["additional_keys"][key][value]["count"] += 1

def update_position_ranges(interval_dict: dict, target_position_str: str, ranges_dict_key: str = "annotation_ranges", range_id: str = "") -> None:
    """Generic function for maintaining continuous position ranges.

    Args:
        interval_dict: Dictionary containing ranges data
        target_position_str: Position to process
        range_id: Identifier for the range group (annotation ID or conservation type)
        ranges_dict_key: Key for the ranges dictionary ("annotation_ranges" or "conservation_ranges")
    """
    position = int(target_position_str)
    ranges_dict = interval_dict[ranges_dict_key]

    if range_id not in ranges_dict and range_id:
        ranges_dict[range_id] = {
            "ranges": [(position, position)],
            "positions": {position}
        }
        return

    ranges = ranges_dict[range_id]["ranges"]
    positions = ranges_dict[range_id]["positions"]

    if position in positions:
        return

    positions.add(position)

    # Find appropriate range to merge with or create new
    for i, (start, end) in enumerate(ranges):
        if position == start - 1:
            ranges[i] = (position, end)
            merge_adjacent_ranges(ranges)
            return
        if position == end + 1:
            ranges[i] = (start, position)
            merge_adjacent_ranges(ranges)
            return
        if start <= position <= end:
            return

    # Position doesn't fit in existing ranges
    ranges.append((position, position))
    ranges.sort(key=lambda x: x[0])

def merge_adjacent_ranges(ranges: list) -> None:
    """Combines overlapping or consecutive position ranges.

    Helper for _update_annotation_ranges. Maintains sorted,
    non-overlapping ranges list by merging adjacent intervals.
    """
    if not ranges:
        return

    ranges.sort(key=lambda x: x[0])
    i = 0
    while i < len(ranges) - 1:
        current_end = ranges[i][1]
        next_start = ranges[i + 1][0]

        if current_end + 1 >= next_start:
            ranges[i] = (ranges[i][0], max(ranges[i][1], ranges[i + 1][1]))
            ranges.pop(i + 1)
        else:
            i += 1

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
    logger: Optional[logging.Logger] = None, # DEBUGGING, delete in production
    entry_annotations: Optional[Dict[str, Any]] = None,
    caller_target_pos_str: Optional[str] = None) -> Dict[str, Any]:
    """
    Processes single or paired annotations into a format requested by
    process_annotation() and validate_paired_annotations() and sent for use in
    add_to_transfer_dict(). Handles evidence validation and extracts special fields for BINDING/ACT_SITE types.

    Args:
        good_eco_codes: List of valid evidence codes to filter by
        entry_*: Source entry identifiers and annotations
        annotation_dict: Raw annotation data to process
        counter_*: Current positions in target/annotation sequences
        index: Position in alignment
        target_amino: Target sequence amino acid
        caller_target_pos_str: For paired annotations, caller's target position

    Returns:
        dict: Processed annotation data containing:
            - annotation: Original annotation dict
            - anno_type: Annotation type string
            - anno_id: Unique identifier (type + | + description)
            - anno_total: Processed data for add_to_transfer_dict()
            - paired_annot_pos_str: Paired position if applicable
    """
    anno_type = ""
    anno_id = ""
    anno_total = None
    annotation = annotation_dict

    # DEBUGGING, still keeping logger and entry_annotations for now, will DELETE in PRODUCTION!
    logger.debug(f"TRANSFER_ANNOTS --- MAKE_ANNO --- Annotation Dict: {annotation_dict}")
    logger.debug(f"TRANSFER_ANNOTS --- MAKE_ANNO --- Entry Annotations List: {entry_annotations.get(counter_annot_pos_str)}")
    logger.debug(f"TRANSFER_ANNOTS --- MAKE_ANNO --- Good ECO codes : {good_eco_codes}")

    paired_annot_pos_str = annotation.get("paired_position", None)
    anno_type = annotation["type"]
    anno_desc = annotation.get("description", None)
    anno_id = f"{anno_type} | {anno_desc}"
    anno_count = 1
    anno_evidence = {entry_mnemo_name: annotation.get("evidence", None)}
    annot_amino_acid = annotation.get("aminoacid", None)

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
            "target_amino_acid": target_amino,
            # Paired annotation flags
            "gapped_paired": False,
            "insert_column_paired": False,
            "paired_target_pos_unreachable": False
        }
        if paired_annot_pos_str and caller_target_pos_str:
            anno_total["paired_target_position"] = caller_target_pos_str
        # Changed to include ACT_SITE, because of ligand_ccd_id in BioLiP sourced annotations
        if anno_type in ["BINDING", "ACT_SITE"]:
            keys_to_exclude = {"type", "description", "count", "evidence",
                               "paired_position", "aminoacid", "entry"}
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
    multi_logger: Callable,
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
    Directs annotations to preview if we should add them to transfer_dict by calling make_anno_total_dict()
    returning anno_total containing type, description, count, evidence, target/annot/index positions
    and target/annot amino acids data. If so, call add_to_transfer_dict().
    We may also add paired position information in much the same way, if it hits
    (map_and_filter_annot_pos() attempts to validate_paired_positions()).

    Args:
        logger: Process logging handler
        multi_logger: Callable for logging to multiple loggers
        good_eco_codes: Valid evidence codes
        target_*: Target data (name, sequence, positions, amino)
        entry_*: Source entry data (name, annotations)
        annotation_*: Current annotation info (dict, key)
        transfer_dict: Output dictionary for results
        offset_*: Alignment boundaries
        counter_*: Position strings being processed
        index: Current position in alignment
        processed_annotations: Set of handled annotations
        res_hit: Whether current position matches
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
    annotation = result_dict["annotation"]
    anno_type = result_dict["anno_type"]
    anno_id = result_dict["anno_id"]
    anno_total = result_dict["anno_total"]
    paired_annot_pos_str = result_dict["paired_annot_pos_str"]

    if anno_total:
        entry_primary_accession = annotation.get("entry", None)
        target_sequence_continuous = ''.join([char.upper() for char in target_sequence if char.isalpha()])
        if anno_type in ["DISULFID", "CROSSLNK", "SITE", "BINDING"] and paired_annot_pos_str is not None:
            paired_annotations = entry_annotations.get(paired_annot_pos_str, [])
            paired_annotation_dict = next(
                (annotation for annotation in paired_annotations if annotation["type"] == anno_type), None
            )
            if paired_annotation_dict:
                paired_position_res_hit, paired_result_dict = map_and_filter_annot_pos(
                    logger=logger,
                    multi_logger=multi_logger,
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
                    paired_annot_pos_str=paired_annot_pos_str,
                    caller_target_pos_str=counter_target_pos_str,
                    paired_annotation_dict=paired_annotation_dict
                )
                if paired_result_dict.get("gapped_paired", False):
                    # Gapped paired position, will only add caller position
                    # in add_to_transfer_dict(), mentioning the missing paired
                    anno_total["gapped_paired"] = True
                    paired_anno_total = None
                    paired_anno_id = None
                    paired_target_position_str = "0"
                elif paired_result_dict.get("insert_column_paired", False):
                    # Insert column in paired position, will only add caller position
                    # in add_to_transfer_dict(), mentioning the missing paired
                    anno_total["insert_column_paired"] = True
                    paired_anno_total = None
                    paired_anno_id = None
                    paired_target_position_str = "0"
                elif paired_result_dict.get("paired_target_pos_unreachable", False):
                    # A position counter hit target sequence's end before a match state column
                    # was found, so we'll add caller position and mention the missing paired
                    anno_total["paired_target_pos_unreachable"] = True
                    paired_anno_total = None
                    paired_anno_id = None
                    paired_target_position_str = "0"
                else:
                    paired_anno_total = paired_result_dict.get("anno_total", None)
                    paired_anno_id = paired_result_dict["anno_id"]
                    paired_target_position_str = paired_anno_total.get("target_position", "0")
            else:
                paired_position_res_hit = False
                paired_anno_total = None
                paired_anno_id = None
                paired_target_position_str = "0"
            anno_total["paired_target_position"] = paired_target_position_str
            logger.debug(
                f"TRANSFER_ANNOTS --- PROCESS_ANNOT --- Paired Position Valid (Present in Entry Annotations) with "
                f"target---caller_pos-hit_status-paired_pos-hit_status {target_name}---{counter_target_pos_str}-{res_hit}-"
                f"{paired_target_position_str}-{paired_position_res_hit} and "
                f"annotated {entry_mnemo_name}"
            )
            try:
                add_to_transfer_dict(
                    hit=res_hit,
                    multi_logger=multi_logger,
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
                if anno_total["gapped_paired"]:
                    paired_annotation_key = (
                        entry_mnemo_name,
                        target_name,
                        paired_annot_pos_str,
                        paired_target_position_str,
                        anno_type,
                        "gapped_target_position"
                    )
                elif anno_total["insert_column_paired"]:
                    paired_annotation_key = (
                        entry_mnemo_name,
                        target_name,
                        paired_annot_pos_str,
                        paired_target_position_str,
                        anno_type,
                        "insert_column"
                    )
                elif anno_total["paired_target_pos_unreachable"]:
                    paired_annotation_key = (
                        entry_mnemo_name,
                        target_name,
                        paired_annot_pos_str,
                        paired_target_position_str,
                        anno_type,
                        "paired_target_pos_unreachable"
                    )
                else:
                    paired_annotation_key = (
                        entry_mnemo_name,
                        target_name,
                        paired_annot_pos_str,
                        paired_target_position_str,
                        anno_type,
                        "match_column"
                    )
                processed_annotations.add(paired_annotation_key)
                log_base = (
                    f"TRANSFER_ANNOTS --- PROCESS_ANNOT --- SUCCESS - ROUTE: Paired - "
                    f"Added to transfer_dict for target {target_name} and "
                    f"annotated {entry_mnemo_name} at target_caller {counter_target_pos_str}"
                )
                if paired_anno_total:
                    logger.debug(
                        f"{log_base} with paired {paired_target_position_str} --- "
                        f"Caller-paired positions: {counter_annot_pos_str}-{paired_anno_total['annot_position']} --- "
                        f"ORIGIN: Type {anno_type} and Evidence {anno_total['evidence']}"
                    )
                else:
                    logger.debug(
                        f"{log_base} --- No valid paired annotation found at position {paired_annot_pos_str} --- "
                        f"ORIGIN: Type {anno_type} and Evidence {anno_total['evidence']}"
                    )
            except Exception as e:
                multi_logger("error",
                "TRANSFER_ANNOTS --- PROCESS_ANNOT --- FAILURE - ROUTE: Paired - Error in add_to_transfer_dict: %s", e)
                traceback.print_exc()
                raise
        else:
            add_to_transfer_dict(
                hit=res_hit,
                multi_logger=multi_logger,
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
                f"TRANSFER_ANNOTS --- PROCESS_ANNOT --- SUCCESS - ROUTE: Single - "
                f"Added to transfer_dict for target {target_name} and "
                f"annotated {entry_mnemo_name} at target {counter_target_pos_str} and "
                f"annotated {counter_annot_pos_str} --- ORIGIN: Type {anno_type} and "
                f"Evidence {anno_total['evidence']}"
            )
    else:
        logger.debug(
            f"TRANSFER_ANNOTS --- PROCESS_ANNOT --- FAILURE - ROUTE: Single - "
            f"Anno Total was None for target {target_name} and annotated {entry_mnemo_name} "
            f"at target {counter_target_pos_str} and annotated {counter_annot_pos_str} "
            f"--- ORIGIN: No evidence of desired type"
        )
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
    processed_annotations: set
) -> tuple[bool, dict]:
    """
    Validates the 2nd member of a paired position by checking
    if the target paired annot position from caller is a valid match.
    Success and failure both mean adding the annotation to the transfer dictionary
    in the process_annotation parent function, difference lies in where (match | miss).
    The caller_target_pos_str is the first member of the pair in the target sequence,
    and, if it was a gap, the paired result dict will have a key "gapped_paired" set to True,
    used in adding to transfer_dict to track this failure.

    Args:
        logger: Process logging handler
        good_eco_codes: Valid evidence codes
        target_*: Target sequence data (sequence, name, hit boundaries)
        offset_*: Alignment boundaries (start, end)
        annot_*: Source annotation data (sequence, entry name)
        paired_*: Paired position info (annotation dict, position str)
        entry_annotations: Position-based annotation dictionary
        counter_*: Current positions in sequences
        caller_target_pos_str: ideally first position (or second if the former was a gap in target_sequence) of pair in target_sequence
        processed_annotations: Set of handled annotations

    Returns:
        tuple: (match_status, paired_result_dict)
    """
    earlier_annotation_key_gapped = (
        entry_mnemo_name,
        target_name,
        paired_annot_pos_str,
        "0",
        paired_annotation_dict["type"],
        "gapped_target_position"
    )

    earlier_annotation_key_insert = (
        entry_mnemo_name,
        target_name,
        paired_annot_pos_str,
        "0",
        paired_annotation_dict["type"],
        "insert_column"
    )

    has_gap_in_first_pos = earlier_annotation_key_gapped in processed_annotations
    first_pos_in_insert_col = earlier_annotation_key_insert in processed_annotations

    paired_result_dict = {
        "gapped_paired": has_gap_in_first_pos,
        "insert_column_paired": first_pos_in_insert_col,
        "paired_target_pos_unreachable": False
    }

    paired_position_res_hit = False

    # Paired pos. came earlier than caller, and paired was a gap in target_sequence
    # or located in an insert column - Sent to processed_annotations(),
    # and previously flagged in validate_annotations()
    if has_gap_in_first_pos or first_pos_in_insert_col:
        return paired_position_res_hit, paired_result_dict

    paired_annot_pos_int = int(paired_annot_pos_str)

    logger.debug(
        f"TRANSFER_ANNOTS --- VAL_PAIRED --- Running for target "
        f"{target_name}/{target_hit_start}-{target_hit_end}"
        f"and annotated {entry_mnemo_name}/{offset_start}-{offset_end} -"
        f"Target Caller Position {caller_target_pos_str}"
        f"and Paired Annot Position {paired_annot_pos_str}"
    )

    for index, counter_annot_pos, counter_target_pos, char_annot, char_target in iterate_aligned_sequences(
        source_sequence=annot_sequence,
        target_sequence=target_sequence,
        source_start=offset_start,
        target_start=target_hit_start,
        source_end=offset_end,
        target_end=target_hit_end
        ):
        if counter_annot_pos is None:
            continue

        if counter_annot_pos == paired_annot_pos_int:
            if char_annot.islower() or char_target.islower():
                paired_result_dict["insert_column_paired"] = True
                return paired_position_res_hit, paired_result_dict

        if counter_annot_pos == paired_annot_pos_int and counter_target_pos is None:
            paired_result_dict["gapped_paired"] = True
            return paired_position_res_hit, paired_result_dict

        if counter_annot_pos == paired_annot_pos_int and counter_target_pos is not None:
            counter_target_pos_str = str(counter_target_pos)

            # DEBUGGING INFO
            logger.debug("\n TRANSFER_ANNOTS--- VAL_PAIRED --- Helpful Info for paired match | miss \n")
            start_index = max(0, index - 3)
            end_index = min(len(annot_sequence), index + 4)
            annot_window = annot_sequence[start_index:end_index]
            target_window = target_sequence[start_index:end_index]
            logger.debug(f"TRANSFER_ANNOTS --- VAL_PAIRED --- Counter Tar Pos {counter_target_pos_str} and amino acid {target_sequence[index]} + Counter annot Pos (annotated) {str(counter_annot_pos)} and amino acid {annot_sequence[index]}")
            logger.debug(f"TRANSFER_ANNOTS --- VAL_PAIRED --- Annot Window: {annot_window} + Target Window: {target_window}")

            paired_target_amino = target_sequence[index]
            paired_position_res_hit = bool(char_target == char_annot)

            paired_result_dict.update(make_anno_total_dict(
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
            ))

            return paired_position_res_hit, paired_result_dict

        if counter_target_pos == target_hit_end or counter_annot_pos == offset_end:
            break

    paired_result_dict["paired_target_pos_unreachable"] = True
    return paired_position_res_hit, paired_result_dict

#@measure_time_and_memory
#@profile
def validate_annotations(
    logger: logging.Logger,
    multi_logger: Callable,
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
    Regardless of hit, send to process_annotation, granted that processed_annotations does not contain the annotation key (entry_mnemo_name, target_name, counter_annot_pos_str, counter_target_pos_str and anno_type), avoiding reprocessing.

    Args:
        logger: Logger for process tracking
        multi_logger: Callable for logging to multiple loggers
        good_eco_codes: List of valid evidence codes
        target_*: Target sequence info (sequence, name, hit_start/end)
        offset_*: Alignment offset positions (start/end)
        annot_*: Annotation source info (sequence, entry_mnemo_name)
        entry_annotations: Dictionary of position-based annotations
        transfer_dict: Output dictionary for storing results
        processed_annotations: Set of already handled annotations
        counter_*: Current position counters (target/annot)
    """
    annotated_annot_pos_list = [int(key) for key in entry_annotations.keys() if key != "0"]
    if not any(pos in range(offset_start, offset_end + 1) for pos in annotated_annot_pos_list):
        logger.debug(f"TRANSFER_ANNOTS --- VAL_ANNOTS --- No annotations found for target {target_name} and annotated {entry_mnemo_name} in offset range {offset_start}-{offset_end}")
        return
    failure_target_pos = "0"

    for index, counter_annot_pos, counter_target_pos, char_annot, char_target in iterate_aligned_sequences(
        source_sequence=annot_sequence,
        target_sequence=target_sequence,
        source_start=offset_start,
        target_start=target_hit_start,
        source_end=offset_end,
        target_end=target_hit_end):

        if counter_target_pos == target_hit_end and counter_annot_pos:
            # Check if there are any remaining paired annotations
            remaining_positions = [
                pos for pos in entry_annotations.keys()
                if pos != "0" and int(pos) > counter_annot_pos
            ]
            for pos in remaining_positions:
                for annotation_dict in entry_annotations[pos]:
                    anno_type = annotation_dict.get("type", None)
                    paired_annot_pos_str = annotation_dict.get("paired_position", None)
                    if anno_type in ["DISULFID", "CROSSLNK", "SITE", "BINDING"] and paired_annot_pos_str is not None:
                        unreachable_key = (
                            entry_mnemo_name,
                            target_name,
                            pos,
                            failure_target_pos,
                            anno_type,
                            "paired_target_pos_unreachable"
                        )
                        processed_annotations.add(unreachable_key)
            break

        if counter_annot_pos is None:
            continue

        counter_annot_pos_str = str(counter_annot_pos)

        if char_annot.islower() or char_target.islower():
            if counter_annot_pos_str in entry_annotations:
                for annotation_dict in entry_annotations[counter_annot_pos_str]:
                    anno_type = annotation_dict.get("type", None)
                    paired_annot_pos_str = annotation_dict.get("paired_position", None)
                    if anno_type in ["DISULFID", "CROSSLNK", "SITE", "BINDING"] and paired_annot_pos_str is not None:
                        insert_col_annotation_key = (
                            entry_mnemo_name,
                            target_name,
                            counter_annot_pos_str,
                            failure_target_pos,
                            anno_type,
                            "insert_column"  # Field to track insert col. failures
                        )
                        processed_annotations.add(insert_col_annotation_key)
                continue

        if counter_annot_pos_str in entry_annotations and counter_target_pos is None:
            for annotation_dict in entry_annotations[counter_annot_pos_str]:
                anno_type = annotation_dict.get("type", None)
                paired_annot_pos_str = annotation_dict.get("paired_position", None)
                if anno_type in ["DISULFID", "CROSSLNK", "SITE", "BINDING"] and paired_annot_pos_str is not None:
                    gapped_annotation_key = (
                        entry_mnemo_name,
                        target_name,
                        counter_annot_pos_str,
                        failure_target_pos,
                        anno_type,
                        "gapped_target_position" # Field to track gapped target failures
                    )
                    processed_annotations.add(gapped_annotation_key)
            continue

        if counter_annot_pos_str in entry_annotations:
            counter_target_pos_str = str(counter_target_pos)

            # DEBUGGING INFO
            logger.debug(f"TRANSFER_ANNOTS--- VAL_ANNOTS --- Helpful Info for Single/Caller Annotation Match \n NAMES: target {target_name} and annot {entry_mnemo_name} \n POSITIONS: target {counter_target_pos_str} and annot {counter_annot_pos_str} \n")
            start_index = max(0, index - 3)
            end_index = min(len(annot_sequence), index + 4)
            annot_window = annot_sequence[start_index:end_index]
            target_window = target_sequence[start_index:end_index]
            logger.debug(f"TRANSFER_ANNOTS --- VAL_ANNOTS --- Counter annot Pos {counter_annot_pos_str} and amino acid {annot_sequence[index]} + Counter Tar Pos {str(counter_target_pos)} and amino acid {target_sequence[index]}")
            logger.debug(f"TRANSFER_ANNOTS --- VAL_ANNOTS --- Annot Window: {annot_window} + Target Window: {target_window}")

            target_amino = target_sequence[index]
            res_hit = bool(char_target == char_annot)
            for annotation_dict in entry_annotations[counter_annot_pos_str]:
                anno_type = annotation_dict.get("type", None)
                match_annotation_key = (
                    entry_mnemo_name,
                    target_name,
                    counter_annot_pos_str,
                    counter_target_pos_str,
                    anno_type,
                    "match_column" # Field to track match successes
                )
                if match_annotation_key in processed_annotations:
                    continue
                try:
                    process_annotation(
                        res_hit=res_hit,
                        logger=logger,
                        multi_logger=multi_logger,
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
                        annotation_key=match_annotation_key
                    )
                except Exception as e:
                    multi_logger("error", "TRANSFER_ANNOTS --- VAL_ANNOTS --- Downstream error in process_annotation: %s", e)
                    traceback.print_exc()
                    raise


def main():
    """Main function, initializes this script"""
    args = parse_arguments()
    dom_align = args.dom_align
    resource_dir = args.resource_dir
    output_dir = args.output_dir
    good_eco_codes = args.eco_codes
    pfam_interpro_map_filepath = os.path.join(resource_dir, "mappings/interpro_pfam_accession_mapping.tsv")
    main_logger, _ = get_logger(args.log, scope="main")
    domain_logger, _ = get_logger(args.log, scope="domain", identifier=args.domain_accession)
    multi_logger = get_multi_logger([main_logger, domain_logger])

    domain_logger.info("TRANSFER_ANNOTS --- MAIN --- Running transfer_annotations.py for %s", dom_align)

    pfam_id = get_pfam_id_from_hmmalign_result(dom_align)
    annotations_filepath, conservations_filepath = get_annotation_filepath(resource_dir, pfam_id)
    hmmalign_lines, annotations = read_files(dom_align, annotations_filepath)

    try:
        if annotations == {"sequence_id": {}}:
            # CONSERVATIONS ONLY PATH
            domain_logger.info("TRANSFER_ANNOTS --- MAIN --- No annotations file found - proceeding with conservations-only mode")
            transfer_dict = setup_for_conservations_only(domain_logger, multi_logger, hmmalign_lines, pfam_id)
        else:
            # ANNOTATIONS + CONSERVATIONS PATH
            domain_logger.debug("TRANSFER_ANNOTS --- MAIN --- Anno. + Cons. mode - Good ECO Codes to Filter by %s", good_eco_codes)
            transfer_dict = find_and_map_annots(domain_logger, multi_logger, hmmalign_lines, annotations, good_eco_codes)
        if not transfer_dict:
            domain_logger.info("TRANSFER_ANNOTS --- MAIN --- Transfer Dict was EMPTY")
        else:
            domain_logger.info("TRANSFER_ANNOTS --- MAIN --- Transfer Dict FILLED")

    except (KeyError, IndexError, AttributeError) as e:
        error_info = traceback.format_exc()
        multi_logger("error", "TRANSFER_ANNOTS --- MAIN --- ERROR transferring annotations for Pfam ID %s: %s\n%s", pfam_id, e, error_info)
        raise

    improved_transfer_dict = cleanup_improve_transfer_dict(
        domain_logger, multi_logger, transfer_dict, pfam_id, hmmalign_lines,
        conservations_filepath, annotations_filepath, output_dir, pfam_interpro_map_filepath
        )
    write_reports(domain_logger, multi_logger, improved_transfer_dict, output_dir)

if __name__ == "__main__":
    main()
