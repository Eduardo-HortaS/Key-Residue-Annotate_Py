"""
transfer_annotations.py

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

This script contains all functions dealing directly with transferring annotations.
Intended usage begins by calling find_and_map_annotations with 2 arguments:
a list of hmmalign result lines and the loaded annotations dict.

"""

import os
import json
import argparse
import logging
import traceback
import copy
from typing import Any, Dict, Optional, Set
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
    parser.add_argument("-iA", "--dom-align", help="Path to domain's hmmalign alignment", required=True, type=str)
    parser.add_argument("-r", "--resource-dir", help="Resource dir path", required=True, type=str)
    parser.add_argument("-o", "--output-dir", help="Output dir path", required=True, type=str)
    parser.add_argument("-e", "--eco-codes", help="List of ECO codes to filter annotations", nargs="*", type=str, default=[])
    parser.add_argument("-l", "--log", help="Log path", \
        required=False, type=str, default="logs/transfer_annotations.log")
    return parser.parse_args()

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
    Constructs and returns the filepaths for annotations,
    seed alignment, and report based on the PFAM ID.
    """
    annotations_filepath = os.path.join(resource_dir, pfam_id, "annotations.json")
    return annotations_filepath

def read_files(hmmalign_result: str, annotations_filepath: str) -> tuple[list[str], Any]:
    """
    Reads and returns the content of the hmmalign result and annotations files.
    """
    with open(hmmalign_result, 'r', encoding="utf-8") as hmmaligned_file:
        hmmalign_lines = hmmaligned_file.readlines()
    with open(annotations_filepath, 'r', encoding="utf-8") as annotations_file:
        annotations = json.load(annotations_file)
    return hmmalign_lines, annotations

#@measure_time_and_memory
#@profile
def find_and_map_annots(logger: logging.Logger, hmmalign_lines: list, annotations: dict, good_eco_codes: list) -> dict:
    """
    Finds target lines and stores data for them in target_info.
    Then, finds sequence lines in hmmalign alignment for which we have annotations
    and calls map_and_filter_annot_pos to map positions with target lines iteratively.
    Makes visited_pos set to keep track of visited positions across calls
    and sends it to map_and_filter_annot_pos.
    """
    transfer_dict = {}
    target_info = {}
    visited_pos = {}
    for line in hmmalign_lines:
        if "target/" in line and not line.startswith("#"):
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
                        if target_name not in visited_pos:
                            visited_pos[target_name] = {}
                        if entry_mnemo_name not in visited_pos[target_name]:
                            visited_pos[target_name][entry_mnemo_name] = {}
                        if target_hit_interval not in visited_pos[target_name][entry_mnemo_name]:
                            visited_pos[target_name][entry_mnemo_name][target_hit_interval] = set()
                        for target_hit_start, target_hit_end, target_hit_sequence in target_info_list:
                            entry_annotations_copy = copy.deepcopy(entry_annotations)
                            # DEBUGGING a specific pair?
                            # if target_name == "sp||PTPRJ_HUMAN" and entry_mnemo_name == "TIE2_HUMAN":
                            map_and_filter_annot_pos(logger, good_eco_codes, target_hit_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_mnemo_name, entry_annotations_copy, transfer_dict, visited_pos[target_name][entry_mnemo_name][target_hit_interval])
                            logger.info(f"Mapped, filtered and possibly added to transfer_dict: target {target_name} and annotated {entry_mnemo_name} at target hit interval {target_hit_interval}")
    return transfer_dict

def write_report(logger: logging.Logger, transfer_dict: dict, pfam_id: str, output_dir: str) -> None:
    """
    Writes the transfer dictionary to a report file in JSON format, per sequence-domain pair.
    The report is written in the sequence's directory, with the name PF*_report.json.
    """
    for target_name in transfer_dict.keys():
        if transfer_dict[target_name]:
            json_data = json.dumps(transfer_dict[target_name], indent=4)
            logger.info(f"---> DEBUG --- WRITE_REPORT --- We got data to write for {target_name}!")
        else:
            json_data = json.dumps({"Annotations": "None"}, indent=4)
            logger.info(f"---> DEBUG --- WRITE_REPORT --- No data to write for {target_name}!")
        target_name = target_name.replace("|", "-")
        pfam_id_report_filepath = os.path.join(output_dir, target_name, pfam_id + "_report.json")
        with open(pfam_id_report_filepath, 'w', encoding="utf-8") as report_file:
            report_file.write(json_data)
        logger.info(f"---> DEBUG --- WRITE_REPORT ---  Wrote Transfer Report for sequence-domain pair, {target_name}-{pfam_id}: {pfam_id_report_filepath}")

#@measure_time_and_memory
##@profile
def map_and_filter_annot_pos(logger: logging.Logger, good_eco_codes: list, target_sequence: str, target_name: str, target_hit_start: int, target_hit_end: int, offset_start: int, offset_end: int, annot_sequence: str, entry_mnemo_name: str, entry_annotations: Dict[str, Any], transfer_dict: dict, visited_pos: Set[str], target_paired_uniprot_pos: str =None, caller_target_pos: int =None, annotation_dict: Optional[Dict[str, Any]] = None) -> (bool | None):
    """
    Asks to map positions between annotated and target sequences.
    If target paired uniprot position and caller target position are provided,
    this is an associated call for a paired position sent to validate_paired_positions.
    Otherwise, it is a call for a single position sent to validate_annotations.
    Makes counter_uniprot_pos and counter_target_pos to track the current position in both sequences.
    """
    counter_uniprot_pos = None
    counter_target_pos = None

    if target_paired_uniprot_pos and caller_target_pos and annotation_dict:
        return validate_paired_positions(logger, good_eco_codes, target_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_mnemo_name, annotation_dict, entry_annotations, transfer_dict, visited_pos, counter_target_pos, counter_uniprot_pos, target_paired_uniprot_pos, caller_target_pos)
    try:
        return validate_annotations(logger, good_eco_codes, target_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_mnemo_name, entry_annotations, transfer_dict, visited_pos, counter_target_pos, counter_uniprot_pos)
    except Exception as e:
        logger.error(f"---> DEBUG --- MAP_AND_FILTER --- Error in validate_annotations: {e}")
        traceback.print_exc()
        raise

#@measure_time_and_memory
#@profile
def add_to_transfer_dict(logger: logging.Logger, transfer_dict: dict, target_name: str, counter_target_pos: int, anno_id: str, anno_total: dict, entry_mnemo_name: str, entry_alphanum_name: str) -> None:
    """
    Adds anno_total data to the transfer dictionary.
    The dictionary is structured by target position, annotation ID (type + description),
    and various possible additional annotation details, mainly for BINDING annotations.
    """
    try:
        additional_keys = {k: v for k, v in anno_total.items() if k not in ['type', 'description', 'count', 'evidence', 'paired_position', 'aminoacid', 'target_position']}
    except AttributeError as ae:
        logger.error("@ @ ------------- AttributeError: Anno_total: %s", anno_total)
        logger.error(f"@ @-------------- AttributeError: {ae} - target_name: {target_name}, entry_mnemo_name: {entry_mnemo_name} \n", traceback.format_exc())
        raise
    entry_codes = entry_alphanum_name + " | " + entry_mnemo_name

    if target_name not in transfer_dict:
        transfer_dict[target_name] = {}

    if counter_target_pos not in transfer_dict[target_name]:
        transfer_dict[target_name][counter_target_pos] = {}

    if anno_id not in transfer_dict[target_name][counter_target_pos]:
        essentials = {
            'type': anno_total.get('type'),
            'description': anno_total.get('description'),
            'count': anno_total.get('count')
        }
        transfer_dict[target_name][counter_target_pos][anno_id] = {}
        transfer_dict[target_name][counter_target_pos][anno_id].setdefault('essentials', essentials)
        evidence_value = anno_total.get('evidence', None)
        paired_position_value = anno_total.get('paired_position', None)
    else:
        transfer_dict[target_name][counter_target_pos][anno_id]['essentials']['count'] += 1
        evidence_value = anno_total.get('evidence', None)
        paired_position_value = anno_total.get('paired_position', None)

    if evidence_value:
        if 'evidence' not in transfer_dict[target_name][counter_target_pos][anno_id]:
            transfer_dict[target_name][counter_target_pos][anno_id].setdefault('evidence', {})
        if evidence_value not in transfer_dict[target_name][counter_target_pos][anno_id]['evidence']:
            transfer_dict[target_name][counter_target_pos][anno_id]['evidence'][evidence_value] = {
                "rep_entry_codes": entry_codes, "count": 1}
        else:
            transfer_dict[target_name][counter_target_pos][anno_id]['evidence'][evidence_value]["count"] += 1

    if paired_position_value:
        if 'paired_position' not in transfer_dict[target_name][counter_target_pos][anno_id]:
            transfer_dict[target_name][counter_target_pos][anno_id].setdefault('paired_position', {})
        if paired_position_value not in transfer_dict[target_name][counter_target_pos][anno_id]['paired_position']:
            transfer_dict[target_name][counter_target_pos][anno_id]['paired_position'][paired_position_value] = {
                "rep_entry_codes": entry_codes, "count": 1}
        else:
            transfer_dict[target_name][counter_target_pos][anno_id]['paired_position'][paired_position_value]["count"] += 1

    if additional_keys:
        if 'additional_keys' not in transfer_dict[target_name][counter_target_pos][anno_id]:
            transfer_dict[target_name][counter_target_pos][anno_id].setdefault('additional_keys', {})
        for key, value in additional_keys.items():
            if anno_total['type'] == 'BINDING':
                if key not in transfer_dict[target_name][counter_target_pos][anno_id]['additional_keys']:
                    transfer_dict[target_name][counter_target_pos][anno_id]['additional_keys'][key] = {
                        "rep_entry_codes": entry_codes, "count": 1, "value": value}
                else:
                    transfer_dict[target_name][counter_target_pos][anno_id]['additional_keys'][key]["count"] += 1

#@measure_time_and_memory
#@profile
def remove_failed_annotations(entry_annotations: Dict[str, Any], counter_uniprot_pos_str: str, paired_position: str, annotation_type: str) -> None:
    """
    Removes the specific failed annotations from entry annotations.
    """
    # Remove the specific annotation type from the counter_uniprot_pos_str
    entry_annotations[counter_uniprot_pos_str] = [
        annotation for annotation in entry_annotations[counter_uniprot_pos_str]
        if annotation['type'] != annotation_type
    ]

    # Remove the specific annotation type from the paired_position
    if paired_position in entry_annotations and entry_annotations[paired_position]:
        entry_annotations[paired_position] = [
            annotation for annotation in entry_annotations[paired_position]
            if annotation['type'] != annotation_type
        ]


#@measure_time_and_memory
#@profile
def make_anno_total_dict(logger: str, good_eco_codes: list, entry_mnemo_name: str, annotation_dict: Dict[str, Any], entry_annotations: Dict[str, Any], counter_target_pos: int, counter_uniprot_pos_str: str, visited_pos: Set[str], caller_target_pos: Optional[int] = None) -> Dict[str, Any]:
    """
    For any given annotation in appropriate format (list with 1 dict inside), paired or not, produces an anno_total dict and
    associated data to return as a dict to be extracted from in accordance to the calling function needs.
    Used by process_annotation and validate_paired_positions.
    """
    annotation = {}
    anno_type = ""
    anno_id = ""
    anno_total = None
    paired_position = None
    # DEBUGGING, still keeping logger and entry_annotations for now, will DELETE in PRODUCTION!
    # logger.info(f"---> DEBUG --- MAKE_ANNO --- Annotation Dict: {annotation_dict}")
    # logger.info(f"---> DEBUG --- MAKE_ANNO --- Entry Annotations List: {entry_annotations.get(counter_uniprot_pos_str)}")
    annotation = annotation_dict
    annotation['target_position'] = counter_target_pos
    if caller_target_pos:
        annotation['paired_position'] = str(caller_target_pos)
        paired_position = str(caller_target_pos)
    else:
        paired_position = annotation.get('paired_position', None)
    anno_type = annotation['type']
    anno_desc = annotation.get('description', None)
    anno_id = f"{anno_type} | {anno_desc}"
    anno_count = 1
    anno_evidence = {entry_mnemo_name: annotation.get('evidence', None)}

    if anno_evidence[entry_mnemo_name]:
        if good_eco_codes:
            evidence_items = anno_evidence[entry_mnemo_name].split(',')
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
        }
        if caller_target_pos:
            anno_total['paired_position'] = paired_position
        if anno_type == 'BINDING':
            keys_to_exclude = {'type', 'description', 'count', 'evidence', 'paired_position', 'aminoacid'}
            if not caller_target_pos:
                keys_to_exclude.add('paired_position')
            for key, value in annotation.items():
                if key not in keys_to_exclude:
                    anno_total[key] = value

    visited_pos.add(counter_uniprot_pos_str)
    return {
        "annotation": annotation,
        "anno_type": anno_type,
        "anno_id": anno_id,
        "anno_total": anno_total,
        "paired_position": paired_position
        }

#@measure_time_and_memory
#@profile
def process_annotation(logger: logging.Logger, good_eco_codes: list, entry_mnemo_name: str, target_name: str, target_hit_start: int, target_hit_end: int, annotation_dict: Dict[str, Any], entry_annotations: Dict[str, Any], transfer_dict: dict, target_sequence: str, offset_start: int, offset_end: int, annot_sequence: str, counter_target_pos: int, counter_uniprot_pos_str: str, visited_pos: Set[str]) -> None:
    """
    Processes annotations to add them to the transfer dictionary by calling make_anno_total dict
    with anno_total containing at least type, description, count, and evidence, plus associated data.
    Additionally, it adds paired position if the other member also hits (successfull map_and_filter_annot_pos for it),
    and adds additional keys for BINDING annotations, if present.
    """
    result_dict = make_anno_total_dict(logger, good_eco_codes, entry_mnemo_name, annotation_dict, entry_annotations, counter_target_pos, counter_uniprot_pos_str, visited_pos)
    annotation = result_dict['annotation']
    anno_type = result_dict['anno_type']
    anno_id = result_dict['anno_id']
    anno_total = result_dict['anno_total']

    paired_position = result_dict['paired_position']
    if anno_total:
        entry_alphanum_name = annotation.get('entry', None)
        if anno_type in ['DISULFID', 'CROSSLNK', 'SITE', 'BINDING'] and paired_position is not None:
            target_paired_annotations = entry_annotations.get(paired_position, [])
            target_paired_annotation_dict = next((annotation for annotation in target_paired_annotations if annotation['type'] == anno_type), None)
            if target_paired_annotation_dict:
                paired_position_valid = map_and_filter_annot_pos(logger, good_eco_codes, target_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_mnemo_name, entry_annotations, transfer_dict, visited_pos, target_paired_uniprot_pos=paired_position, caller_target_pos=counter_target_pos, annotation_dict=target_paired_annotation_dict)
            else:
                paired_position_valid = False
            if paired_position_valid:
                logger.info(f"---> DEBUG --- PROCESS_ANNOT --- Paired Position Valid for target {target_name} and annotated {entry_mnemo_name}")
                paired_target_position = str(target_paired_annotation_dict['target_position'])
                annotation['paired_position'] = paired_target_position
                anno_total['paired_position'] = paired_target_position
                try:
                    add_to_transfer_dict(logger, transfer_dict, target_name, counter_target_pos, anno_id, anno_total, entry_mnemo_name, entry_alphanum_name)
                except Exception as e:
                    logger.error(f"---> DEBUG --- PROCESS_ANNOT --- Error in add_to_transfer_dict: {e}")
                    traceback.print_exc()
                    raise
                visited_pos.add(counter_uniprot_pos_str)
            else:
                remove_failed_annotations(entry_annotations, counter_uniprot_pos_str, paired_position, anno_type)
                visited_pos.add(counter_uniprot_pos_str)
        else:
            add_to_transfer_dict(logger, transfer_dict, target_name, counter_target_pos, anno_id, anno_total, entry_mnemo_name, entry_alphanum_name)
            visited_pos.add(counter_uniprot_pos_str)
    else:
        logger.info(f"---> DEBUG --- PROCESS_ANNOT --- NO ANNOTATION FOR SINGLE --- Anno Total was None for target {target_name} and annotated {entry_mnemo_name} at target {counter_target_pos} and annotated {counter_uniprot_pos_str}")
        visited_pos.add(counter_uniprot_pos_str)

#@measure_time_and_memory
#@profile
def validate_paired_positions(logger: logging.Logger, good_eco_codes: list, target_sequence: str, target_name: str, target_hit_start: int, target_hit_end: int, offset_start: int, offset_end: int, annot_sequence: str, entry_mnemo_name: str, annotation_dict: Dict[str, Any], entry_annotations: Dict[str, Any], transfer_dict: dict, visited_pos: Set[str], counter_target_pos: Optional[int], counter_uniprot_pos: Optional[int], target_paired_uniprot_pos: str, caller_target_pos: int) -> bool:
    """
    Validates the 2nd member of a paired position by checking
    if the target paired uniprot position from caller is a valid match.
    Success means adding the annotation to the transfer dictionary.
    Failure means returning paired_position_valid as False and
    asking to remove both positions from entry_annotations
    in the process_annotation parent function.
    """

    paired_position_valid = False
    int_target_paired_uniprot_pos = int(target_paired_uniprot_pos)
    last_target_pos = False

    logger.info(f"---> DEBUG --- VAL_PAIRED --- Running for target {target_name}/{target_hit_start}-{target_hit_end} and annotated {entry_mnemo_name}/{offset_start}-{offset_end}")
    while last_target_pos is False:
        for index, char in enumerate(annot_sequence):
            if char.isalpha():
                counter_uniprot_pos = int(offset_start) if counter_uniprot_pos is None else counter_uniprot_pos + 1

            if target_sequence[index].isalpha():
                counter_target_pos = int(target_hit_start) if counter_target_pos is None else counter_target_pos + 1

            if counter_target_pos is None and counter_uniprot_pos is None:
                continue

            if counter_uniprot_pos == int_target_paired_uniprot_pos and target_sequence[index] == annot_sequence[index]:
                # Calculate the start and end indices for the window
                logger.info("\n --- DEBUG --- VAL_PAIRED --- Helpful Info for Paired Annotation Match \n")
                start_index = max(0, index - 3)
                end_index = min(len(annot_sequence), index + 4)  # end_index is exclusive

                # Extract the window of amino acids
                annot_window = annot_sequence[start_index:end_index]
                target_window = target_sequence[start_index:end_index]

                logger.info(f"---> DEBUG --- VAL_PAIRED --- Counter Tar Pos {str(counter_target_pos)} and amino acid {target_sequence[index]} + Counter Uni Pos (annotated) {str(counter_uniprot_pos)} and amino acid {annot_sequence[index]}")
                logger.info(f"---> DEBUG --- VAL_PAIRED --- Annot Window: {annot_window} + Target Window: {target_window}")

                result_dict = make_anno_total_dict(logger, good_eco_codes, entry_mnemo_name, annotation_dict, entry_annotations, counter_target_pos, counter_uniprot_pos_str=target_paired_uniprot_pos, visited_pos=visited_pos, caller_target_pos=caller_target_pos)
                anno_id = result_dict['anno_id']
                anno_total = result_dict['anno_total']
                annotation = result_dict['annotation']
                if anno_total is None:
                    counter_uniprot_pos_str = str(counter_uniprot_pos)
                    logger.info(f"---> DEBUG --- VAL_PAIRED --- NO ANNOTATION FOR PAIR --- Anno Total was None - FAILED - for target {target_name} and annotated {entry_mnemo_name} at target {counter_target_pos} and annotated {counter_uniprot_pos_str}")
                    return paired_position_valid
                entry_alphanum_name = annotation.get('entry', None)
                add_to_transfer_dict(logger, transfer_dict, target_name, counter_target_pos, anno_id, anno_total, entry_mnemo_name, entry_alphanum_name)
                visited_pos.add(target_paired_uniprot_pos)
                paired_position_valid = True
                return paired_position_valid
            else:
                visited_pos.add(target_paired_uniprot_pos)

            if counter_target_pos == target_hit_end or counter_uniprot_pos == offset_end:
                last_target_pos = True
                break
    return paired_position_valid

#@measure_time_and_memory
#@profile
def validate_annotations(logger: logging.Logger, good_eco_codes: list, target_sequence: str, target_name: str, target_hit_start: int, target_hit_end: int, offset_start: int, offset_end: int, annot_sequence: str, entry_mnemo_name: str, entry_annotations: Dict[str, Any], transfer_dict: dict, visited_pos: Set[str], counter_target_pos: Optional[int], counter_uniprot_pos: Optional[int]) -> None:
    """
    Iterate on both annot and target sequence to find annotated positions that present the exact same amino acid in both sequences.
    These are sent to processing in process_annotation. Regardless, keep score of where we've been by adding to visited_pos set.
    """

    annotated_uni_pos_list = [int(key) for key in entry_annotations.keys()]
    if not any(pos in range(offset_start, offset_end + 1) for pos in annotated_uni_pos_list):
        return
    last_target_pos = False

    while last_target_pos is False:
        for index, char in enumerate(annot_sequence):
            if char.isalpha():
                counter_uniprot_pos = int(offset_start) if counter_uniprot_pos is None else counter_uniprot_pos + 1
            if target_sequence[index].isalpha():
                counter_target_pos = int(target_hit_start) if counter_target_pos is None else counter_target_pos + 1
            counter_uniprot_pos_str = str(counter_uniprot_pos)

            if counter_target_pos is None and counter_uniprot_pos is None:
                continue

            if counter_uniprot_pos_str not in visited_pos:
                visited_pos.add(counter_uniprot_pos_str)

                if counter_uniprot_pos_str in entry_annotations and target_sequence[index] == annot_sequence[index]:
                    counter_target_pos_str = str(counter_target_pos)

                    # DEBUGGING INFO
                    logger.info(f"\n --- DEBUG --- VAL_ANNOTS --- Helpful Info for Single/Caller Annotation Match \n NAMES: target {target_name} and annot {entry_mnemo_name} \n POSITIONS: target {counter_target_pos_str} and annot {counter_uniprot_pos_str} \n")
                    start_index = max(0, index - 3)
                    end_index = min(len(annot_sequence), index + 4)
                    annot_window = annot_sequence[start_index:end_index]
                    target_window = target_sequence[start_index:end_index]
                    logger.info(f"---> DEBUG --- VAL_ANNOTS --- Counter Uni Pos {str(counter_uniprot_pos)} and amino acid {annot_sequence[index]} + Counter Tar Pos {str(counter_target_pos)} and amino acid {target_sequence[index]}")
                    logger.info(f"---> DEBUG --- VAL_ANNOTS --- Annot Window: {annot_window} + Target Window: {target_window}")

                    for annotation_dict in entry_annotations[counter_uniprot_pos_str]:
                        try:
                            process_annotation(logger, good_eco_codes, entry_mnemo_name, target_name, target_hit_start, target_hit_end, annotation_dict, entry_annotations, transfer_dict, target_sequence, offset_start, offset_end, annot_sequence, counter_target_pos, counter_uniprot_pos_str, visited_pos)
                        except Exception as e:
                            logger.info(f"---> DEBUG --- VAL_ANNOTS --- Error in validate_annotations: {e}")
                            traceback.print_exc()
                            raise

            if counter_target_pos == target_hit_end or counter_uniprot_pos == offset_end:
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
    annotations_filepath = get_annotation_filepath(resource_dir, pfam_id)
    hmmalign_lines, annotations = read_files(dom_align, annotations_filepath)

    try:
        transfer_dict = find_and_map_annots(logger, hmmalign_lines, annotations, good_eco_codes)
        if not transfer_dict:
            logger.info("---> MAIN --- Transfer Dict was EMPTY")
        else:
            logger.info("---> MAIN --- Transfer Dict FILLED")
    except (KeyError, IndexError, AttributeError) as e:
        error_info = traceback.format_exc()
        logger.error("---> MAIN ---  Error transferring annotations for Pfam ID %s: %s\n%s", pfam_id, e, error_info)
        raise
    write_report(logger, transfer_dict, pfam_id, output_dir)

if __name__ == '__main__':
    outer_args = parse_arguments()
    outer_logger = configure_logging(outer_args.log)
    main(outer_logger)
