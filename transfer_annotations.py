
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
Intense usage begins by calling find_and_map_annotations with 2 arguments:
a list of hmmalign result lines and the loaded dict of annotations.

"""

import os
import json
import argparse
import logging
import traceback
# from modules.decorators import measure_time_and_memory
# from memory_profiler import profile


def parse_arguments():
    """Parse command-line arguments for running hmmsearch
    using a sequence database against target HMMs,
    aiming to generate a hmmsearch_per_domain_output.json.
    Also includes an output directory and a log file path.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description=
    'Generates a temporary multifasta for running hmmalign using a hits per domain JSON.')
    parser.add_argument("-iA", "--dom-align", help="Path to domain's hmmalign alignment", required=True, type=str)
    parser.add_argument("-r", "--resource-dir", help="Resource dir path", required=True, type=str)
    parser.add_argument("-o", "--output-dir", help="Output dir path", required=True, type=str)
    parser.add_argument("-l", "--log", help="Log path", \
        required=False, type=str, default="logs/transfer_annotations.log")
    return parser.parse_args()

def setup_logging(log_path):
    """Set up logging for the script."""
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(filename=log_path, level=logging.DEBUG, filemode='w', \
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def get_pfam_id_from_hmmalign_result(hmmalign_result):
    """
    Extracts the PFAM ID from the hmmalign result filename.
    """
    return os.path.splitext(os.path.basename(hmmalign_result))[0].split('_')[-2]

def get_annotation_filepath(resource_dir, pfam_id):
    """
    Constructs and returns the filepaths for annotations,
    seed alignment, and report based on the PFAM ID.
    """
    annotations_filepath = os.path.join(resource_dir, pfam_id, "annotations.json")
    return annotations_filepath

def read_files(hmmalign_result, annotations_filepath):
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
def find_and_map_annots(hmmalign_lines, annotations):
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

    for entry_name, entry_annotations in annotations.items():
        for line in hmmalign_lines:
            if line.startswith("#") or 'target/' in line:
                continue
            splitted = line.split('/')
            entry_name_in_line = splitted[0]
            if entry_name == entry_name_in_line:
                offset_start = int(splitted[1].split('-')[0])
                offset_end = int(splitted[1].split('-')[1].split()[0])
                annot_sequence = splitted[1].split()[1]
                for target_name, target_per_interval in target_info.items():
                    for target_hit_interval, target_info_list in target_per_interval.items():
                        if target_name not in visited_pos:
                            visited_pos[target_name] = {}
                        if entry_name not in visited_pos[target_name]:
                            visited_pos[target_name][entry_name] = {}
                        if target_hit_interval not in visited_pos[target_name][entry_name]:
                            visited_pos[target_name][entry_name][target_hit_interval] = set()
                        for target_hit_start, target_hit_end, target_hit_sequence in target_info_list:
                            try:
                                map_and_filter_annot_pos(target_hit_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_name, entry_annotations, transfer_dict, visited_pos[target_name][entry_name][target_hit_interval])
                                print(f"Mapped and filtered annotation for {entry_name} and {target_name}")
                            except Exception as e:
                                error_info = traceback.format_exc()
                                print(f"********** --- Error in map_and_filter_annot_pos: {e} \n Error Info: {error_info}")
    return transfer_dict

# NOTE: Bugfix pending, 30 07
def write_report(transfer_dict, pfam_id, output_dir):
    """
    Writes the transfer dictionary to a report file in JSON format.
    """
    for target_name in transfer_dict.keys():
        print(transfer_dict[target_name])
        if transfer_dict[target_name]:
            json_data = json.dumps(transfer_dict[target_name], indent=4)
            print(f"We got data for {target_name}!")
        else:
            json_data = json.dumps({"Annotations": "None"}, indent=4)
            print(f"No data for {target_name}!")
        target_name = target_name.replace("|", "-")
        pfam_id_report_filepath = os.path.join(output_dir, target_name, pfam_id + "_report.json")
        with open(pfam_id_report_filepath, 'w', encoding="utf-8") as report_file:
            report_file.write(json_data)
        print(f"**************************** Transfer Report for sequence-domain pair, {target_name}-{pfam_id}: {pfam_id_report_filepath}")

#@measure_time_and_memory
##@profile
def map_and_filter_annot_pos(target_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_name, entry_annotations, transfer_dict, visited_pos, target_paired_uniprot_pos=None, caller_target_pos=None):
    """
    Asks to map positions between annotated and target sequences.
    If target paired uniprot position and caller target position are provided,
    this is an associated call for a paired position sent to validate_paired_position.
    Otherwise, it is a call for a single position sent to validate_annotations.
    Makes counter_uniprot_pos and counter_target_pos to track the current position in both sequences.
    """
    counter_uniprot_pos = None
    counter_target_pos = None

    if target_paired_uniprot_pos and caller_target_pos:
        return validate_paired_position(target_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_name, entry_annotations, transfer_dict, visited_pos, counter_target_pos, counter_uniprot_pos, target_paired_uniprot_pos, caller_target_pos)
    return validate_annotations(target_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_name, entry_annotations, transfer_dict, visited_pos, counter_target_pos, counter_uniprot_pos)

#@measure_time_and_memory
#@profile
# NOTE: Bugfix pending, 30 07
def add_to_transfer_dict(transfer_dict, target_name, counter_target_pos, anno_id, anno_total, entry_name):
    """
    Adds anno_total data to the transfer dictionary.
    The dictionary is structured by target position, annotation ID (type + description),
    and various possible additional annotation details, mainly for BINDING annotations.
    """
    additional_keys = {k: v for k, v in anno_total.items() if k not in ['type', 'description', 'count', 'evidence', 'paired_position', 'aminoacid', 'target_position']}

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
                "rep_entry_name": entry_name, "count": 1}
        else:
            transfer_dict[target_name][counter_target_pos][anno_id]['evidence'][evidence_value]["count"] += 1

    if paired_position_value:
        if 'paired_position' not in transfer_dict[target_name][counter_target_pos][anno_id]:
            transfer_dict[target_name][counter_target_pos][anno_id].setdefault('paired_position', {})
        if paired_position_value not in transfer_dict[target_name][counter_target_pos][anno_id]['paired_position']:
            transfer_dict[target_name][counter_target_pos][anno_id]['paired_position'][paired_position_value] = {
                "rep_entry_name": entry_name, "count": 1}
        else:
            transfer_dict[target_name][counter_target_pos][anno_id]['paired_position'][paired_position_value]["count"] += 1

    if additional_keys:
        if 'additional_keys' not in transfer_dict[target_name][counter_target_pos][anno_id]:
            transfer_dict[target_name][counter_target_pos][anno_id].setdefault('additional_keys', {})
        for key, value in additional_keys.items():
            if anno_total['type'] == 'BINDING':
                if key not in transfer_dict[target_name][counter_target_pos][anno_id]['additional_keys']:
                    transfer_dict[target_name][counter_target_pos][anno_id]['additional_keys'][key] = {
                        "rep_entry_name": entry_name, "count": 1, "value": value}
                else:
                    transfer_dict[target_name][counter_target_pos][anno_id]['additional_keys'][key]["count"] += 1
    print(transfer_dict)

#@measure_time_and_memory
#@profile
def remove_failed_annotations(entry_annotations, counter_uniprot_pos_str, paired_position):
    """
    Removes failed annotations from entry annotations.
    """
    del entry_annotations[counter_uniprot_pos_str][0]
    if paired_position in entry_annotations and entry_annotations[paired_position]:
        del entry_annotations[paired_position][0]

#@measure_time_and_memory
#@profile
def make_anno_total_dict(entry_name, target_name, entry_annotations, counter_target_pos, counter_uniprot_pos_str, visited_pos, caller_target_pos=None):
    """
    For any given annotation in appropriate format (list with 1 dict inside), paired or not, produces an anno_total dict and
    associated data to return as a dict to be extracted from in accordance to the calling function needs.
    Used by process_annotation and validate_paired_position.
    """
    good_eco_codes = ["ECO:0000269", "ECO:0000303", "ECO:0000305", "ECO:0000312", "ECO:0007744"]
    print(f"Entry Annotations in make_anno_total_dict --- {entry_annotations}")
    if isinstance(entry_annotations[counter_uniprot_pos_str], list) and entry_annotations[counter_uniprot_pos_str]:
        anno_total = None
        annotation = entry_annotations[counter_uniprot_pos_str][0]
        annotation['target_position'] = counter_target_pos
        if caller_target_pos:
            annotation['paired_position'] = str(caller_target_pos)
            paired_position = str(caller_target_pos)
            anno_paired_position = {entry_name: paired_position}
            paired_position = None
        else:
            paired_position = annotation.get('paired_position', None)
        anno_type = annotation['type']
        anno_desc = annotation.get('description', None)
        anno_id = f"{anno_type} | {anno_desc}"
        anno_count = 1
        anno_evidence = {entry_name: annotation.get('evidence', None)}

        if anno_evidence[entry_name]:
            evidence_items = anno_evidence[entry_name].split(',')
            # print(f"Evidence Items Before Filtering: {evidence_items}")
            valid_evidence_items = []
            for item in evidence_items:
                code, _, _ = item.strip().partition('|')  # Split by the first occurrence of '|'
                if code.strip() in good_eco_codes:
                    valid_evidence_items.append(item.strip())  # Preserve the whole item if the code is valid
            # print(f"Valid Evidence Items: {valid_evidence_items}")
            anno_evidence[entry_name] = ', '.join(valid_evidence_items)
        if anno_evidence[entry_name]:
            anno_total = {
                "type": anno_type,
                "description": anno_desc,
                "count": anno_count,
                "evidence": anno_evidence[entry_name],
            }
            # print(f"Look at anno_total's evidence: {anno_total['evidence']}")
            if caller_target_pos:
                anno_total['paired_position'] = anno_paired_position

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
    else:
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
def process_annotation(entry_name, target_name, target_hit_start, target_hit_end, entry_annotations, transfer_dict, target_sequence, offset_start, offset_end, annot_sequence, counter_target_pos, counter_uniprot_pos_str, visited_pos):
    """
    Processes annotations to add them to the transfer dictionary by calling make_anno_total dict
    with anno_total containing at least type, description, count, and evidence, plus associated data.
    Additionally, it adds paired position if the other member also hits (successfull map_and_filter_annot_pos for it),
    and adds additional keys for BINDING annotations, if present.
    """
    result_dict = make_anno_total_dict(entry_name, target_name, entry_annotations, counter_target_pos, counter_uniprot_pos_str, visited_pos)
    annotation = result_dict['annotation']
    anno_type = result_dict['anno_type']
    anno_id = result_dict['anno_id']
    anno_total = result_dict['anno_total']
    paired_position = result_dict['paired_position']
    if anno_total:
        if anno_type in ['DISULFID', 'CROSSLNK', 'SITE', 'BINDING'] and paired_position is not None:
            paired_position_valid = map_and_filter_annot_pos(target_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_name, entry_annotations, transfer_dict, visited_pos, target_paired_uniprot_pos=paired_position, caller_target_pos=counter_target_pos)
            if paired_position_valid:
                paired_target_position = str(entry_annotations[paired_position][0]['target_position'])
                annotation['paired_position'] = paired_target_position
                anno_total['paired_position'] = paired_target_position
                add_to_transfer_dict(transfer_dict, target_name, counter_target_pos, anno_id, anno_total, entry_name)
                visited_pos.add(counter_uniprot_pos_str)
            else:
                remove_failed_annotations(entry_annotations, counter_uniprot_pos_str, paired_position)
                visited_pos.add(counter_uniprot_pos_str)
        else:
            add_to_transfer_dict(transfer_dict, target_name, counter_target_pos, anno_id, anno_total, entry_name)
            visited_pos.add(counter_uniprot_pos_str)
    else:
        visited_pos.add(counter_uniprot_pos_str)

#@measure_time_and_memory
#@profile
def validate_paired_position(target_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_name, entry_annotations, transfer_dict, visited_pos, counter_target_pos, counter_uniprot_pos, target_paired_uniprot_pos, caller_target_pos):
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

    print(f" --- VALIDATE PAIRED --- Running for {target_name}/{target_hit_start}-{target_hit_end} and {entry_name}/{offset_start}-{offset_end}")
    while last_target_pos is False:
        for index, char in enumerate(annot_sequence):
            if char.isalpha():
                counter_uniprot_pos = int(offset_start) if counter_uniprot_pos is None else counter_uniprot_pos + 1

            if target_sequence[index].isalpha():
                counter_target_pos = int(target_hit_start) if counter_target_pos is None else counter_target_pos + 1

            if counter_target_pos is None and counter_uniprot_pos is None:
                continue

            if counter_uniprot_pos == int_target_paired_uniprot_pos and target_sequence[index] == annot_sequence[index]:
                result_dict = make_anno_total_dict(entry_name, target_name, entry_annotations, counter_target_pos, target_paired_uniprot_pos, visited_pos, caller_target_pos)
                anno_id = result_dict['anno_id']
                anno_total = result_dict['anno_total']
                add_to_transfer_dict(transfer_dict, target_name, counter_target_pos, anno_id, anno_total, entry_name)
                visited_pos.add(target_paired_uniprot_pos)
                paired_position_valid = True
                return paired_position_valid
            else:
                visited_pos.add(target_paired_uniprot_pos)
            if counter_target_pos == target_hit_end or counter_uniprot_pos == offset_end:
                print(f"Ran out!")
                print(f"Visited pos brother {visited_pos}")
                last_target_pos = True
                break


    return paired_position_valid

#@measure_time_and_memory
#@profile
def validate_annotations(target_sequence, target_name, target_hit_start, target_hit_end, offset_start, offset_end, annot_sequence, entry_name, entry_annotations, transfer_dict, visited_pos, counter_target_pos, counter_uniprot_pos):
    """
    Iterate on both annot and target sequence to find annotated positions that present the exact same amino acid in both sequences.
    These are sent to processing in process_annotation. Regardless, keep score of where we've been by adding to visited_pos set.
    """

    annotated_uni_pos_list = [int(key) for key in entry_annotations.keys()]
    if not any(pos in range(offset_start, offset_end + 1) for pos in annotated_uni_pos_list):
        return

    print(f" --- VAL_ANNOTS --- Running for {target_name}/{target_hit_start}-{target_hit_end} and {entry_name}/{offset_start}-{offset_end}")
    print(f"Entry Annotations --- {entry_annotations}")

    last_target_pos = False
    # print(f"Visited Pos in validate annotations --- {visited_pos}")

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
                print(f"Counter Uni Pos {counter_uniprot_pos_str} = * {annot_sequence[index]} * \nCounter Tar Pos {str(counter_target_pos)} = ! {target_sequence[index]} !")

                if counter_uniprot_pos_str in entry_annotations and target_sequence[index] == annot_sequence[index]:
                    print(f" --- VAL_ANNOTS --- Went inside!!! For {target_name} and {entry_name}!")
                    counter_target_pos_str = str(counter_target_pos)
                    print(f" --- VAL_ANNOTS --- We found a pos match for {target_name} and {entry_name} at target {counter_target_pos_str} and annot {counter_uniprot_pos_str}!")
                    process_annotation(entry_name, target_name, target_hit_start, target_hit_end, entry_annotations, transfer_dict, target_sequence, offset_start, offset_end, annot_sequence, counter_target_pos_str, counter_uniprot_pos_str, visited_pos)

            if counter_target_pos == target_hit_end or counter_uniprot_pos == offset_end:
                print(f"Ran out!")
                print(f"Visited pos brother {visited_pos}")
                last_target_pos = True
                break

def main():
    """Main function, initializes this script"""
    args = parse_arguments()
    dom_align = args.dom_align
    resource_dir = args.resource_dir
    output_dir = args.output_dir
    setup_logging(args.log)

    pfam_id = get_pfam_id_from_hmmalign_result(dom_align)
    annotations_filepath = get_annotation_filepath(resource_dir, pfam_id)
    hmmalign_lines, annotations = read_files(dom_align, annotations_filepath)
    try:
        transfer_dict = find_and_map_annots(hmmalign_lines, annotations)
        print("!!Transferidooooo!!")
        print(transfer_dict)
    except (KeyError, IndexError) as e:
        error_info = traceback.format_exc()
        print(f"Error transferring annotations for Pfam ID {pfam_id}: {e}\n{error_info}")
        transfer_dict = {}
    write_report(transfer_dict, pfam_id, output_dir)

if __name__ == '__main__':
    main()
