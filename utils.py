"""
utils.py

Copyright 2024 Eduardo Horta Santos <GitHub: Eduardo-HortaS>

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

This script contains all utility functions used by KRA.

"""

import argparse
import logging
import os
import pyhmmer.easel
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
from typing import Iterator
from datetime import datetime
from Bio import SeqIO
from collections import defaultdict
from typing import TypeVar, List, Callable, Optional, Literal, Any

Scope = Literal["main", "domain", "sequence", "seq_batch"]

def get_logger(log_path: str, scope: Scope = 'main', identifier: str = None) -> tuple[logging.Logger, str]:
    """
    Creates and configures a logger that writes to a timestamped log file.

    Parameters:
    log_path (str): The base path for the log file.
    scope (str): The scope of the logger, which can be 'main', 'domain', or 'sequence'. Defaults to 'main'.
    identifier (str, optional): An optional identifier to include in the log file name. Defaults to None.

    Returns:
    tuple[logging.Logger, str]: A tuple containing the configured logger and the final log file path.
    """
    # 1. Separate directory and filename
    log_dir = os.path.dirname(log_path)
    base_name = os.path.basename(log_path)
    base_stem, base_ext = os.path.splitext(base_name)

    # 2. Check for existing timestamp at the END of the filename
    timestamp = None
    if "_" in base_stem:
        # Split into parts and check only the LAST TWO segments
        parts = base_stem.split("_")
        if len(parts) >= 2:
            candidate_date = parts[-2]
            candidate_time = parts[-1]
            if (len(candidate_date) == 6 and candidate_date.isdigit() and
                len(candidate_time) == 4 and candidate_time.isdigit()):
                timestamp = f"{candidate_date}_{candidate_time}"
                base_stem = "_".join(parts[:-2])  # Rebuild base without timestamp

    # 3. Generate fresh timestamp if none found
    if not timestamp:
        timestamp = datetime.now().strftime("%y%m%d_%H%M")

    # 4. Build final filename
    if scope in ("domain", "sequence", "seq_batch") and identifier:
        clean_id = str(identifier).replace("|", "-").replace(" ", "_")[:64]
        filename = f"{base_stem}_{timestamp}_{clean_id}{base_ext}"
    else:
        filename = f"{base_stem}_{timestamp}{base_ext}"

    final_log = os.path.join(log_dir, filename)

    # 5. Configure logger
    os.makedirs(log_dir, exist_ok=True)
    logger = logging.getLogger(f"{scope}_{identifier}" if identifier else "main")
    if not logger.handlers:
        handler = logging.FileHandler(final_log)
        handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)

    return logger, final_log

LogLevel = Literal["debug", "info", "warning", "error", "critical"]

def get_multi_logger(loggers: List[logging.Logger]) -> Callable[[LogLevel, str, Any], None]:
    """Helper to write same message to multiple loggers.

    Args:
        loggers: List of logger instances to write to

    Returns:
        Function that takes:
            - level: Must be one of "debug", "info", "warning", "error", "critical"
            - message: Format string
            - args: Values for format string

    Example:
        >>> log = get_multi_logger([domain_logger, main_logger])
        >>> log("info", "Processing %s", domain_id)  # IDE suggests valid levels
    """
    def log(level: LogLevel, message: str, *args: Any) -> None:
        for logger in loggers:
            getattr(logger, level)(message, *args)
    return log

def get_querynames(fasta: str) -> List[str]:
    """Parses a FASTA file and returns a list of sequence "query names"."""
    sequences = []
    seq = SeqIO.parse(fasta, "fasta")
    for record in seq:
        queryname = record.id.replace("|", "-")
        sequences.append(queryname)
    return sequences

def make_dirs_and_write_fasta(sequences: Iterator[SeqRecord], base_dir: str) -> None:
    """Creates a directory for each sequence and writes the sequence to a fasta file."""
    for record in sequences:
        used_queryname = record.id.replace("|", "-")
        os.makedirs(os.path.join(base_dir, used_queryname), exist_ok=True)
        SeqIO.write(record, os.path.join(base_dir, used_queryname, "sequence.fasta"), "fasta")

def translate_sequence(seq_record: SeqRecord, logger: logging.Logger) -> SeqRecord:
    """Translates a nucleotide sequence using pyHMMER's translation table."""
    # Create digital sequence using pyHMMER
    digital_seq = pyhmmer.easel.TextSequence(
        name=seq_record.id.encode(),
        sequence=str(seq_record.seq).encode()
    ).digitize(pyhmmer.easel.Alphabet.dna())

    # Translate using pyHMMER
    translated = digital_seq.translate()
    logger.debug("Translated %s to amino acids.", seq_record.id)
    logger.debug("Original sequence: %s", seq_record.seq)
    logger.debug("Translated sequence: %s", translated.sequence.decode())

    # Create new SeqRecord with translated sequence
    return SeqRecord(
        Seq(translated.sequence.decode()),
        id=seq_record.id,
        description=seq_record.description
    )

def seqrecord_yielder(fasta: str, is_nucleotide: bool = False, logger: logging.Logger = None) -> Iterator[SeqRecord]:
    """Yields SeqRecord objects from a FASTA file, translating if nucleotide."""
    for record in SeqIO.parse(fasta, "fasta"):
        if is_nucleotide:
            yield translate_sequence(record, logger)
        else:
            yield record

def convert_lists_to_original_types(data):
    """Reverses convert_sets_and_tuples_to_lists operation.

    Note: This function is intended for downstream scripts.
    Consider moving to a shared utilities module.

    Special handling for annotation_ranges structure where:
        - positions (list) -> set
        - ranges (list of lists) -> list of tuples

    Args:
        data: Input data structure containing lists

    Returns:
        Converted data structure with appropriate sets/tuples
    """
    if isinstance(data, dict):
        # Special handling for annotation_ranges structure
        if 'positions' in data and 'ranges' in data:
            return {
                'positions': set(data['positions']),
                'ranges': [tuple(r) for r in data['ranges']]
            }
        return {k: convert_lists_to_original_types(v) for k, v in data.items()}
    elif isinstance(data, list):
        return data  # Keep as list by default
    else:
        return data

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
        if 'positions' in data and 'ranges' in data:
            return {
                'positions': sorted(list(data['positions'])),
                'ranges': [list(r) for r in data['ranges']]
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

def convert_defaultdict_to_dict(obj):
    """Recursively convert defaultdicts to dicts."""
    if isinstance(obj, defaultdict):
        obj = dict(obj)
    if isinstance(obj, dict):
        return {k: convert_defaultdict_to_dict(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_defaultdict_to_dict(item) for item in obj]
    return obj
