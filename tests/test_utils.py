# import logging
# import json
import sys
import os
# import copy
# import pandas as pd
# from io import StringIO
# from importlib import reload
# from argparse import Namespace
# from unittest.mock import patch, ANY, call, mock_open, MagicMock
# from tempfile import TemporaryDirectory

# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from utils import (
    get_querynames,
    make_dirs_and_write_fasta,
    translate_sequence,
    seqrecord_yielder,
    parallelize,
    convert_lists_to_original_types,
    convert_sets_and_tuples_to_lists,
    parse_arguments,
    setup_logging
)

import pytest

### Fixtures

@pytest.fixture
def transfer_dict_populated_disulfid_Q9NU22():
    """Populated transfer_dict with DISULFID annotations for a single pair, 333 and 373 in Q9NU22.
    Note that sets and tuples will be converted to lists before JSON serialization."""
    return {
    "DOMAIN": {
        "sequence_id": {
            "sp|Q9NU22|MDN1_HUMAN": {
                "hit_intervals": {
                    "325-451": {
                        "sequence": "VLLEGPIGCGKTSLVEYLAAVTGRTKPPQLLKVQLGDQTDSKMLLGMYCCTDVPGEFVWQPGTLTQAATMGHWILLEDIDYAPLDVVSVLIPLLENGELLIPGRGDCLKVAPGFQFFATRRLLSCGG",
                        "length": 127,
                        "hit_start": 325,
                        "hit_end": 451,
                        "annotations": {
                            "positions": {
                                "333": {
                                    "DISULFID | Intrachain (with C-246); in linked form": {
                                        "essentials": {
                                            "type": "DISULFID",
                                            "description": "Intrachain (with C-246); in linked form",
                                            "count": 1,
                                            "annot_amino_acid": "C",
                                            "target_amino_acid": "C"
                                        },
                                        "evidence": {
                                            "ECO:0000269|PubMed:12345678": {
                                                "rep_primary_accession": "P15005",
                                                "rep_mnemo_name": "MCRB_ECOLI",
                                                "count": 1
                                            }
                                        },
                                        "paired_position": {
                                            "373": {
                                                "rep_primary_accession": "P15005",
                                                "rep_mnemo_name": "MCRB_ECOLI",
                                                "count": 1
                                            }
                                        },
                                        "hit" : True,
                                        "additional_keys": {
                                            "annot_position": {
                                                "205": {
                                                    "rep_primary_accession": "P15005",
                                                    "rep_mnemo_name": "MCRB_ECOLI",
                                                    "count": 1
                                                }
                                            }
                                        }
                                    }
                                },
                                "373": {
                                    "DISULFID | Intrachain (with C-205); in linked form": {
                                        "essentials": {
                                            "type": "DISULFID",
                                            "description": "Intrachain (with C-205); in linked form",
                                            "count": 1,
                                            "annot_amino_acid": "C",
                                            "target_amino_acid": "C"
                                        },
                                        "evidence": {
                                            "ECO:0000269|PubMed:12345678": {
                                                "rep_primary_accession": "P15005",
                                                "rep_mnemo_name": "MCRB_ECOLI",
                                                "count": 1
                                            }
                                        },
                                        "paired_position": {
                                            "333": {
                                                "rep_primary_accession": "P15005",
                                                "rep_mnemo_name": "MCRB_ECOLI",
                                                "count": 1
                                            }
                                        },
                                        "hit" : True,
                                        "additional_keys": {
                                            "annot_position": {
                                                "246": {
                                                    "rep_primary_accession": "P15005",
                                                    "rep_mnemo_name": "MCRB_ECOLI",
                                                    "count": 1
                                                }
                                            }
                                        }
                                    }
                                }
                            },
                            "indices": {
                                "matches": {"333", "373"}, # set literal
                                "misses": set()
                            }
                        },
                        "conservations": {
                            "positions": {},
                            "indices": {
                                "matches": set(),
                                "misses": set()
                            }
                        },
                        "position_conversion": {
                            "target_to_aln": {
                                "333": "18",
                                "373": "72",
                                },
                            "aln_to_target": {
                                # Used to be set literals, not needed since no target-ali mapping repetition within an interval
                                "18": "333",
                                "72": "373",
                                }
                            },
                            "annotation_ranges": {
                                "DISULFID | Intrachain (with C-246); in linked form": {
                                    "positions": {333},
                                    "ranges": [(333, 333)]
                                },
                                "DISULFID | Intrachain (with C-205); in linked form": {
                                    "positions": {373},
                                    "ranges": [(373, 373)]
                                }
                            }
                        }
                    }
                }
            }
        }
    }

@pytest.fixture
def transfer_dict_populated_disulfid_list_original_structure_Q9NU22(transfer_dict_populated_disulfid_Q9NU22):
    """Convert fixture to list version to test conversion back, using original structure.
    transfer_dict[DOMAIN][sequence_id][sequence_name][hit_intervals][interval]..."""
    return convert_sets_and_tuples_to_lists(transfer_dict_populated_disulfid_Q9NU22)


###T convert_lists_to_original_types

def test_convert_lists_to_original_types_annotation_ranges():
    """Test conversion of annotation_ranges structure back to sets and tuples"""
    test_data = {
        "positions": [333],
        "ranges": [[333, 333]]
    }

    result = convert_lists_to_original_types(test_data)

    assert isinstance(result["positions"], set)
    assert result["positions"] == {333}
    assert isinstance(result["ranges"][0], tuple)
    assert result["ranges"] == [(333, 333)]

def test_convert_lists_to_original_types_indices():
    """Test conversion of indices matches/misses back to sets"""
    test_data = {
        "indices": {
            "matches": ["333", "373"],
            "misses": []
        }
    }

    result = convert_lists_to_original_types(test_data)

    # Lists should stay as lists since they're not in annotation_ranges structure
    assert isinstance(result["indices"]["matches"], list)
    assert isinstance(result["indices"]["misses"], list)

def test_convert_lists_to_original_types_real_data(transfer_dict_populated_disulfid_list_original_structure_Q9NU22):
    """Test conversion of real transfer dict data back to original types"""
    result = convert_lists_to_original_types(transfer_dict_populated_disulfid_list_original_structure_Q9NU22)

    # Check annotation_ranges conversion
    anno_ranges = result["DOMAIN"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]["annotation_ranges"]
    for key in anno_ranges:
        assert isinstance(anno_ranges[key]["positions"], set)
        assert isinstance(anno_ranges[key]["ranges"][0], tuple)

    # Check indices remain as lists
    indices = result["DOMAIN"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]["annotations"]["indices"]
    assert isinstance(indices["matches"], list)
    assert isinstance(indices["misses"], list)

