import logging
import json
import sys
import os
from argparse import Namespace
from collections import defaultdict
from unittest.mock import patch, ANY, MagicMock
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from make_view_jsons import (
    parse_arguments,
    merge_nested_data,
    aggregate_range_positions,
    transform_to_ranges,
    process_sequence_report,
    write_range_views,
    main
)

import pytest

### Fixtures

@pytest.fixture
def logger():
    """Create a mock logger that also prints log messages to the console."""
    mock_logger = MagicMock()

    # Define side effect functions to print log messages
    def info_side_effect(msg, *args):
        formatted_msg = msg % args if args else msg
        print(f"INFO: {formatted_msg}")

    def error_side_effect(msg, *args):
        formatted_msg = msg % args if args else msg
        print(f"ERROR: {formatted_msg}")

    def warning_side_effect(msg, *args):
        formatted_msg = msg % args if args else msg
        print(f"WARNING: {formatted_msg}")

    def debug_side_effect(msg, *args):
        formatted_msg = msg % args if args else msg
        print(f"DEBUG: {formatted_msg}")

    # Assign side effects to the mock logger methods
    mock_logger.info.side_effect = info_side_effect
    mock_logger.error.side_effect = error_side_effect
    mock_logger.warning.side_effect = warning_side_effect
    mock_logger.debug.side_effect = debug_side_effect

    return mock_logger

@pytest.fixture
def multi_logger():
    """Create a mock multi_logger that matches the get_multi_logger behavior"""
    def multi_log(level: str, message: str, *args):
        # This simulates the behavior of the actual multi_logger
        # which formats the message with args and logs to multiple loggers
        return getattr(logging.getLogger(), level)(message, *args)

    mock = MagicMock(side_effect=multi_log)
    return mock

@pytest.fixture
def mock_aggregated_report_data():
    return {
        "sp|Q9NU22|MDN1_HUMAN": {
            "PF07728": {
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
                                            "annot_residue": "C",
                                            "target_residue": "C"
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
                                        "hit": True,
                                        "additional_keys": {
                                            "annot_position": {
                                                "205": {
                                                    "rep_primary_accession": "P15005",
                                                    "rep_mnemo_name": "MCRB_ECOLI",
                                                    "count": 1
                                                }
                                            }
                                        },
                                        "GO": {
                                            "MCRB_ECOLI": {
                                                "terms": {
                                                    "BP": {},
                                                    "MF": {
                                                        "GO:0005524": "ATP binding",
                                                        "GO:0016887": "ATP hydrolysis activity"
                                                    }
                                                },
                                                "status": "normal",
                                                "wang_sem_sim_bma_bp": 0.0,
                                                "wang_sem_sim_bma_mf": 0.541
                                            }
                                        }
                                    }
                                }
                            },
                            "indices": {
                                "matches": ["333"],
                                "misses": []
                            }
                        },
                        "conservations": {
                            "positions": {
                                # Mocked a continuous range to test
                                "364": {"conservation": 0.90, "cons_residue": "G", "target_residue": "G", "hit": True},
                                "365": {"conservation": 0.90, "cons_residue": "G", "target_residue": "G", "hit": True},
                                "366": {"conservation": 0.91, "cons_residue": "G", "target_residue": "T", "hit": False},
                                "367": {"conservation": 0.92, "cons_residue": "G", "target_residue": "G", "hit": True},
                                "368": {"conservation": 0.93, "cons_residue": "G", "target_residue": "G", "hit": True},
                                "369": {"conservation": 0.94, "cons_residue": "G", "target_residue": "T", "hit": True},
                                "370": {"conservation": 0.95, "cons_residue": "G", "target_residue": "G", "hit": False},
                                "371": {"conservation": 0.96, "cons_residue": "G", "target_residue": "G", "hit": True},
                                "372": {"conservation": 0.97, "cons_residue": "G", "target_residue": "G", "hit": True},
                                "373": {"conservation": 0.98, "cons_residue": "G", "target_residue": "G", "hit": True},
                                "374": {"conservation": 0.99, "cons_residue": "G", "target_residue": "G", "hit": True},
                            },
                            "indices": {
                                "matches": ["364", "365", "367", "368", "369", "371", "372", "373", "374"],
                                "misses": ["366", "370"]
                            }
                        },
                        "position_conversion": {
                            "target_to_aln": {
                                "333": "18",
                                "364": "63",
                                "365": "64",
                                "366": "65",
                                "367": "66",
                                "368": "67",
                                "369": "68",
                                "370": "69",
                                "371": "70",
                                "372": "71",
                                "373": "72",
                                "374": "73",
                                },
                            "aln_to_target": {
                                "18": "333",
                                "63": "364",
                                "64": "365",
                                "65": "366",
                                "66": "367",
                                "67": "368",
                                "68": "369",
                                "69": "370",
                                "70": "371",
                                "71": "372",
                                "72": "373",
                                "73": "374",
                                }
                        },
                        "annotation_ranges": {
                            "DISULFID | Intrachain (with C-246); in linked form": {
                                "positions": [333],
                                "ranges": [[333, 333]]
                            }
                        },
                        "conservation_ranges": {
                            "conserved_positions": {
                                "positions": [364],
                                "ranges": [[364, 374]]
                            }
                        }
                    }
                }
            }
        }
    }

@pytest.fixture
def mock_transformed_data():
    return {
        "sequence_id": "sp|Q9NU22|MDN1_HUMAN",
            "domain": {
                "PF07728": {
                    "325-451": {
                        "DISULFID | Intrachain (with C-246); in linked form": {
                            "333-333": {
                                "hit": False,
                                "essentials": {
                                "count": {
                                    "333_DISULFID | Intrachain (with C-246); in linked form": 1
                                },
                                "type": "DISULFID",
                                "description": "Intrachain (with C-246); in linked form",
                                "annot_residue": {
                                    "333": "C"
                                },
                                "target_residue": {
                                    "333": "C"
                                },
                                },
                                "evidence": {
                                "ECO:0000269|PubMed:12345678": {
                                    "count": {
                                    "333_DISULFID | Intrachain (with C-246); in linked form": 1
                                    },
                                    "rep_primary_accession": "P15005",
                                    "rep_mnemo_name": "MCRB_ECOLI"
                                }
                                },
                                "paired_position": {
                                "373": {
                                    "count": {
                                    "333_DISULFID | Intrachain (with C-246); in linked form": 1
                                    },
                                    "rep_primary_accession": "P15005",
                                    "rep_mnemo_name": "MCRB_ECOLI"
                                }
                                },
                                "additional_keys": {
                                "annot_position": {
                                    "205": {
                                    "count": {
                                        "333_DISULFID | Intrachain (with C-246); in linked form": 1
                                    },
                                    "rep_primary_accession": "P15005",
                                    "rep_mnemo_name": "MCRB_ECOLI"
                                    }
                                }
                                }
                            }
                            },
                            "annotations": {
                            "indices": {
                                "matches": [
                                "333"
                                ],
                                "misses": []
                            }
                            },
                            "conserved_positions": {
                            "364-374": {
                                "conservation": 0.9409090909090909,
                                "hit": False
                            }
                            },
                            "conservations": {
                            "indices": {
                                "matches": [
                                "364",
                                "365",
                                "367",
                                "368",
                                "369",
                                "371",
                                "372",
                                "373",
                                "374"
                                ],
                                "misses": [
                                "366",
                                "370"
                                ]
          }
        }
      }
    }
  }
}


@pytest.fixture
def mock_annotation_subset():
    """Subset of mock_aggregated_report_data focusing on annotations"""
    return {
        "annotations": {
            "positions": {
                "333": {
                    "DISULFID | Intrachain (with C-246); in linked form": {
                        "essentials": {
                            "type": "DISULFID",
                            "description": "Intrachain (with C-246); in linked form",
                            "count": 1,
                            "annot_residue": "C",
                            "target_residue": "C"
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
                        "hit": True,
                        "additional_keys": {
                            "annot_position": {
                                "205": {
                                    "rep_primary_accession": "P15005",
                                    "rep_mnemo_name": "MCRB_ECOLI",
                                    "count": 1
                                }
                            }
                        },
                        "GO": {
                            "MCRB_ECOLI": {
                                "terms": {
                                    "BP": {},
                                    "MF": {
                                        "GO:0005524": "ATP binding",
                                        "GO:0016887": "ATP hydrolysis activity"
                                    }
                                },
                                "status": "normal",
                                "wang_sem_sim_bma_bp": 0.0,
                                "wang_sem_sim_bma_mf": 0.541
                            }
                        }
                    }
                }
            }
        },
        "annotation_ranges": {
            "DISULFID | Intrachain (with C-246); in linked form": {
                "positions": [333],
                "ranges": [[333, 333]]
            }
        }
    }

@pytest.fixture
def mock_conservation_subset():
    """Subset of mock_aggregated_report_data focusing on conservations"""
    return {
        "conservations": {
            "positions": {
                "364": {"conservation": 0.90, "cons_residue": "G", "target_residue": "G", "hit": True},
                "365": {"conservation": 0.90, "cons_residue": "G", "target_residue": "G", "hit": True}
            },
            "indices": {
                "matches": ["364", "365"],
                "misses": []
            }
        },
        "conservation_ranges": {
            "conserved_positions": {
                "positions": [364, 365],
                "ranges": [[364, 365]]
            }
        }
    }

#### Mock variables
# sequence_dir_mock = "/home/user/results/sequence1/"
# log_filepath_mock = "/home/user/logs/make_views.log"
# sequence_id_mock = "sp|Q9NU22|MDN1_HUMAN"
# clean_sequence_id_mock = "sp-Q9NU22-MDN1_HUMAN"

###T parse_arguments - Unit tests

def test_parse_arguments_required():
    sequence_dir_mock = "/home/user/results/sequence1/"
    test_args = [
        "-sD", sequence_dir_mock,
        "-s", "sp|Q9NU22|MDN1_HUMAN"
    ]

    with pytest.raises(SystemExit):
        parse_arguments()

    with pytest.MonkeyPatch().context() as mp:
        mp.setattr("sys.argv", ["script_name.py"] + test_args)
        args = parse_arguments()

    expected = Namespace(
        sequence_dir="/home/user/results/sequence1/",
        sequence="sp|Q9NU22|MDN1_HUMAN",
        log="logs/make_views_jsons.log"
    )
    assert vars(args) == vars(expected)

def test_parse_arguments_optional():
    sequence_dir_mock = "/home/user/results/sequence1/"
    log_filepath_mock = "/home/user/logs/make_views.log"

    test_args = [
        "-sD", sequence_dir_mock,
        "-s", "sp|Q9NU22|MDN1_HUMAN",
        "-l", log_filepath_mock
    ]

    with pytest.MonkeyPatch().context() as mp:
        mp.setattr("sys.argv", ["script_name.py"] + test_args)
        args = parse_arguments()

    assert args.sequence_dir == sequence_dir_mock
    assert args.log == log_filepath_mock

###T merge_nested_data - Integration Tests (+ track_position_data())

def test_merge_nested_data(mock_aggregated_report_data):
    """Test nested data merging with position tracking"""
    source = mock_aggregated_report_data["sp|Q9NU22|MDN1_HUMAN"]["PF07728"]["hit_intervals"]["325-451"]["annotations"]["positions"]["333"]["DISULFID | Intrachain (with C-246); in linked form"]
    target = defaultdict(lambda: defaultdict(dict))

    merge_nested_data(source, target, "333")
    print(json.dumps(target, indent=2))
    assert target["essentials"]["type"] == "DISULFID"
    assert "count" in target["essentials"]
    assert isinstance(target["essentials"]["count"], dict)
    assert target["essentials"]["count"]["333"] == 1
    assert target["GO"]["333"]["MCRB_ECOLI"]["terms"]["MF"]["GO:0016887"] == "ATP hydrolysis activity"
    assert target["GO"]["333"]["MCRB_ECOLI"]["terms"]["MF"]["GO:0005524"] == "ATP binding"
    assert target["GO"]["333"]["MCRB_ECOLI"]["status"] == "normal"
    assert target["GO"]["333"]["MCRB_ECOLI"]["wang_sem_sim_bma_bp"] == 0.0
    assert target["GO"]["333"]["MCRB_ECOLI"]["wang_sem_sim_bma_mf"] == 0.541

def test_merge_nested_data_empty():
    """Test merging with empty source"""
    source = {}
    target = defaultdict(lambda: defaultdict(dict))
    merge_nested_data(source, target, "pos")
    assert dict(target) == {}


###T aggregate_range_positions - Integration Tests (+ merge_nested_data())

def test_aggregate_range_positions_annotations(mock_annotation_subset):
    """Test annotation data aggregation"""
    result = aggregate_range_positions(
        mock_annotation_subset,
        "DISULFID | Intrachain (with C-246); in linked form",
        (333, 333),
        "annotations"
    )

    print(json.dumps(result, indent=2))
    assert result["essentials"]["type"] == "DISULFID"
    assert "count" in result["essentials"]
    assert result["essentials"]["count"]["333"] == 1
    assert result["GO"]["333"]["MCRB_ECOLI"]["wang_sem_sim_bma_bp"] == 0.0


def test_aggregate_range_positions_conservations(mock_conservation_subset):
    """Test conservation data aggregation"""
    result = aggregate_range_positions(
        mock_conservation_subset,
        "conserved_positions",
        (364, 365),
        "conservations"
    )

    assert result["conservation"]["364"] == 0.90
    assert result["conservation"]["365"] == 0.90

    assert result["hit"]["364"] is True
    assert result["hit"]["365"] is True

    assert result["cons_residue"]["364"] == "G"
    assert result["target_residue"]["364"] == "G"
    assert result["cons_residue"]["365"] == "G"
    assert result["target_residue"]["364"] == "G"


def test_aggregate_range_positions_invalid_type(mock_aggregated_report_data):
    """Test with invalid data type"""
    interval_data = mock_aggregated_report_data["sp|Q9NU22|MDN1_HUMAN"]["PF07728"]["hit_intervals"]["325-451"]
    result = aggregate_range_positions(
        interval_data,
        "invalid_range",
        (333, 333),
        "invalid_type"
    )
    assert result == {}

###T transform_to_ranges - Integration Tests (+ aggregate_range_positions())

def test_transform_to_ranges_integration(mock_aggregated_report_data, multi_logger):
    """Test full range transformation"""
    interval = mock_aggregated_report_data["sp|Q9NU22|MDN1_HUMAN"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval, multi_logger)
    disulfid_data = result["data"]["annotations"]["DISULFID | Intrachain (with C-246); in linked form"]["333-333"]
    assert disulfid_data["essentials"]["type"] == "DISULFID"
    assert disulfid_data["essentials"]["count"]["333"] == 1
    assert disulfid_data["hit"]["333"] is True

    cons_data = result["data"]["conservations"]["conserved_positions"]["364-374"]

    assert cons_data["conservation"]["364"] == 0.90
    assert cons_data["conservation"]["374"] == 0.99

    assert cons_data["hit"]["364"] is True
    assert cons_data["hit"]["366"] is False  # Known miss position

    assert cons_data["cons_residue"]["364"] == "G"
    assert cons_data["cons_residue"]["374"] == "G"

    assert "364" in result["indices"]["conservations"]["matches"]
    assert "366" in result["indices"]["conservations"]["misses"]

def test_transform_to_ranges_full(mock_aggregated_report_data, multi_logger):
    """Integration test for range transformation"""
    interval_data = mock_aggregated_report_data["sp|Q9NU22|MDN1_HUMAN"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval_data, multi_logger)

    disulfid_id = "DISULFID | Intrachain (with C-246); in linked form"
    assert disulfid_id in result["data"]["annotations"]
    assert "333-333" in result["data"]["annotations"][disulfid_id]

    assert "conserved_positions" in result["data"]["conservations"]
    assert "364-374" in result["data"]["conservations"]["conserved_positions"]

def test_transform_to_ranges_basic(mock_aggregated_report_data, multi_logger):
    """Test basic range transformation functionality"""
    interval = mock_aggregated_report_data["sp|Q9NU22|MDN1_HUMAN"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval, multi_logger)

    disulfid_id = "DISULFID | Intrachain (with C-246); in linked form"
    assert disulfid_id in result["data"]["annotations"]
    assert "333-333" in result["data"]["annotations"][disulfid_id]

    anno_data = result["data"]["annotations"][disulfid_id]["333-333"]
    assert anno_data["essentials"]["type"] == "DISULFID"
    assert anno_data["hit"]["333"] is True

def test_transform_to_ranges_conservation(mock_aggregated_report_data, multi_logger):
    """Test conservation data transformation"""
    interval = mock_aggregated_report_data["sp|Q9NU22|MDN1_HUMAN"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval, multi_logger)

    assert "conserved_positions" in result["data"]["conservations"]
    assert "364-374" in result["data"]["conservations"]["conserved_positions"]

    cons_data = result["data"]["conservations"]["conserved_positions"]["364-374"]

    assert cons_data["conservation"]["364"] == 0.90
    assert cons_data["conservation"]["374"] == 0.99

    assert cons_data["hit"]["364"] is True
    assert cons_data["hit"]["366"] is False  # Known miss position

    assert cons_data["cons_residue"]["364"] == "G"
    assert cons_data["cons_residue"]["374"] == "G"

    assert "364" in result["indices"]["conservations"]["matches"]
    assert "366" in result["indices"]["conservations"]["misses"]

def test_transform_to_ranges_metadata(mock_aggregated_report_data, multi_logger):
    """Test metadata handling in transformation"""
    interval = mock_aggregated_report_data["sp|Q9NU22|MDN1_HUMAN"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval, multi_logger)

    assert "annotations" in result["data"]
    assert "333" in result["indices"]["annotations"]["matches"]

    assert "conservations" in result["data"]
    assert "conservations" in result["indices"]
    assert "364" in result["indices"]["conservations"]["matches"]
    assert "366" in result["indices"]["conservations"]["misses"]

def test_transform_to_ranges_error_handling(mock_aggregated_report_data, multi_logger):
    """Test error handling in range transformation"""
    interval = mock_aggregated_report_data["sp|Q9NU22|MDN1_HUMAN"]["PF07728"]["hit_intervals"]["325-451"]

    # Corrupt ranges data
    interval["annotation_ranges"]["DISULFID | Intrachain (with C-246); in linked form"]["ranges"] = "invalid"

    result = transform_to_ranges(interval, multi_logger)

    # Verify warning logged
    multi_logger.assert_called_with(
    'warning',
    'MAKE_VIEW - TRANSFORM_2_RANGES - %s invalid format for %s',
    'annotation_ranges', 'DISULFID | Intrachain (with C-246); in linked form'
)

    # Verify empty result for corrupted range
    assert "DISULFID | Intrachain (with C-246); in linked form" not in result

def test_transform_to_ranges_empty(multi_logger):
    """Test transformation with empty input"""
    result = transform_to_ranges({}, multi_logger)
    assert result == {}

###T process_sequence_report - Integration Tests (
# + convert_lists_to_original_types() + transform_to_ranges() + convert_sets_and_tuples_to_lists() + convert_default_dict_to_dict()
# )

def test_process_sequence_report_file_not_found(tmp_path, logger, multi_logger):
    """Test handling of missing input file"""
    sequence_id_mock = "sp|Q9NU22|MDN1_HUMAN"
    non_existent_path = tmp_path / "non_existent.json"

    with pytest.raises(FileNotFoundError):
        process_sequence_report(str(non_existent_path), sequence_id_mock, logger, multi_logger)

    assert multi_logger.call_count == 1
    call_args = multi_logger.call_args[0]
    assert call_args[0] == "error"  # First arg is log level
    assert call_args[1] == "MAKE_VIEW - PROC_SEQ_REP - Failed to open report file %s: %s"  # Message format
    assert call_args[2].endswith("/non_existent.json")  # File path ends correctly
    assert isinstance(call_args[3], FileNotFoundError)  # Error type is correct

def test_process_sequence_report_basic(tmp_path, mock_aggregated_report_data, logger, multi_logger):
    """Test full sequence report processing"""
    sequence_id_mock = "sp|Q9NU22|MDN1_HUMAN"
    report_path = tmp_path / "PF07728_report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(mock_aggregated_report_data, f)

    result = process_sequence_report(str(report_path), sequence_id_mock, logger, multi_logger)

    assert result["sequence_id"] == "sp|Q9NU22|MDN1_HUMAN"
    assert "PF07728" in result["domain"]
    assert "325-451" in result["domain"]["PF07728"]
    ranges = result["domain"]["PF07728"]["325-451"]
    assert "DISULFID | Intrachain (with C-246); in linked form" in ranges["data"]["annotations"]
    assert "conserved_positions" in ranges["data"]["conservations"]

###T write_range_views - Integration tests (+ convert_sets_and_tuples_to_lists())

def test_write_range_views_debug_logging(tmp_path, mock_transformed_data, logger, multi_logger):
    """Test debug logs for directory and file path info."""
    sequence_dir = tmp_path / "output" / "sp-Q9NU22-MDN1_HUMAN"
    clean_sequence_id_mock = "sp-Q9NU22-MDN1_HUMAN"
    sequence_dir.mkdir(parents=True)

    write_range_views(mock_transformed_data, str(sequence_dir), clean_sequence_id_mock, logger, multi_logger)

    debug_msgs = [str(call) for call in logger.debug.call_args_list]
    assert any("Writing data to:" in msg for msg in debug_msgs)


def test_write_range_views_permission_error(tmp_path, mock_transformed_data, logger, multi_logger):
    """Test handling of permission errors during writing"""
    # Create a readonly directory
    output_dir = tmp_path / "readonly"
    output_dir.mkdir()
    os.chmod(output_dir, 0o444)  # Read-only
    clean_sequence_id = "sp-Q9NU22-MDN1_HUMAN"

    with pytest.raises(PermissionError):
        write_range_views(mock_transformed_data, str(output_dir), clean_sequence_id, logger, multi_logger)

    # Check that multi_logger was called with error
    multi_logger.assert_called_with(
        "error",
        "MAKE VIEW - Failed to write range view for %s: %s",
        "sp-Q9NU22-MDN1_HUMAN",
        ANY  # Uses ANY to match the PermissionError which may have varying text
    )

def test_write_range_views_file_paths(tmp_path, mock_transformed_data, logger, multi_logger):
    """Test file path construction in write_range_views"""
    sequence_dir = tmp_path / "output" / "sp-Q9NU22-MDN1_HUMAN"
    sequence_dir.mkdir(parents=True)
    clean_sequence_id = "sp-Q9NU22-MDN1_HUMAN"

    write_range_views(mock_transformed_data, str(sequence_dir), clean_sequence_id, logger, multi_logger)

    range_file = sequence_dir / "PF07728_ranges.json"
    assert range_file.exists(), "Range file not created"

    logger.debug(f"Sequence dir: {sequence_dir}")
    logger.debug(f"Range file: {range_file}")

def test_write_range_views_json_content(tmp_path, mock_transformed_data, logger, multi_logger):
    """Test JSON content written by write_range_views"""
    sequence_dir = tmp_path / "output" / "sp-Q9NU22-MDN1_HUMAN"
    sequence_dir.mkdir(parents=True)
    clean_sequence_id = "sp-Q9NU22-MDN1_HUMAN"


    write_range_views(mock_transformed_data, str(sequence_dir), clean_sequence_id, logger, multi_logger)

    range_file = sequence_dir / "PF07728_ranges.json"

    with open(range_file, "r", encoding="utf-8") as f:
        content = json.load(f)

    # Check data structure preservation
    assert content["sequence_id"] == mock_transformed_data["sequence_id"]
    assert "domain" in content
    assert "PF07728" in content["domain"]

    # Log actual content for comparison
    logger.debug(f"Written content: {json.dumps(content, indent=2)}")
    logger.debug(f"Expected content: {json.dumps(mock_transformed_data, indent=2)}")

def test_write_range_views_data_transformation(tmp_path, logger, multi_logger):
    """Test data transformation during write"""
    # Create basic test data with sets/tuples
    test_data = {
        "sequence_id": "test",
        "domain": {
            "PF00001": {
                "test_set": {"a", "b"},
                "test_tuple": (1, 2)
            }
        }
    }

    sequence_dir = tmp_path / "output" / "sp-Q9NU22-MDN1_HUMAN"
    sequence_dir.mkdir(parents=True)
    clean_sequence_id = "sp-Q9NU22-MDN1_HUMAN"


    write_range_views(test_data, str(sequence_dir), clean_sequence_id, logger, multi_logger)
    range_file = sequence_dir / "PF00001_ranges.json"

    with open(range_file, "r", encoding="utf-8") as f:
        content = json.load(f)

    # Check conversion
    assert isinstance(content["domain"]["PF00001"]["test_set"], list)
    assert isinstance(content["domain"]["PF00001"]["test_tuple"], list)

    logger.debug(f"Transformed content: {json.dumps(content, indent=2)}")

def test_write_range_views(tmp_path, mock_transformed_data, logger, multi_logger):
    sequence_dir = tmp_path / "output" / "sp-Q9NU22-MDN1_HUMAN"
    sequence_dir.mkdir(parents=True)
    clean_sequence_id = "sp-Q9NU22-MDN1_HUMAN"


    write_range_views(mock_transformed_data, str(sequence_dir),  clean_sequence_id, logger, multi_logger)

    expected_file = sequence_dir / "PF07728_ranges.json"
    assert expected_file.exists()

    with open(expected_file, "r", encoding="utf-8") as f:
        written_data = json.load(f)
        assert written_data["sequence_id"] == "sp|Q9NU22|MDN1_HUMAN"

def test_process_and_write_integration(tmp_path, mock_aggregated_report_data, logger, multi_logger):
    """Test full workflow from processing to file writing"""
    # Setup
    sequence_id = "sp|Q9NU22|MDN1_HUMAN"
    clean_sequence_id = "sp-Q9NU22-MDN1_HUMAN"

    # Create sequence directory
    sequence_dir = tmp_path / "sp-Q9NU22-MDN1_HUMAN"
    os.makedirs(sequence_dir)

    # Write input report
    report_path = sequence_dir / "aggregated_report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(mock_aggregated_report_data, f)

    # Process
    transformed_data = process_sequence_report(str(report_path), sequence_id, logger, multi_logger)

    # Write views - Use sequence_dir instead of tmp_path
    write_range_views(transformed_data, str(sequence_dir), clean_sequence_id, logger, multi_logger)

    # Verify
    output_file = sequence_dir / "PF07728_ranges.json"
    assert output_file.exists(), f"Output file not found at {output_file}"

    with open(output_file, "r", encoding="utf-8") as f:
        result = json.load(f)

    # Verify structure
    assert "domain" in result
    assert "PF07728" in result["domain"]
    assert "325-451" in result["domain"]["PF07728"]

    # Verify content
    ranges = result["domain"]["PF07728"]["325-451"]
    assert "DISULFID | Intrachain (with C-246); in linked form" in ranges["data"]["annotations"]
    assert "conserved_positions" in ranges["data"]["conservations"]

    # Verify logging
    logger.info.assert_called()

###T main - Integration tests + everything

def test_complete_workflow_integration(tmp_path, mock_aggregated_report_data, logger, multi_logger):
    """Test complete workflow from processing to file writing with comprehensive verification"""
    # Setup
    sequence_id = "sp|Q9NU22|MDN1_HUMAN"
    clean_sequence_id = "sp-Q9NU22-MDN1_HUMAN"
    sequence_dir = tmp_path / "output" / clean_sequence_id
    os.makedirs(sequence_dir)

    # Write test data
    report_path = sequence_dir / "aggregated_report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(mock_aggregated_report_data, f)

    # Process the report
    transformed_data = process_sequence_report(str(report_path), sequence_id, logger, multi_logger)

    # Verify processing logs
    assert logger.info.call_args_list, "No logging calls were made"
    assert any("Starting to process report" in str(call) for call in logger.info.call_args_list)

    # Write the views
    write_range_views(transformed_data, str(sequence_dir), clean_sequence_id, logger, multi_logger)

    # Verify file creation
    output_file = sequence_dir / "PF07728_ranges.json"
    assert output_file.exists(), f"Output file not found at {output_file}"

    # Verify content
    with open(output_file, "r", encoding="utf-8") as f:
        result = json.load(f)

    # print("\n Here's the content written to JSON:\n", json.dumps(result, indent=2))

    # Structure verification
    assert result["sequence_id"] == sequence_id
    assert "PF07728" in result["domain"]
    assert "325-451" in result["domain"]["PF07728"]

    # Detailed content verification
    domain_data = result["domain"]["PF07728"]["325-451"]
    assert "DISULFID | Intrachain (with C-246); in linked form" in domain_data["data"]["annotations"]
    assert "conserved_positions" in domain_data["data"]["conservations"]
    assert "364-374" in domain_data["data"]["conservations"]["conserved_positions"]

    # Specific data point verification
    disulfid = domain_data["data"]["annotations"]["DISULFID | Intrachain (with C-246); in linked form"]["333-333"]
    assert disulfid["essentials"]["type"] == "DISULFID"
    assert "333" in disulfid["essentials"]["count"]

    # Log verification
    write_logs = [str(call) for call in logger.info.call_args_list]
    assert any("Writing range views for sequence" in log for log in write_logs)
    assert any("Successfully wrote range view" in log for log in write_logs)

def test_main_missing_report(tmp_path, multi_logger):
    """Test main() handling of missing aggregated_report.json"""
    sequence_dir = tmp_path / "sp-Q9NU22-MDN1_HUMAN"
    sequence_dir.mkdir()

    test_args = [
        "-sD", str(sequence_dir),
        "-s", "sp|Q9NU22|MDN1_HUMAN"
    ]

    with patch("make_view_jsons.get_multi_logger", return_value=multi_logger):
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("sys.argv", ["script_name.py"] + test_args)
            with pytest.raises(SystemExit) as exc_info:
                main()

    assert exc_info.value.code == 1

    # Get all calls made to multi_logger
    call_args_list = [call[0] for call in multi_logger.call_args_list]

    expected_call = (
        "error",
        "MAKE VIEW - No aggregated_report.json found at %s",
        str(sequence_dir / "aggregated_report.json")
    )
    assert expected_call in call_args_list

def test_main_empty_report(tmp_path, multi_logger):
    """Test main() handling of empty report"""
    sequence_dir = tmp_path / "test_sequence"
    sequence_dir.mkdir()

    # Create empty report
    report_path = sequence_dir / "aggregated_report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump({}, f)

    test_args = [
        "-sD", str(sequence_dir),
        "-s", "sp|Q9NU22|MDN1_HUMAN"
    ]

    with patch("make_view_jsons.get_multi_logger", return_value=multi_logger):
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("sys.argv", ["script_name.py"] + test_args)
            with pytest.raises(SystemExit) as exc_info:
                main()

    assert exc_info.value.code == 0

    logger_calls = multi_logger.call_args_list
    assert any(
        args[0] == ("warning", "MAKE_VIEW - PROC_SEQ_REP - Empty report for sequence %s - no domains were found", "sp|Q9NU22|MDN1_HUMAN")
        for args in logger_calls
    )

def test_main_invalid_json(tmp_path, multi_logger):
    """Test main() handling of invalid JSON"""
    sequence_dir = tmp_path / "test_sequence"
    sequence_dir.mkdir()

    # Create invalid JSON
    report_path = sequence_dir / "aggregated_report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("invalid json")

    test_args = [
        "-sD", str(sequence_dir),
        "-s", "sp|Q9NU22|MDN1_HUMAN"
    ]

    with patch("make_view_jsons.get_multi_logger", return_value=multi_logger):
        with pytest.MonkeyPatch().context() as mp:
            mp.setattr("sys.argv", ["script_name.py"] + test_args)
            with pytest.raises(SystemExit) as exc_info:
                main()

    assert exc_info.value.code == 1

    logger_calls = multi_logger.call_args_list
    # Need to check if any of the calls match since there may be multiple log messages
    assert any(
        args[0] == ("error", "Error processing aggregated report: %s", "Expecting value: line 1 column 1 (char 0)")
        for args in logger_calls
    )