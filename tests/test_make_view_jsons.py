import logging
import json
import sys
import os
import copy
import pandas as pd
from io import StringIO
from importlib import reload
from argparse import Namespace
from collections import defaultdict
from unittest.mock import patch, ANY, call, mock_open, MagicMock
from tempfile import TemporaryDirectory
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
    def info_side_effect(msg, *args, **kwargs):
        formatted_msg = msg % args if args else msg
        print(f"INFO: {formatted_msg}")

    def error_side_effect(msg, *args, **kwargs):
        formatted_msg = msg % args if args else msg
        print(f"ERROR: {formatted_msg}")

    def warning_side_effect(msg, *args, **kwargs):
        formatted_msg = msg % args if args else msg
        print(f"WARNING: {formatted_msg}")

    def debug_side_effect(msg, *args, **kwargs):
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
def mock_report_data():
    return {
        "sequence_id": "sp|Q9NU22|MDN1_HUMAN",
        "domain": {
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
                                        "hit": True,
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
                                "364": {"conservation": 0.90, "hit": True},
                                "365": {"conservation": 0.90, "hit": True},
                                "366": {"conservation": 0.91, "hit": False},
                                "367": {"conservation": 0.92, "hit": True},
                                "368": {"conservation": 0.93, "hit": True},
                                "369": {"conservation": 0.94, "hit": True},
                                "370": {"conservation": 0.95, "hit": False},
                                "371": {"conservation": 0.96, "hit": True},
                                "372": {"conservation": 0.97, "hit": True},
                                "373": {"conservation": 0.98, "hit": True},
                                "374": {"conservation": 0.99, "hit": True},
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
def mock_annotation_subset():
    """Subset of mock_report_data focusing on annotations"""
    return {
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
                            "ECO:0000269|PubMed:12345678": {"count": 1}
                        },
                        "hit": True
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
    """Subset of mock_report_data focusing on conservations"""
    return {
        "conservations": {
            "positions": {
                "364": {"conservation": 0.90, "hit": True},
                "365": {"conservation": 0.90, "hit": True}
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

###T parse_arguments

def test_parse_arguments_required():
    output_dir_mock = "/home/user/results/views/"
    test_args = [
        "-o", output_dir_mock,
        "-s", "sp|Q9NU22|MDN1_HUMAN"
    ]

    with pytest.raises(SystemExit):
        parse_arguments()

    with pytest.MonkeyPatch().context() as mp:
        mp.setattr("sys.argv", ["script_name.py"] + test_args)
        args = parse_arguments()

    expected = Namespace(
        output_dir="/home/user/results/views/",
        sequence="sp|Q9NU22|MDN1_HUMAN",
        log="logs/make_views_jsons.log"
    )
    assert vars(args) == vars(expected)

def test_parse_arguments_optional():
    output_dir_mock = "/home/user/results/views/"
    log_filepath_mock = "/home/user/logs/make_views.log"

    test_args = [
        "-o", output_dir_mock,
        "-s", "sp|Q9NU22|MDN1_HUMAN",
        "-l", log_filepath_mock
    ]

    with pytest.MonkeyPatch().context() as mp:
        mp.setattr("sys.argv", ["script_name.py"] + test_args)
        args = parse_arguments()

    assert args.output_dir == output_dir_mock
    assert args.log == log_filepath_mock

###T merge_nested_data - Unit Tests

def test_merge_nested_data(mock_report_data):
    """Test nested data merging with position tracking"""
    source = mock_report_data["domain"]["PF07728"]["hit_intervals"]["325-451"]["annotations"]["positions"]["333"]["DISULFID | Intrachain (with C-246); in linked form"]
    target = defaultdict(lambda: defaultdict(dict))

    merge_nested_data(source, target, "333", "DISULFID | Intrachain (with C-246); in linked form")

    assert target["essentials"]["type"] == "DISULFID"
    assert "count" in target["essentials"]
    assert isinstance(target["essentials"]["count"], dict)
    assert target["essentials"]["count"]["333"] == 1

def test_merge_nested_data_empty():
    """Test merging with empty source"""
    source = {}
    target = defaultdict(lambda: defaultdict(dict))
    merge_nested_data(source, target, "pos", "range")
    assert dict(target) == {}


###T aggregate_range_positions - Unit Tests

def test_aggregate_range_positions_annotations(mock_annotation_subset):
    """Test annotation data aggregation"""
    result, metadata = aggregate_range_positions(
        mock_annotation_subset,
        "DISULFID | Intrachain (with C-246); in linked form",
        (333, 333),
        "annotations"
    )

    # Verify data structure
    assert result["essentials"]["type"] == "DISULFID"
    assert "count" in result["essentials"]
    assert result["essentials"]["count"]["333"] == 1

def test_aggregate_range_positions_conservations(mock_conservation_subset):
    """Test conservation data aggregation"""
    result, metadata = aggregate_range_positions(
        mock_conservation_subset,
        "conserved_positions",
        (364, 365),
        "conservations"
    )

    # Verify aggregation
    assert result["conservation"] == 0.90  # Average of [0.90, 0.90]
    assert result["hit"] is True  # All positions are hits

def test_aggregate_range_positions_invalid_type(mock_report_data):
    """Test with invalid data type"""
    interval_data = mock_report_data["domain"]["PF07728"]["hit_intervals"]["325-451"]
    result, metadata = aggregate_range_positions(
        interval_data,
        "invalid_range",
        (333, 333),
        "invalid_type"
    )
    assert result == {}
    assert metadata == {}

###T transform_to_ranges

def test_transform_to_ranges_integration(mock_report_data, multi_logger):
    """Test full range transformation"""
    interval = mock_report_data["domain"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval, multi_logger)

    print(json.dumps(result, indent=2))
    # Check annotation ranges
    disulfid_data = result["DISULFID | Intrachain (with C-246); in linked form"]["333-333"]
    assert disulfid_data["essentials"]["type"] == "DISULFID"
    assert "333_DISULFID | Intrachain (with C-246); in linked form" in disulfid_data["essentials"]["count"]

    # Check conservation ranges
    cons_data = result["conserved_positions"]["364-374"]
    assert 0.90 <= cons_data["conservation"] <= 0.99
    assert not cons_data["hit"]  # Some positions are misses

def test_transform_to_ranges_full(mock_report_data, multi_logger):
    """Integration test for range transformation"""
    interval_data = mock_report_data["domain"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval_data, multi_logger)

    # Check DISULFID annotation
    disulfid_id = "DISULFID | Intrachain (with C-246); in linked form"
    assert disulfid_id in result
    assert "333-333" in result[disulfid_id]

    # Check conservation data
    assert "conserved_positions" in result
    assert "364-374" in result["conserved_positions"]

    # Verify metadata placement
    assert "annotations" in result
    assert "conservations" in result


def test_transform_to_ranges_basic(mock_report_data, multi_logger):
    """Test basic range transformation functionality"""
    interval = mock_report_data["domain"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval, multi_logger)

    # Test correct annotation ID
    disulfid_id = "DISULFID | Intrachain (with C-246); in linked form"
    assert disulfid_id in result
    assert "333-333" in result[disulfid_id]

    # Test annotation content
    anno_data = result[disulfid_id]["333-333"]
    assert anno_data["essentials"]["type"] == "DISULFID"
    assert anno_data["hit"] is True


def test_transform_to_ranges_conservation(mock_report_data, multi_logger):
    """Test conservation data transformation"""
    interval = mock_report_data["domain"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval, multi_logger)

    # Test conservation range presence
    assert "conserved_positions" in result
    assert "364-374" in result["conserved_positions"]

    # Test conservation data
    cons_data = result["conserved_positions"]["364-374"]
    assert 0.90 <= cons_data["conservation"] <= 0.99
    assert cons_data["hit"] is False  # Contains misses

def test_transform_to_ranges_metadata(mock_report_data, multi_logger):
    """Test metadata handling in transformation"""
    interval = mock_report_data["domain"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval, multi_logger)

    # Test metadata placement
    assert "annotations" in result
    assert "333" in result["annotations"]["indices"]["matches"]

    assert "conservations" in result
    assert "indices" in result["conservations"]
    assert "364" in result["conservations"]["indices"]["matches"]
    assert "366" in result["conservations"]["indices"]["misses"]

def test_transform_to_ranges_error_handling(mock_report_data, multi_logger):
    """Test error handling in range transformation"""
    interval = mock_report_data["domain"]["PF07728"]["hit_intervals"]["325-451"]

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

###T process_sequence_report

def test_process_sequence_report_file_not_found(tmp_path, logger, multi_logger):
    """Test handling of missing input file"""
    non_existent_path = tmp_path / "non_existent.json"

    with pytest.raises(FileNotFoundError):
        process_sequence_report(str(non_existent_path), logger, multi_logger)

    logger.error.assert_called()

def test_process_sequence_report_integration(tmp_path, mock_report_data, logger, multi_logger):
    """Test full sequence report processing"""
    # Setup
    report_path = tmp_path / "PF07728_report.json"
    with open(report_path, "w") as f:
        json.dump(mock_report_data, f)

    result = process_sequence_report(str(report_path), logger, multi_logger)

    # Verify structure
    assert result["sequence_id"] == "sp|Q9NU22|MDN1_HUMAN"
    assert "PF07728" in result["domain"]
    assert "325-451" in result["domain"]["PF07728"]

    # Verify transformations
    ranges = result["domain"]["PF07728"]["325-451"]
    assert "DISULFID | Intrachain (with C-246); in linked form" in ranges
    assert "conserved_positions" in ranges

def test_process_sequence_report(tmp_path, mock_report_data, logger, multi_logger):
    report_path = tmp_path / "PF07728_report.json"
    with open(report_path, "w") as f:
        json.dump(mock_report_data, f)

    result = process_sequence_report(str(report_path), logger, multi_logger)

    assert result["sequence_id"] == "sp|Q9NU22|MDN1_HUMAN"
    assert "PF07728" in result["domain"]
    assert "325-451" in result["domain"]["PF07728"]


###T write_range_views

def test_write_range_views_debug_logging(tmp_path, mock_report_data, logger, multi_logger):
    """Test debug logs for directory and file path info."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    write_range_views(mock_report_data, str(output_dir), logger, multi_logger)

    debug_msgs = [str(call) for call in logger.debug.call_args_list]
    assert any("Created or verified directory:" in msg for msg in debug_msgs)
    assert any("Writing data to:" in msg for msg in debug_msgs)


def test_write_range_views_permission_error(tmp_path, mock_report_data, logger, multi_logger):
    """Test handling of permission errors during writing"""
    # Create a readonly directory
    output_dir = tmp_path / "readonly"
    output_dir.mkdir()
    os.chmod(output_dir, 0o444)  # Read-only

    with pytest.raises(PermissionError):
        write_range_views(mock_report_data, str(output_dir), logger, multi_logger)

    logger.error.assert_called()

def test_write_range_views_file_paths(tmp_path, mock_report_data, logger, multi_logger):
    """Test file path construction in write_range_views"""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    write_range_views(mock_report_data, str(output_dir), logger, multi_logger)

    # Check directory structure
    sequence_dir = output_dir / "sp|Q9NU22|MDN1_HUMAN"
    assert sequence_dir.exists(), "Sequence directory not created"

    # Check file creation
    range_file = sequence_dir / "PF07728_ranges.json"
    assert range_file.exists(), "Range file not created"

    # Log actual paths for debugging
    logger.debug(f"Output dir: {output_dir}")
    logger.debug(f"Sequence dir: {sequence_dir}")
    logger.debug(f"Range file: {range_file}")

def test_write_range_views_json_content(tmp_path, mock_report_data, logger, multi_logger):
    """Test JSON content written by write_range_views"""
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    write_range_views(mock_report_data, str(output_dir), logger, multi_logger)
    range_file = output_dir / "sp|Q9NU22|MDN1_HUMAN" / "PF07728_ranges.json"

    with open(range_file) as f:
        content = json.load(f)

    # Check data structure preservation
    assert content["sequence_id"] == mock_report_data["sequence_id"]
    assert "domain" in content
    assert "PF07728" in content["domain"]

    # Log actual content for comparison
    logger.debug(f"Written content: {json.dumps(content, indent=2)}")
    logger.debug(f"Expected content: {json.dumps(mock_report_data, indent=2)}")

def test_write_range_views_data_transformation(tmp_path, mock_report_data, logger, multi_logger):
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

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    write_range_views(test_data, str(output_dir), logger, multi_logger)
    range_file = output_dir / "test" / "PF00001_ranges.json"

    with open(range_file) as f:
        content = json.load(f)

    # Check conversion
    assert isinstance(content["domain"]["PF00001"]["test_set"], list)
    assert isinstance(content["domain"]["PF00001"]["test_tuple"], list)

    logger.debug(f"Transformed content: {json.dumps(content, indent=2)}")

def test_write_range_views(tmp_path, mock_report_data, logger, multi_logger):
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    write_range_views(mock_report_data, str(output_dir), logger, multi_logger)

    expected_file = output_dir / "sp|Q9NU22|MDN1_HUMAN" / "PF07728_ranges.json"
    assert expected_file.exists()

    with open(expected_file) as f:
        written_data = json.load(f)
        assert written_data["sequence_id"] == "sp|Q9NU22|MDN1_HUMAN"

def test_process_and_write_integration(tmp_path, mock_report_data, logger, multi_logger):
    """Test full workflow from processing to file writing"""
    # Setup
    report_path = tmp_path / "test_target" / "PF07728_report.json"
    os.makedirs(report_path.parent)
    with open(report_path, "w") as f:
        json.dump(mock_report_data, f)

    # Process
    transformed = process_sequence_report(str(report_path), logger, multi_logger)
    write_range_views(transformed, str(tmp_path), logger, multi_logger)

    # Verify
    output_file = tmp_path / "sp|Q9NU22|MDN1_HUMAN" / "PF07728_ranges.json"
    assert output_file.exists()

    with open(output_file) as f:
        result = json.load(f)

    # Verify structure
    assert "domain" in result
    assert "PF07728" in result["domain"]
    assert "325-451" in result["domain"]["PF07728"]

    # Verify content
    ranges = result["domain"]["PF07728"]["325-451"]
    assert "DISULFID | Intrachain (with C-246); in linked form" in ranges
    assert "conserved_positions" in ranges

    # Verify logging
    logger.info.assert_called()

###T main

def test_main_workflow_integration(tmp_path, mock_report_data, logger, multi_logger):
    """Test entire workflow from processing to file writing"""
    # Setup directory structure
    sequence_dir = tmp_path / "sp|Q9NU22|MDN1_HUMAN"
    os.makedirs(sequence_dir)
    report_path = sequence_dir / "PF07728_report.json"
    with open(report_path, "w") as f:
        json.dump(mock_report_data, f)

    # Run workflow
    transformed = process_sequence_report(str(report_path), logger, multi_logger)
    write_range_views(transformed, str(tmp_path), logger, multi_logger)

    # Verify output
    output_file = tmp_path / "sp|Q9NU22|MDN1_HUMAN" / "PF07728_ranges.json"
    assert output_file.exists()

    with open(output_file) as f:
        result = json.load(f)

    # Verify content structure and accuracy
    range_data = result["domain"]["PF07728"]["325-451"]
    assert "DISULFID | Intrachain (with C-246); in linked form" in range_data
    assert "conserved_positions" in range_data

    # Verify specific data points
    disulfid = range_data["DISULFID | Intrachain (with C-246); in linked form"]["333-333"]
    assert disulfid["essentials"]["type"] == "DISULFID"
    assert "333_DISULFID | Intrachain (with C-246); in linked form" in disulfid["essentials"]["count"]

def test_end_to_end_workflow(tmp_path, mock_report_data, logger, multi_logger):
    """Test the entire workflow with detailed logging verification"""
    # Setup
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    sequence_dir = tmp_path / "test_sequence"
    os.makedirs(sequence_dir)
    report_path = sequence_dir / "test_report.json"

    # Write test data
    with open(report_path, "w") as f:
        json.dump(mock_report_data, f)

    # Process the report
    transformed_data = process_sequence_report(str(report_path), logger, multi_logger)

    # Verify processing logs
    assert logger.info.call_args_list, "No logging calls were made"
    assert any("Starting to process report" in str(call) for call in logger.info.call_args_list)

    # Write the views
    write_range_views(transformed_data, str(output_dir), logger, multi_logger)

    # Check output file
    expected_file = output_dir / transformed_data["sequence_id"] / f"{list(transformed_data['domain'].keys())[0]}_ranges.json"
    assert expected_file.exists(), f"Output file not found at {expected_file}"

    # Verify file content
    with open(expected_file) as f:
        written_data = json.load(f)

    # Verify structure
    assert written_data["sequence_id"] == mock_report_data["sequence_id"]
    assert written_data["domain"] == transformed_data["domain"]

    # Verify writing logs
    write_logs = [str(call) for call in logger.info.call_args_list]
    assert any("Writing range views for sequence" in log for log in write_logs)
    assert any("Successfully wrote range view" in log for log in write_logs)