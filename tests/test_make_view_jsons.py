import logging
import json
import sys
import os
import copy
import pandas as pd
from io import StringIO
from importlib import reload
from argparse import Namespace
from unittest.mock import patch, ANY, call, mock_open, MagicMock
from tempfile import TemporaryDirectory
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from make_view_jsons import (
    parse_arguments,
    configure_logging,
    aggregate_range_annotations,
    transform_to_ranges,
    process_sequence_report,
    write_range_views,
    main
)

import pytest

### Fixtures

@pytest.fixture
def logger():
    """Create a properly configured mock logger with all necessary methods"""
    mock_logger = MagicMock()
    mock_logger.error = MagicMock()
    mock_logger.info = MagicMock()
    mock_logger.warning = MagicMock()
    mock_logger.debug = MagicMock()
    return mock_logger

@pytest.fixture
def mock_report_data():
    return {
        "sequence_id": "sp|Q9NU22|MDN1_HUMAN",
        "domain": {
            "PF07728": {
                "hit_intervals": {
                    "325-451": {
                        "sequence": "VLLEGPIGCGKTSLVEYLAAVTGRT...",
                        "length": 127,
                        "hit_start": 325,
                        "hit_end": 451,
                        "annotations": {
                            "positions": {
                                "333": {
                                    "DISULFID | Intrachain": {
                                        "essentials": {
                                            "type": "DISULFID",
                                            "description": "Intrachain",
                                            "count": 1,
                                            "annot_amino_acid": "C"
                                        },
                                        "evidence": {
                                            "ECO:0000269": {
                                                "count": 1
                                            }
                                        },
                                        "hit": True
                                    }
                                }
                            }
                        },
                        "conservations": {
                            "positions": {
                                "329": {"conservation": 0.9853, "hit": True}
                            },
                            "indices": {
                                "matches": ["329"],
                                "misses": []
                            }
                        },
                        "annotation_ranges": {
                            "DISULFID | Intrachain": {
                                "positions": [333],
                                "ranges": [[333, 333]]
                            }
                        },
                        "conservation_ranges": {
                            "conserved_positions": {
                                "positions": [329],
                                "ranges": [[329, 329]]
                            }
                        }
                    }
                }
            }
        }
    }


###T parse_arguments

def test_parse_arguments_required():
    output_dir_mock = "/home/user/results/views/"
    test_args = [
        "-o", output_dir_mock
    ]

    with pytest.raises(SystemExit):
        parse_arguments()

    with pytest.MonkeyPatch().context() as mp:
        mp.setattr("sys.argv", ["script_name.py"] + test_args)
        args = parse_arguments()

    expected = Namespace(
        output_dir="/home/user/results/views/",
        log="logs/make_views_jsons.log"
    )
    assert vars(args) == vars(expected)

def test_parse_arguments_optional():
    output_dir_mock = "/home/user/results/views/"
    log_filepath_mock = "/home/user/logs/make_views.log"

    test_args = [
        "-o", output_dir_mock,
        "-l", log_filepath_mock
    ]

    with pytest.MonkeyPatch().context() as mp:
        mp.setattr("sys.argv", ["script_name.py"] + test_args)
        args = parse_arguments()

    assert args.output_dir == output_dir_mock
    assert args.log == log_filepath_mock

###T configure_logging

def test_configure_logging(tmp_path):
    reload(logging)
    # Temporary log path for testing
    log_path = tmp_path / "test.log"
    logger = configure_logging(str(log_path))

    # Confirm the logger is set up
    assert logger is not None

    # Write something to create the file
    logger.info("Test log entry")

    # Now check if file exists
    assert log_path.exists(), f"Log file not found at {log_path}"

    # Verify contents
    with open(log_path, "r", encoding="utf-8") as log_file:
        log_contents = log_file.read()
    assert "Test log entry" in log_contents, "Log entry not found in log file"

###T aggregate_range_annotations

def test_aggregate_range_annotations():
    interval_dict = {
        "annotations": {
            "positions": {
                "333": {
                    "DISULFID | Intrachain": {
                        "essentials": {
                            "type": "DISULFID",
                            "count": 1
                        },
                        "evidence": {
                            "ECO:0000269": {"count": 1}
                        },
                        "hit": True
                    }
                }
            }
        }
    }

    result = aggregate_range_annotations(interval_dict, "DISULFID | Intrachain", (333, 333))

    assert result["essentials"]["type"] == "DISULFID"
    assert result["evidence"]["ECO:0000269"]["count"] == 1
    assert result["hit"] is True

###T transform_to_ranges

def test_transform_to_ranges(mock_report_data):
    interval = mock_report_data["domain"]["PF07728"]["hit_intervals"]["325-451"]
    result = transform_to_ranges(interval)

    # Test annotation ranges
    assert "DISULFID | Intrachain" in result
    assert "333-333" in result["DISULFID | Intrachain"]

    # Test conservation ranges
    assert "conserved_positions" in result
    assert "329-329" in result["conserved_positions"]


###T process_sequence_report

def test_process_sequence_report(tmp_path, mock_report_data):
    report_path = tmp_path / "PF07728_report.json"
    with open(report_path, "w") as f:
        json.dump(mock_report_data, f)

    result = process_sequence_report(str(report_path))

    assert result["sequence_id"] == "sp|Q9NU22|MDN1_HUMAN"
    assert "PF07728" in result["domain"]
    assert "325-451" in result["domain"]["PF07728"]


###T write_range_views

def test_write_range_views(tmp_path, mock_report_data, logger):
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    write_range_views(mock_report_data, str(output_dir), logger)

    expected_file = output_dir / "sp|Q9NU22|MDN1_HUMAN" / "PF07728_ranges.json"
    assert expected_file.exists()

    with open(expected_file) as f:
        written_data = json.load(f)
        assert written_data["sequence_id"] == "sp|Q9NU22|MDN1_HUMAN"

###T main
