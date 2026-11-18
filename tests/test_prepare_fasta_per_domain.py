"""
Unit tests for prepare_fasta_per_domain.py
"""

import os
import sys
import json
from argparse import Namespace
from tempfile import TemporaryDirectory
from unittest.mock import patch, MagicMock
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from prepare_fasta_per_domain import (
    prep_domain_fasta,
    ExtraDonorsError,
    main
)

import pytest

domain_accession_mock = "PF07728"
target_header_mock = "sp|Q9NU22|MDN1_HUMANtarget/325-451"
extra_donors_fasta_mock = ">EXTRA1_HUMAN/10-50\nACDEFGHIK\n"


### Fixtures
@pytest.fixture
def logger():
    return MagicMock()

@pytest.fixture
def multi_logger():
    return MagicMock()

@pytest.fixture
def hits_per_domain_Q9NU22():
    """One hit of PF07728 in sp|Q9NU22|MDN1_HUMAN"""
    return {
        "PF07728": {
            "sp|Q9NU22|MDN1_HUMAN": [
                {"subseq": "VLLEGPIG", "ali_range": "325-451", "target_seq_name": "sp|Q9NU22|MDN1_HUMAN"}
            ]
        }
    }

def write_inputs(tmp_dir: str, hits_per_domain: dict, extra_donors_fasta: str = None) -> tuple:
    """Writes the hits per domain JSON and, if given, the domain's extra_donors.fasta.
    Returns the paths to the JSON, the resource dir and the output dir."""
    per_dom_json = os.path.join(tmp_dir, "hits_per_domain.json")
    resource_dir = os.path.join(tmp_dir, "resources")
    output_dir = os.path.join(tmp_dir, "out")
    os.makedirs(os.path.join(resource_dir, domain_accession_mock))
    os.makedirs(output_dir)
    with open(per_dom_json, "w", encoding="utf-8") as f:
        json.dump(hits_per_domain, f)
    if extra_donors_fasta is not None:
        with open(os.path.join(resource_dir, domain_accession_mock, "extra_donors.fasta"), "w", encoding="utf-8") as f:
            f.write(extra_donors_fasta)
    return per_dom_json, resource_dir, output_dir


###T prep_domain_fasta

def test_prep_domain_fasta_without_extra_donors(logger, multi_logger, hits_per_domain_Q9NU22):
    """Test that extra_donors.fasta is ignored when extra donors are not requested"""
    with TemporaryDirectory() as tmp_dir:
        per_dom_json, resource_dir, output_dir = write_inputs(tmp_dir, hits_per_domain_Q9NU22, extra_donors_fasta_mock)

        result = prep_domain_fasta(per_dom_json, domain_accession_mock, output_dir, logger, multi_logger,
                                   resource_dir=resource_dir, extra_donors=False)

        pure_path = os.path.join(output_dir, "PF07728", "PF07728_hits.fasta")
        assert result == pure_path
        with open(pure_path, encoding="utf-8") as f:
            assert f.read() == f">{target_header_mock}\nVLLEGPIG\n"
        assert not os.path.exists(os.path.join(output_dir, "PF07728", "PF07728_hits_with_donors.fasta"))

def test_prep_domain_fasta_with_extra_donors(logger, multi_logger, hits_per_domain_Q9NU22):
    """Test that the extra donors go to a separate FASTA and the hits FASTA keeps only the targets"""
    with TemporaryDirectory() as tmp_dir:
        per_dom_json, resource_dir, output_dir = write_inputs(tmp_dir, hits_per_domain_Q9NU22, extra_donors_fasta_mock)

        result = prep_domain_fasta(per_dom_json, domain_accession_mock, output_dir, logger, multi_logger,
                                   resource_dir=resource_dir, extra_donors=True)

        pure_path = os.path.join(output_dir, "PF07728", "PF07728_hits.fasta")
        with_donors_path = os.path.join(output_dir, "PF07728", "PF07728_hits_with_donors.fasta")
        assert result == with_donors_path
        with open(pure_path, encoding="utf-8") as f:
            assert f.read() == f">{target_header_mock}\nVLLEGPIG\n"
        with open(with_donors_path, encoding="utf-8") as f:
            assert f.read() == f">{target_header_mock}\nVLLEGPIG\n" + extra_donors_fasta_mock

def test_prep_domain_fasta_extra_donors_missing_newline(logger, multi_logger, hits_per_domain_Q9NU22):
    """Test that an extra_donors.fasta without a final newline still gives a well-formed FASTA"""
    with TemporaryDirectory() as tmp_dir:
        per_dom_json, resource_dir, output_dir = write_inputs(tmp_dir, hits_per_domain_Q9NU22, ">EXTRA1_HUMAN/10-50\nACDEFGHIK")

        result = prep_domain_fasta(per_dom_json, domain_accession_mock, output_dir, logger, multi_logger,
                                   resource_dir=resource_dir, extra_donors=True)

        with open(result, encoding="utf-8") as f:
            assert f.read() == f">{target_header_mock}\nVLLEGPIG\n>EXTRA1_HUMAN/10-50\nACDEFGHIK\n"

@pytest.mark.parametrize("extra_donors_fasta", [None, "", "\n  \n"])
def test_prep_domain_fasta_extra_donors_absent_or_empty(logger, multi_logger, hits_per_domain_Q9NU22, extra_donors_fasta):
    """Test that a domain with no extra donors (no file, or an empty one) uses the hits FASTA, without error"""
    with TemporaryDirectory() as tmp_dir:
        per_dom_json, resource_dir, output_dir = write_inputs(tmp_dir, hits_per_domain_Q9NU22, extra_donors_fasta)

        result = prep_domain_fasta(per_dom_json, domain_accession_mock, output_dir, logger, multi_logger,
                                   resource_dir=resource_dir, extra_donors=True)

        assert result == os.path.join(output_dir, "PF07728", "PF07728_hits.fasta")
        assert not os.path.exists(os.path.join(output_dir, "PF07728", "PF07728_hits_with_donors.fasta"))
        multi_logger.assert_not_called()

def test_prep_domain_fasta_extra_donors_write_failure(logger, multi_logger, hits_per_domain_Q9NU22):
    """Test that failing to write the FASTA with donors aborts instead of returning the hits FASTA.
    The failure is real: a directory sits where the file must be written."""
    with TemporaryDirectory() as tmp_dir:
        per_dom_json, resource_dir, output_dir = write_inputs(tmp_dir, hits_per_domain_Q9NU22, extra_donors_fasta_mock)
        os.makedirs(os.path.join(output_dir, "PF07728", "PF07728_hits_with_donors.fasta"))

        with pytest.raises(ExtraDonorsError) as error:
            prep_domain_fasta(per_dom_json, domain_accession_mock, output_dir, logger, multi_logger,
                              resource_dir=resource_dir, extra_donors=True)

        assert "PF07728" in str(error.value)
        assert "extra_donors.fasta" in str(error.value)
        level, message = multi_logger.call_args[0][0], multi_logger.call_args[0][1]
        assert level == "error"
        assert "FATAL" in message

@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0, reason="root can read files without permission")
def test_prep_domain_fasta_extra_donors_read_failure(logger, multi_logger, hits_per_domain_Q9NU22):
    """Test that failing to read extra_donors.fasta aborts. The failure is real: the file has no read permission."""
    with TemporaryDirectory() as tmp_dir:
        per_dom_json, resource_dir, output_dir = write_inputs(tmp_dir, hits_per_domain_Q9NU22, extra_donors_fasta_mock)
        extra_donors_path = os.path.join(resource_dir, "PF07728", "extra_donors.fasta")
        os.chmod(extra_donors_path, 0)
        try:
            with pytest.raises(ExtraDonorsError):
                prep_domain_fasta(per_dom_json, domain_accession_mock, output_dir, logger, multi_logger,
                                  resource_dir=resource_dir, extra_donors=True)
        finally:
            os.chmod(extra_donors_path, 0o644)


###T main

def test_main_extra_donors_points_domain_info_to_fasta_with_donors(hits_per_domain_Q9NU22):
    """Test that domain_info.json sends hmmalign to the FASTA with donors"""
    with TemporaryDirectory() as tmp_dir:
        per_dom_json, resource_dir, output_dir = write_inputs(tmp_dir, hits_per_domain_Q9NU22, extra_donors_fasta_mock)
        for filename in ("domain.hmm", "alignment.seed", "conservations.json"):
            open(os.path.join(resource_dir, "PF07728", filename), "w", encoding="utf-8").close()
        args = Namespace(per_dom_json=per_dom_json, domain_accession=domain_accession_mock,
                         resource_dir=resource_dir, output_dir=output_dir,
                         log=os.path.join(tmp_dir, "test.log"), extra_donors=True)

        with patch("prepare_fasta_per_domain.parse_arguments", return_value=args):
            main()

        with open(os.path.join(output_dir, "PF07728", "domain_info.json"), encoding="utf-8") as f:
            domain_info = json.load(f)
        assert domain_info["dom_fasta"] == os.path.join(output_dir, "PF07728", "PF07728_hits_with_donors.fasta")

def test_main_exits_on_extra_donors_failure():
    """Test that main exits with status 1 when the extra donors cannot be prepared"""
    args = Namespace(per_dom_json="hits.json", domain_accession=domain_accession_mock,
                     resource_dir="/resources", output_dir="/out", log="/logs/test.log", extra_donors=True)

    with patch("prepare_fasta_per_domain.parse_arguments", return_value=args), \
         patch("prepare_fasta_per_domain.get_logger", return_value=(MagicMock(), None)), \
         patch("prepare_fasta_per_domain.get_multi_logger", return_value=MagicMock()), \
         patch("prepare_fasta_per_domain.can_run_hmmalign", return_value={"can_align": True}), \
         patch("prepare_fasta_per_domain.prep_domain_fasta", side_effect=ExtraDonorsError("PF07728: failed")):
        with pytest.raises(SystemExit) as exit_info:
            main()

    assert exit_info.value.code == 1
