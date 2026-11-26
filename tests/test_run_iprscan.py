"""
Unit tests for run_iprscan.py
"""

import os
import sys
import shlex
from unittest.mock import patch, MagicMock
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from run_iprscan import run_interproscan, parse_arguments

import pytest

iprscan_path_mock = "/opt/interproscan/interproscan.sh"
input_fasta_mock = "batch_1.fasta"
output_basefile_mock = "batch_1"


def command_for(**kwargs) -> list:
    """The InterProScan command run_interproscan would run, split into its arguments."""
    arguments = {"iprscan_path": iprscan_path_mock, "enable_precalc": True, "disable_res": False,
                 "input_fasta": input_fasta_mock, "output_basefile": output_basefile_mock, "output_format": "TSV",
                 "analyses": "pfam,gene3d", "cpu_cores": 4, "multi_logger": MagicMock()}
    arguments.update(kwargs)
    with patch("run_iprscan.subprocess.run", return_value=MagicMock(returncode=0)) as mock_run:
        run_interproscan(**arguments)
    return shlex.split(mock_run.call_args[0][0])


###T run_interproscan

def test_run_interproscan_residue_annotations_by_default():
    """Test that by default residue annotations are requested (-etra) and not also disabled (-dra)"""
    command = command_for()
    assert "-etra" in command
    assert "-dra" not in command

def test_run_interproscan_disable_residue_annotations():
    """Test that disable_res gives -dra alone"""
    command = command_for(disable_res=True)
    assert "-dra" in command
    assert "-etra" not in command

def test_run_interproscan_several_output_formats():
    """Test that several formats reach InterProScan as a single comma-separated -f value, TSV included"""
    command = command_for(output_format="JSON, XML")
    assert command[command.index("-f") + 1] == "json,xml,tsv"

def test_run_interproscan_precalculated_lookup():
    """Test that -dp (no precalculated lookup) is passed only when the lookup is not enabled"""
    assert "-dp" not in command_for(enable_precalc=True)
    assert "-dp" in command_for(enable_precalc=False)

def test_run_interproscan_base_arguments():
    command = command_for()
    assert command[:7] == [iprscan_path_mock, "-i", input_fasta_mock, "-b", output_basefile_mock, "-f", "tsv"]
    assert command[command.index("-appl") + 1] == "pfam,gene3d"
    assert command[command.index("--cpu") + 1] == "4"
    assert "-goterms" in command and "-iprlookup" in command

def test_run_interproscan_returns_formats_with_tsv():
    """Test that the formats returned, which main splits per sequence, include the TSV added for GO terms"""
    with patch("run_iprscan.subprocess.run", return_value=MagicMock(returncode=0)):
        formats = run_interproscan(iprscan_path_mock, True, False, input_fasta_mock, output_basefile_mock,
                                   "JSON, XML", "pfam", 4, MagicMock())
    assert formats == ["json", "xml", "tsv"]


###T parse_arguments

def test_parse_arguments_single_mode():
    """Test that single mode parses without the batch-only arguments"""
    argv = ["run_iprscan.py", "-iPr", iprscan_path_mock, "single", "-f", input_fasta_mock, "-s", "seq1"]
    with patch.object(sys, "argv", argv):
        args = parse_arguments()
    assert args.mode == "single" and args.fasta == input_fasta_mock

def test_parse_arguments_batch_mode_requires_parent_dir():
    argv = ["run_iprscan.py", "-iPr", iprscan_path_mock, "batch", "-sB", "seq1", "-sBi", "1"]
    with patch.object(sys, "argv", argv), pytest.raises(SystemExit):
        parse_arguments()
