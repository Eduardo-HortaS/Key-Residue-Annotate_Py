"""
Integration tests running the whole pipeline on tests/fixture/

The fixture is small but real: four UniProt sequences and five Pfam 37.0 families picked so
that every case the pipeline must handle is present (Swiss-Prot and TrEMBL targets, a target
that is also a donor, a family without annotations.json, a conservations-only family,
donors with more than one region). See tests/fixture/SOURCES.txt.

InterProScan is not run: its per-sequence output is in the fixture and is placed where
executor.py expects it (run_iprscan.done marks the step as complete); the Nextflow pipeline
receives the same files through --iprscan_precomputed.
"""

import os
import sys
import json
import shutil
import subprocess

import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(TESTS_DIR)
FIXTURE_DIR = os.path.join(TESTS_DIR, "fixture")
FIXTURE_FASTA = os.path.join(FIXTURE_DIR, "proteome.fasta")
FIXTURE_HMM = os.path.join(FIXTURE_DIR, "Pfam-A.hmm")
FIXTURE_RESOURCES = os.path.join(FIXTURE_DIR, "resources")
FIXTURE_IPRSCAN = os.path.join(FIXTURE_DIR, "iprscan_by_sequence")

# Tools are looked up next to the running interpreter first (the conda env), then on PATH.
ENV_BIN = os.path.dirname(sys.executable)
RUN_ENV = dict(os.environ, PATH=ENV_BIN + os.pathsep + os.environ.get("PATH", ""))
HMMER_AVAILABLE = shutil.which("hmmalign", path=RUN_ENV["PATH"]) is not None
NEXTFLOW_AVAILABLE = shutil.which("nextflow", path=RUN_ENV["PATH"]) is not None

# Domains hit per sequence, with the fixture's Pfam-A.hmm and gathering cut-offs.
expected_hits_per_sequence = {
    "sp|Q32ZE1|POLG_ZIKV": ["PF01570", "PF20907"],
    "tr|A7K8Y4|A7K8Y4_9PHYC": ["PF01607"],
    "tr|A7K940|A7K940_9PHYC": ["PF01607"],
    "tr|A7K9F8|A7K9F8_9PHYC": ["PF01541", "PF07460"],
}

# Files executor.py keeps in its output tree and the Nextflow pipeline does not publish:
# its own copies of the inputs and the step sentinels.
executor_only_files = {
    "sequence.fasta", "iprscan.tsv",
    "run_iprscan.done", "prepare_fasta_per_domain.done", "run_hmmalign.done",
    "transfer_annotations.done", "merge_sequences.done", "make_view_jsons.done",
}

config_template = """[Inputs]
fasta = {fasta}
hmm = {hmm}

[Paths]
iprscan_path = /not/used/interproscan.sh
resource_dir = {resource_dir}
output_dir = {output_dir}
python = {python}
log = {log}

[Parameters]
output_format_iprscan = TSV
cpu_cores_iprscan = 3
number_jobs_iprscan = 1
seq_batch_size_iprscan = 2000
analyses_iprscan = pfam
enable_precalc_iprscan = True
disable_res_iprscan = False
threads = 4
total_memory = 6
nucleotide = false
eco_codes =
"""


def run_executor_on_fixture(run_dir: str) -> str:
    """Runs executor.py on the fixture inside run_dir and returns its output directory."""
    output_dir = os.path.join(run_dir, "out")
    resource_dir = os.path.join(run_dir, "resources")
    shutil.copytree(FIXTURE_RESOURCES, resource_dir)

    # Precomputed InterProScan output, where run_iprscan.py would have written it
    for sequence_dir in os.listdir(FIXTURE_IPRSCAN):
        os.makedirs(os.path.join(output_dir, sequence_dir))
        shutil.copyfile(os.path.join(FIXTURE_IPRSCAN, sequence_dir, "iprscan.tsv"),
                        os.path.join(output_dir, sequence_dir, "iprscan.tsv"))
    open(os.path.join(output_dir, "run_iprscan.done"), "w", encoding="utf-8").close()

    config_path = os.path.join(run_dir, "config.ini")
    with open(config_path, "w", encoding="utf-8") as f:
        f.write(config_template.format(
            fasta=FIXTURE_FASTA, hmm=FIXTURE_HMM, resource_dir=resource_dir + os.sep,
            output_dir=output_dir + os.sep, python=sys.executable,
            log=os.path.join(run_dir, "executor.log")))

    result = subprocess.run([sys.executable, "executor.py", "-c", config_path], cwd=REPO_DIR,
                            env=RUN_ENV, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    return output_dir


def run_nextflow_on_fixture(run_dir: str) -> str:
    """Runs main.nf on the fixture inside run_dir and returns its --outdir."""
    output_dir = os.path.join(run_dir, "results")
    resource_dir = os.path.join(run_dir, "resources")
    shutil.copytree(FIXTURE_RESOURCES, resource_dir)
    command = [
        "nextflow", "run", os.path.join(REPO_DIR, "main.nf"), "-ansi-log", "false",
        "-w", os.path.join(run_dir, "work"),
        "--fasta", FIXTURE_FASTA, "--hmm", FIXTURE_HMM, "--resources", resource_dir,
        "--iprscan_precomputed", FIXTURE_IPRSCAN, "--outdir", output_dir,
        "--threads", "1", "--max_cpus", "4", "--max_memory", "6.GB",
    ]
    result = subprocess.run(command, cwd=run_dir, env=RUN_ENV, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-2000:]
    return output_dir


def output_files(root: str) -> dict:
    """Maps ("top"|"domain"|"sequence", relative path) to the file's absolute path.

    executor.py writes <out>/PF*/ and <out>/<sequence>/; the Nextflow pipeline publishes
    PF*/ at the top and the sequences under sequences/.
    """
    sequences_root = os.path.join(root, "sequences")
    has_sequences_root = os.path.isdir(sequences_root)
    files = {}
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in ("pipeline_info", "batches")]
        for filename in filenames:
            path = os.path.join(dirpath, filename)
            parts = os.path.relpath(path, root).split(os.sep)
            if parts[0].startswith("PF"):
                files[("domain", os.path.join(*parts))] = path
            elif has_sequences_root and parts[0] == "sequences" and len(parts) > 2:
                files[("sequence", os.path.join(*parts[1:]))] = path
            elif not has_sequences_root and len(parts) > 1:
                files[("sequence", os.path.join(*parts))] = path
            else:
                files[("top", os.path.join(*parts))] = path
    return files


def same_file_paths(paths_a: dict, paths_b: dict) -> bool:
    """True when two domain_info.json name the same files, given as absolute or relative paths."""
    last_two = lambda v: "/".join(v.split("/")[-2:]) if isinstance(v, str) and "/" in v else v
    return set(paths_a) == set(paths_b) and all(last_two(paths_a[k]) == last_two(paths_b[k]) for k in paths_a)


def compare_outputs(executor_out: str, nextflow_out: str) -> dict:
    """Pairs and compares both output trees: JSON as data, everything else byte for byte."""
    executor_files = output_files(executor_out)
    nextflow_files = output_files(nextflow_out)
    comparison = {"identical": [], "different": [], "only_executor": [], "only_nextflow": []}
    for key in sorted(set(executor_files) | set(nextflow_files)):
        if key not in nextflow_files:
            comparison["only_executor"].append(key)
            continue
        if key not in executor_files:
            comparison["only_nextflow"].append(key)
            continue
        path_a, path_b = executor_files[key], nextflow_files[key]
        if path_a.endswith(".json"):
            with open(path_a, encoding="utf-8") as fa, open(path_b, encoding="utf-8") as fb:
                data_a, data_b = json.load(fa), json.load(fb)
            same = data_a == data_b
            if not same and os.path.basename(path_a) == "domain_info.json":
                same = same_file_paths(data_a, data_b)
        else:
            with open(path_a, "rb") as fa, open(path_b, "rb") as fb:
                same = fa.read() == fb.read()
        comparison["identical" if same else "different"].append(key)
    return comparison


### Fixtures
@pytest.fixture(scope="module")
def executor_output_dir(tmp_path_factory):
    """executor.py output on the fixture, shared by the tests in this module."""
    if not HMMER_AVAILABLE:
        pytest.skip("HMMER not found")
    return run_executor_on_fixture(str(tmp_path_factory.mktemp("executor")))


@pytest.fixture(scope="module")
def nextflow_output_dir(tmp_path_factory):
    """Nextflow pipeline output on the fixture, shared by the tests in this module."""
    if not (HMMER_AVAILABLE and NEXTFLOW_AVAILABLE):
        pytest.skip("HMMER or Nextflow not found")
    return run_nextflow_on_fixture(str(tmp_path_factory.mktemp("nextflow")))


### Tests
def test_executor_fixture_hits(executor_output_dir):
    """Every fixture sequence hits the domains it was chosen for, and nothing else"""
    with open(os.path.join(executor_output_dir, "hmmsearch_per_domain.json"), encoding="utf-8") as f:
        hits_per_domain = json.load(f)

    hits_per_sequence = {}
    for domain, sequences in hits_per_domain.items():
        for sequence in sequences:
            hits_per_sequence.setdefault(sequence, []).append(domain)

    assert {k: sorted(v) for k, v in hits_per_sequence.items()} == expected_hits_per_sequence


def test_executor_fixture_reports(executor_output_dir):
    """A report and a ranges file per hit domain, in every sequence's directory"""
    for sequence, domains in expected_hits_per_sequence.items():
        sequence_dir = os.path.join(executor_output_dir, sequence.replace("|", "-"))
        for domain in domains:
            assert os.path.isfile(os.path.join(sequence_dir, f"{domain}_report.json")), (sequence, domain)
            assert os.path.isfile(os.path.join(sequence_dir, f"{domain}_ranges.json")), (sequence, domain)
        assert os.path.isfile(os.path.join(sequence_dir, "aggregated_report.json")), sequence


def test_nextflow_reproduces_executor_on_fixture(executor_output_dir, nextflow_output_dir):
    """The Nextflow pipeline gives the same outputs as executor.py on the same inputs"""
    comparison = compare_outputs(executor_output_dir, nextflow_output_dir)

    assert comparison["different"] == []
    assert comparison["only_nextflow"] == []
    assert {os.path.basename(rel) for _, rel in comparison["only_executor"]} <= executor_only_files
    # Every domain and sequence output was actually compared
    compared = {rel for _, rel in comparison["identical"]}
    for sequence, domains in expected_hits_per_sequence.items():
        for domain in domains:
            assert os.path.join(sequence.replace("|", "-"), f"{domain}_report.json") in compared
            assert os.path.join(domain, f"{domain}_hmmalign.sth") in compared
