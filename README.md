# Key-Residue-Annotate

A Python pipeline for genomic-scale protein function annotation, assigning both general and residue-specific information to novel proteins, while also leveraging conservation and GO term set comparisons to facilitate manual curation as a necessary follow-up step. The pipeline is aimed at assisting curators annotate and deposit proteins, providing a wide array of information that can be more easily verified than by standard procedures such as extensive unguided literature review.

## Installation

### Prerequisites

Download intermediary files from Zenodo, install InterProScan 5.73-104.0 and recreate conda env from YML in this repo.
- [Intermediary_Files-Zenodo](https://zenodo.org/records/15171019)
- [InterProScan](https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html)
- Recreate conda env from key_residue_annotate.yml

Alternatively:
- [Docker](https://www.docker.com/)
TO-DO: Container not yet made (as of 07/04/25)

### Install Docker

```bash
curl -fsSL https://get.docker.com | sudo sh
```

### Setup Docker

```bash
sudo groupadd docker
sudo usermod -aG docker $USER
docker login
```

## Usage

You may run the pipeline by defining a config.ini file or through command line, with the same arguments in either case. It is highly recommended to use a INI file, such as in the example for the H. sapiens reference proteome UP000005640:

```ini
[Inputs]
fasta = /home/eduardohorta/KRA/data/fasta/proteomes/HUMAN_UP000005640_9606_31_03_2025.fasta
hmm = /home/eduardohorta/KRA/resources/hmm/Pfam-A.hmm

[Paths]
iprscan_path = /home/eduardohorta/my_interproscan/interproscan-5.73-104.0/interproscan.sh
resource_dir = /home/eduardohorta/KRA/resources/
output_dir = /home/eduardohorta/KRA/git_repos/KRA/results/validation/HUMAN/
python = /home/eduardohorta/anaconda3/envs/key_residue_annotate/bin/python3
log = /home/eduardohorta/KRA/git_repos/KRA/logs/validation/HUMAN/executor_human.log

[Parameters]
output_format_iprscan = TSV
cpu_cores_iprscan = 11
number_jobs_iprscan = 1
seq_batch_size_iprscan = 2000
analyses_iprscan = panther,pfam,smart,gene3d,superfamily,prositepatterns,prositeprofiles,pirsf
enable_precalc_iprscan = True # Should be False in actual use with novel proteins
disable_res_iprscan = False
threads = 11
total_memory = 14
nucleotide = false
eco_codes =
```

Note that resource_dir should point to where you are keeping the intermediary files from Zenodo. Also from Zenodo, the pipeline will require both base Pfam-A.hmm and HMMPress-derived files (Pfam-A.hmm and Pfam-A.hmm.h3{p,m,i,f}).

Pleas consider that, while nucleotide FASTA input is supported (indicated by the nucleotide flag), it will be much slower than the expected amino acid input.

### Running with Nextflow

The same pipeline is available as a Nextflow workflow (main.nf), which is the recommended way to run it from now on. It gives the same outputs as executor.py on the same inputs. The main parameters mirror the INI file:

```bash
nextflow run main.nf -profile conda \
    --fasta proteome.fasta --hmm Pfam-A.hmm --resources resources/ --outdir results/ \
    --iprscan_path /path/to/interproscan.sh --cpu_cores_iprscan 8 --max_cpus 12 --max_memory 14.GB
```

Instead of running InterProScan, precomputed results can be given with `--iprscan_precomputed <dir>` (one `<sequence>/iprscan.tsv` per sequence, sequence names with `|` replaced by `-`), or InterProScan can be skipped with `--skip_iprscan true` (no GO comparison is made then). Other parameters: `--eco_codes`, `--bit_cutoffs` (gathering, trusted or noise), `--trim`, `--nucleotide` and `--extra_donors` (below).

### Extra donors

By default, annotations come only from the Swiss-Prot proteins in each domain's seed alignment. With `extra_donors = true` in the INI file (`--extra-donors` for executor.py, `--extra_donors` for Nextflow), Swiss-Prot reviewed proteins that carry the domain but are not in its seed are also used as donors, from each domain's `extra_donors.fasta` and `extra_donors.json` in resource_dir. This finds more annotations at a cost in precision. Every transferred annotation records which kind of donor supported it (`provenance`: seed, non_seed), and the PF*_ranges.json files carry a seed-only index (`annotations_seed`), so results can be filtered afterwards; annotations supported by both kinds of donor are the most reliable.

## Overview

executor.py: controller script for all scripts in the pipeline. The Nextflow workflow (main.nf, with its modules in modules/local/) runs the same scripts and is the recommended way to run the pipeline.

run_hmmsearch.py: runs PyHMMER's hmmsearch with the input FASTA. It translates nucleotides if needed, but at a heavy price in performance.

seq_and_batch_prep.py: creates a mapping JSON linking batches and sequence IDs. Creates individual directories for each of the latter and writes a FASTA containing the respective sequence in each directory. It also translates individual sequences from nucleotides, with the same performance cost. Both this and the preceding use the same translation method from PyHMMER.

run_iprscan.py: runs InterProScan in successive runs using batches delimited in the previous step. Represents an important connection point to other existing workflows that use InterProScan. We only use the GO terms from the TSV files internally, but the user may leverage this and other outputs (JSON, XML, GFF3) in downstream analyses.

prepare_fasta_per_domain.py: checks for intermediary files for each domain in resource_dir, if so, prepare a multifasta with all protein subsequence hits to it and put it in its directory.

run_hmmalign.py: runs HMMER3's hmmalign for each domain.

transfer_annotations.py: transfer annotations per domain from source/seed sequences from the domain's origin MSA to all novel protein subsequences that were hits to that domain. Concentrates the bulk of our custom processing.

merge_reports_in_sequences.py: aggregates all a sequence's PF*_report.json into a single aggregated_report.json in the same directory.

make_view_jsons.py: reformats the JSON structures for each sequence to be used by Nightingale and React (a pending task as of 07/04/2025).

utils.py: contains utility functions used throughout the pipeline, such as those involved in logging.

decorators.py: contains decorators used for future optimizations.

## Output

Currently, the output files include:

PF*_report.json: transferred annotations and metadata per domain for a sequence, containing only data for subsequences that were hits for that domain. Each annotation records its supporting donors and their provenance (seed or extra donors).

aggregated_report.json: groups the same data as the files above, but for the whole sequence and all hits across domains.

PF*_ranges.json: versions of the PF*_report.json files for consumption by Nightingale in a pending visualization step (as of 07/04/2025).

iprscan.tsv: internally used, but may be consumed by downstream analyses. The same can be said of any other requested formats (JSON, XML, GFF3) from InterProScan.



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15858571.svg)](https://doi.org/10.5281/zenodo.15858571)
