"""
Unit tests for transfer_annotations.py
"""

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
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from transfer_annotations import (
    parse_arguments,
    get_pfam_id_from_hmmalign_result,
    get_annotation_filepath,
    read_files,
    iterate_aligned_sequences,
    extract_target_info_from_hmmalign,
    setup_for_conservations_only,
    find_and_map_annots,
    read_conservations_and_annotations,
    parse_go_annotations,
    check_interval_overlap,
    gather_go_terms_for_target,
    get_alignment_sequences,
    populate_conservation,
    populate_go_data_for_annotations,
    cleanup_improve_transfer_dict,
    convert_sets_and_tuples_to_lists,
    write_reports,
    map_and_filter_annot_pos,
    add_to_transfer_dict,
    _add_single_annotation,
    update_position_ranges,
    merge_adjacent_ranges,
    make_anno_total_dict,
    validate_annotations,
    process_annotation,
    validate_paired_annotations,
    main
)

from utils import get_logger

import pytest


hmmalign_result_mock = "/home/user/results/human/PF07728_hmmalign.sth"
hmmalign_result_content_mock = """# STOCKHOLM 1.0

#=GS VWA8_HUMAN/105-261                              AC A3KMH1.2
#=GS Q8ZSL8_PYRAE/20-171                             AC Q8ZSL8.1

VWA8_HUMAN/105-261                                      .........DVFLIGPPGPLRRSIAM.QYLELT..............KREVEYIALSR.DT..TETDLKQRREIR............AGTAFYIDQCAVRAAT..................EGRTLILEGLEKAE.R........N....VLP....V......LNN.LLENR.E...MQLEDGRFLMSAERYD.kLLRDhtkkelds............wkivrvsenFRVIALGLPVP........rYSGNPLDPPLRSRF............................
Q19346_CAEEL/78-234                                     .........DVFLLGVPGKIRLELVL.RYLEAT..............NREFEYLPITR.DT..TEADIKQRREIR............DGTAYYTDLCAVRAAL..................KGRVLVIDGVERAE.R........N....VLP....I......LNN.LLENR.E...MQLDDGRFLMKHDKYD.eLKVKydeatlkk............mgmervsenFHVIALGLPVP........rFPGNSLDPPFRSRF............................
MCRB_ECOLI/196-350                                      .........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYCPN............GVGFRRKDGIFYNFCQqak...........eqpeKKYIFIIDEINRANlS........K....VFG....E......VMM.LMEHD.K...RGENWSVPLTYSENDE.eRFYVp..........................enVYIIGLMNTAD.........RSLAVVDYALRRRF............................
sp|Q9NU22|MDN1_HUMANtarget//325-451                     .........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................FQFFAT-----.........--------------rrllscgg....................
#=GR sp|Q9NU22|MDN1_HUMANtarget//325-451             PP ..........79**************.******98..........77778999*****.**..*************..........*****************..................**************.*........*....***....*......***.*****.*...**9999999999996....................................666665............................56678999....................
sp|Q9NU22|MDN1_HUMANtarget//672-755                     .........PVLLVGETGTGKTSTIQ.YLAHIT..............GHRLRVVNMNQ.QS..DTADLLGGYKP-............-VDHKLIWLPLREAFE..................--------------.-........-....---....-......---.-----.-...----------------..----.............................-----------.........--------------elfaqtfskkqnftflghiqtc......
#=GR sp|Q9NU22|MDN1_HUMANtarget//672-755             PP .........69***************.******..............***********.**..*********96..............555555555555555................................................................................................................................................................5666666666666777777888......
tr|J3QRW1|J3QRW1_HUMANtarget//108-135                   .........GVLLYGPPGTGKTLLAR.AVAHHT..............DCT--------.--..------------............----------------..................--------------.-........-....---....-......---.-----.-...----------------..----.............................-----------.........--------------fi..........................
#=GR tr|J3QRW1|J3QRW1_HUMANtarget//108-135           PP .........59***************.****99..............665.....................................................................................................................................................................................................................69..........................
#=GC PP_cons                                            .........6799*************.******..............88889999999.99..899988888856............6888888888887888..................999********999.9........9....898....8......888.88888.8...8899988888888866..6555.............................99999999889.........78888888888887............................
#=GC RF                                                 .........xxxxxxxxxxxxxxxxx.xxxxxx..............xxxxxxxxxxx.xx..xxxxxxxxxxxx............xxxxxxxxxxxxxxxx..................xxxxxxxxxxxxxx.x........x....xxx....x......xxx.xxxxx.x...xxxxxxxxxxxxxxxx..xxxx.............................xxxxxxxxxxx.........xxxxxxxxxxxxxx............................
//
"""

resource_dir_mock = "/home/user/resources/"
domain_accession_mock = "PF07728"
output_dir_mock = "/home/user/results/human/PF07728/"
go_terms_mock = "/home/user/results/human/"
good_eco_codes_mock = ["ECO:0000269", "ECO:0000255", "ECO:0000313", "ECO:0007669"]
log_filepath_mock = "/home/user/logs/transfer_annotations.log"


### Fixtures
@pytest.fixture
def annotations_content_binding_fixture_Q9NU22_PF07728():
    """Base annotations with BINDING annotations and GO terms"""
    return {
        "MCRB_ECOLI": {
            "201": [
                {
                    "type": "BINDING",
                    "description": "Interacts with GTP",
                    "ligand_id": "ChEBI:CHEBI:37565",
                    "evidence": "ECO:0000255",
                    "entry": "P15005",
                    "aminoacid": "G"
                }
            ],
            "202": [
                {
                    "type": "BINDING",
                    "description": "Interacts with GTP",
                    "ligand_id": "ChEBI:CHEBI:37565",
                    "evidence": "ECO:0000255",
                    "entry": "P15005",
                    "aminoacid": "P"
                }
            ],
            "0": {
                "Molecular Function": {
                    "GO:0005524": "ATP binding",
                    "GO:0016887": "ATP hydrolysis activity",
                    "GO:0003677": "DNA binding",
                    "GO:0010385": "double-stranded methylated DNA binding",
                    "GO:0004519": "endonuclease activity",
                    "GO:0005525": "GTP binding",
                    "GO:0003924": "GTPase activity",
                    "GO:0044729": "hemi-methylated DNA-binding",
                    "GO:0042802": "identical protein binding",
                    "GO:0015666": "restriction endodeoxyribonuclease activity"
                },
                "Biological Process": {
                    "GO:0006308": "DNA catabolic process",
                    "GO:0009307": "DNA restriction-modification system"
                }
            }
        }
    }

@pytest.fixture
# NOTE: Intrachain doesn't really appear in actual descriptions, added for clarity in my tests.
# Instead, you might see Interchain and, in these cases, there'd be no paired_position key.
def annotations_content_disulfid_fixture_Q9NU22_PF07728():
    """Separate fixture for DISULFID pair annotations with GO terms"""
    return {
        "MCRB_ECOLI": {
            "205": [
                {
                    "type": "DISULFID",
                    "description": "Intrachain (with C-246); in linked form",
                    "evidence": "ECO:0000269|PubMed:12345678",
                    "entry": "P15005",
                    "aminoacid": "C",
                    "paired_position": "246"
                }
            ],
            "246": [
                {
                    "type": "DISULFID",
                    "description": "Intrachain (with C-205); in linked form",
                    "evidence": "ECO:0000269|PubMed:12345678",
                    "entry": "P15005",
                    "aminoacid": "C",
                    "paired_position": "205"
                }
            ],
            "0": {
                "Molecular Function": {
                    "GO:0005524": "ATP binding",
                    "GO:0016887": "ATP hydrolysis activity",
                    "GO:0003677": "DNA binding",
                    "GO:0010385": "double-stranded methylated DNA binding",
                    "GO:0004519": "endonuclease activity",
                    "GO:0005525": "GTP binding",
                    "GO:0003924": "GTPase activity",
                    "GO:0044729": "hemi-methylated DNA-binding",
                    "GO:0042802": "identical protein binding",
                    "GO:0015666": "restriction endodeoxyribonuclease activity"
                },
                "Biological Process": {
                    "GO:0006308": "DNA catabolic process",
                    "GO:0009307": "DNA restriction-modification system"
                }
            }
        }
    }

@pytest.fixture
def annotations_content_all_types_fixture_H0YB80_PF00244():
    """Fixture for all annotations for PF00244"""
    return {
    "14338_ARATH": {
        "70": [
            {
                "type": "MOD_RES",
                "description": "Phosphoserine",
                "evidence": "ECO:0000250|UniProtKB:P48349",
                "entry": "P48348",
                "aminoacid": "S"
            }
        ],
        "0": {
            "Biological Process": {
                "GO:0008104": "protein localization",
                "GO:0019222": "regulation of metabolic process",
                "GO:0050826": "response to freezing",
                "GO:0007165": "signal transduction"
            }
        },
        "112": [
            {
                "type": "MOD_RES",
                "description": "Phosphoserine",
                "evidence": "ECO:0000250|UniProtKB:P48349",
                "entry": "P48348",
                "aminoacid": "S"
            }
        ],
        "193": [
            {
                "type": "MOD_RES",
                "description": "Phosphoserine",
                "evidence": "ECO:0000250|UniProtKB:P48349",
                "entry": "P48348",
                "aminoacid": "S"
            }
        ],
        "214": [
            {
                "type": "MOD_RES",
                "description": "Phosphothreonine",
                "evidence": "ECO:0000250|UniProtKB:P48349",
                "entry": "P48348",
                "aminoacid": "T"
            }
        ]
    },
    "Q5CUW0_CRYPI": {
        "76": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:21853016",
                "entry": "Q5CUW0",
                "aminoacid": "K"
            }
        ],
        "0": {
            "Biological Process": {
                "GO:0007165": "signal transduction"
            }
        },
        "83": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:21853016",
                "entry": "Q5CUW0",
                "aminoacid": "R"
            }
        ],
        "156": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:21853016",
                "entry": "Q5CUW0",
                "aminoacid": "R"
            }
        ],
        "157": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:21853016",
                "entry": "Q5CUW0",
                "aminoacid": "Y"
            }
        ],
        "200": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:21853016",
                "entry": "Q5CUW0",
                "aminoacid": "L"
            }
        ],
        "201": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:21853016",
                "entry": "Q5CUW0",
                "aminoacid": "N"
            }
        ],
        "252": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:21853016",
                "entry": "Q5CUW0",
                "aminoacid": "N"
            }
        ],
        "256": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:21853016",
                "entry": "Q5CUW0",
                "aminoacid": "W"
            }
        ],
        "208": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:21853016",
                "entry": "Q5CUW0",
                "aminoacid": "E"
            }
        ]
    },
    "1433_GIAIC": {
        "135": [
            {
                "type": "BINDING",
                "description": "Interacts with O-phospho-L-serine",
                "ligand_id": "ChEBI:CHEBI:57524",
                "evidence": "ECO:0000269|PubMed:26551337, ECO:0007744|PDB:4ZQ0",
                "entry": "E2RU97",
                "aminoacid": "R",
                "paired_position": "136"
            }
        ],
        "0": {
            "Molecular Function": {
                "GO:0003779": "actin binding",
                "GO:0042802": "identical protein binding",
                "GO:0019900": "kinase binding",
                "GO:0051219": "phosphoprotein binding",
                "GO:0050815": "phosphoserine residue binding",
                "GO:0042803": "protein homodimerization activity"
            },
            "Biological Process": {
                "GO:0030036": "actin cytoskeleton organization",
                "GO:1990051": "activation of protein kinase C activity",
                "GO:0030010": "establishment of cell polarity",
                "GO:0000165": "MAPK cascade",
                "GO:0051495": "positive regulation of cytoskeleton organization",
                "GO:0051289": "protein homotetramerization",
                "GO:0070207": "protein homotrimerization",
                "GO:0008104": "protein localization",
                "GO:0097298": "regulation of nucleus size",
                "GO:0007165": "signal transduction"
            }
        },
        "136": [
            {
                "type": "BINDING",
                "description": "Interacts with O-phospho-L-serine",
                "ligand_id": "ChEBI:CHEBI:57524",
                "evidence": "ECO:0000269|PubMed:26551337, ECO:0007744|PDB:4ZQ0",
                "entry": "E2RU97",
                "aminoacid": "Y",
                "paired_position": "135"
            }
        ],
        "53": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:26551337",
                "entry": "E2RU97",
                "aminoacid": "K"
            },
            {
                "type": "SITE",
                "description": "Interaction with phosphoserine",
                "evidence": "ECO:0000269|PubMed:26551337, ECO:0007744|PDB:4ZQ0",
                "entry": "E2RU97",
                "aminoacid": "K"
            },
            {
                "type": "MUTAGEN",
                "description": "K->E: Loss or strongly decreased binding to synthetic human RAF1 phosphopeptides. Loss of binding to difopein.",
                "evidence": "ECO:0000269|PubMed:16368691, ECO:0000269|PubMed:19733174",
                "entry": "E2RU97",
                "aminoacid": "K"
            }
        ],
        "60": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:26551337",
                "entry": "E2RU97",
                "aminoacid": "R"
            },
            {
                "type": "SITE",
                "description": "Interaction with phosphoserine",
                "evidence": "ECO:0000269|PubMed:26551337, ECO:0007744|PDB:4ZQ0",
                "entry": "E2RU97",
                "aminoacid": "R"
            }
        ],
        "179": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:26551337",
                "entry": "E2RU97",
                "aminoacid": "L"
            }
        ],
        "180": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:26551337",
                "entry": "E2RU97",
                "aminoacid": "N"
            }
        ],
        "231": [
            {
                "type": "BINDING",
                "description": "Interacts with peptide",
                "ligand_ccd_id": "peptide",
                "evidence": "ECO:0000269|PubMed:26551337",
                "entry": "E2RU97",
                "aminoacid": "N"
            }
        ],
        "214": [
            {
                "type": "MOD_RES",
                "description": "Phosphothreonine",
                "evidence": "ECO:0000269|PubMed:16368691, ECO:0000269|PubMed:19733174, ECO:0000269|PubMed:24147113, ECO:0000269|PubMed:24658679",
                "entry": "E2RU97",
                "aminoacid": "T"
            },
            {
                "type": "MUTAGEN",
                "description": "T->A: Loss of phosphorylation by a protein kinase. No effect on subcellular localization. Dramatic decrease in the number of encysting parasites and cysts, but a large increase in the number of trophozoites. In encysting cells of 12 hours, significantly slower cyst conversion rate compared to the wild-type. No effect on binding to difopein. Decreased binding to a number of synthetic phosphopeptides.",
                "evidence": "ECO:0000269|PubMed:16368691, ECO:0000269|PubMed:19733174",
                "entry": "E2RU97",
                "aminoacid": "T"
            }
        ],
        "183": [
            {
                "type": "MUTAGEN",
                "description": "V->D: Loss of binding to difopein.",
                "evidence": "ECO:0000269|PubMed:19733174",
                "entry": "E2RU97",
                "aminoacid": "V"
            }
        ],
        "200": [
            {
                "type": "MUTAGEN",
                "description": "R->K: Increased oligomerization.",
                "evidence": "ECO:0000269|PubMed:24658679",
                "entry": "E2RU97",
                "aminoacid": "R"
            }
        ],
        "208": [
            {
                "type": "MUTAGEN",
                "description": "T->A: Slightly decreased oligomerization.",
                "evidence": "ECO:0000269|PubMed:24658679",
                "entry": "E2RU97",
                "aminoacid": "T"
            }
        ]
    }
}

@pytest.fixture
def conservations_content_Q9NU22_PF07728():
    return {
    "Q7US48_RHOBA/138-284": {
        "143": {
            "conservation": 0.9853,
            "amino_acid": "G"
        },
        "146": {
            "conservation": 0.8806,
            "amino_acid": "G"
        },
        "148": {
            "conservation": 0.9084,
            "amino_acid": "G"
        },
        "149": {
            "conservation": 0.9562,
            "amino_acid": "K"
        },
        "267": {
            "conservation": 0.8021,
            "amino_acid": "N"
        },
        "283": {
            "conservation": 0.9596,
            "amino_acid": "R"
        },
        "284": {
            "conservation": 0.8453,
            "amino_acid": "F"
        }
    }
}

@pytest.fixture
def conservations_content_H0YB80_PF00244():
    return {
        "C1E2K1_MICCC/19-248": {
            "20": {
                "conservation": 0.9802880970432145,
                "amino_acid": "G"
            },
            "21": {
                "conservation": 0.8013100436681223,
                "amino_acid": "K"
            },
            "22": {
                "conservation": 0.8363,
                "amino_acid": "L"
            },
            "23": {
                "conservation": 0.8231,
                "amino_acid": "A"
            },
            "24": {
                "conservation": 0.8863,
                "amino_acid": "E"
            },
            "25": {
                "conservation": 0.8405,
                "amino_acid": "Q"
            },
            "26": {
                "conservation": 0.8197,
                "amino_acid": "A"
            },
            "27": {
                "conservation": 0.8922,
                "amino_acid": "E"
            },
            "28": {
                "conservation": 0.95,
                "amino_acid": "R"
            },
            "29": {
                "conservation": 0.8626,
                "amino_acid": "Y"
            },
            "32": {
                "conservation": 0.9008,
                "amino_acid": "M"
            },
            "36": {
                "conservation": 0.8561,
                "amino_acid": "M"
            },
            "51": {
                "conservation": 0.9382,
                "amino_acid": "L"
            },
            "54": {
                "conservation": 0.836,
                "amino_acid": "E"
            },
            "55": {
                "conservation": 0.9296,
                "amino_acid": "E"
            },
            "56": {
                "conservation": 0.9482,
                "amino_acid": "R"
            },
            "57": {
                "conservation": 0.8818,
                "amino_acid": "N"
            },
            "58": {
                "conservation": 0.9157,
                "amino_acid": "L"
            },
            "59": {
                "conservation": 0.8792,
                "amino_acid": "L"
            },
            "60": {
                "conservation": 0.9061,
                "amino_acid": "S"
            },
            "61": {
                "conservation": 0.8446,
                "amino_acid": "V"
            },
            "62": {
                "conservation": 0.8378,
                "amino_acid": "A"
            },
            "63": {
                "conservation": 0.896,
                "amino_acid": "Y"
            },
            "64": {
                "conservation": 0.9528,
                "amino_acid": "K"
            },
            "65": {
                "conservation": 0.8865,
                "amino_acid": "N"
            },
            "68": {
                "conservation": 0.8355,
                "amino_acid": "G"
            },
            "70": {
                "conservation": 0.8069,
                "amino_acid": "R"
            },
            "71": {
                "conservation": 0.97,
                "amino_acid": "R"
            },
            "74": {
                "conservation": 0.8662,
                "amino_acid": "W"
            },
            "75": {
                "conservation": 0.9293,
                "amino_acid": "R"
            },
            "81": {
                "conservation": 0.8133,
                "amino_acid": "E"
            },
            "83": {
                "conservation": 0.8135,
                "amino_acid": "K"
            },
            "101": {
                "conservation": 0.8402,
                "amino_acid": "Y"
            },
            "106": {
                "conservation": 0.8636,
                "amino_acid": "E"
            },
            "108": {
                "conservation": 0.9349,
                "amino_acid": "E"
            },
            "109": {
                "conservation": 0.8949,
                "amino_acid": "L"
            },
            "113": {
                "conservation": 0.9358,
                "amino_acid": "C"
            },
            "124": {
                "conservation": 0.9197,
                "amino_acid": "L"
            },
            "132": {
                "conservation": 0.8321,
                "amino_acid": "E"
            },
            "135": {
                "conservation": 0.8855,
                "amino_acid": "V"
            },
            "136": {
                "conservation": 0.9213,
                "amino_acid": "F"
            },
            "137": {
                "conservation": 0.8509,
                "amino_acid": "Y"
            },
            "139": {
                "conservation": 0.943,
                "amino_acid": "K"
            },
            "140": {
                "conservation": 0.8856,
                "amino_acid": "M"
            },
            "141": {
                "conservation": 0.8162,
                "amino_acid": "K"
            },
            "142": {
                "conservation": 0.8751,
                "amino_acid": "G"
            },
            "143": {
                "conservation": 0.986,
                "amino_acid": "D"
            },
            "144": {
                "conservation": 0.9316,
                "amino_acid": "Y"
            },
            "146": {
                "conservation": 0.9776,
                "amino_acid": "R"
            },
            "147": {
                "conservation": 0.9815,
                "amino_acid": "Y"
            },
            "149": {
                "conservation": 0.8356,
                "amino_acid": "A"
            },
            "150": {
                "conservation": 0.9682,
                "amino_acid": "E"
            },
            "168": {
                "conservation": 0.9736,
                "amino_acid": "Y"
            },
            "171": {
                "conservation": 0.9657,
                "amino_acid": "A"
            },
            "175": {
                "conservation": 0.812,
                "amino_acid": "A"
            },
            "182": {
                "conservation": 0.8392,
                "amino_acid": "T"
            },
            "184": {
                "conservation": 0.9312,
                "amino_acid": "P"
            },
            "186": {
                "conservation": 0.9043,
                "amino_acid": "R"
            },
            "187": {
                "conservation": 0.9743,
                "amino_acid": "L"
            },
            "188": {
                "conservation": 0.9191,
                "amino_acid": "G"
            },
            "189": {
                "conservation": 0.9228,
                "amino_acid": "L"
            },
            "190": {
                "conservation": 0.8745,
                "amino_acid": "A"
            },
            "191": {
                "conservation": 0.978,
                "amino_acid": "L"
            },
            "192": {
                "conservation": 0.9648,
                "amino_acid": "N"
            },
            "194": {
                "conservation": 0.9051,
                "amino_acid": "S"
            },
            "195": {
                "conservation": 0.9377,
                "amino_acid": "V"
            },
            "196": {
                "conservation": 0.9346,
                "amino_acid": "F"
            },
            "198": {
                "conservation": 0.8826,
                "amino_acid": "Y"
            },
            "199": {
                "conservation": 0.9215,
                "amino_acid": "E"
            },
            "200": {
                "conservation": 0.8339,
                "amino_acid": "I"
            },
            "207": {
                "conservation": 0.9722,
                "amino_acid": "A"
            },
            "208": {
                "conservation": 0.8229,
                "amino_acid": "C"
            },
            "211": {
                "conservation": 0.9375,
                "amino_acid": "A"
            },
            "214": {
                "conservation": 0.8853,
                "amino_acid": "A"
            },
            "215": {
                "conservation": 0.87,
                "amino_acid": "F"
            },
            "216": {
                "conservation": 0.8009,
                "amino_acid": "D"
            },
            "218": {
                "conservation": 0.909,
                "amino_acid": "A"
            },
            "222": {
                "conservation": 0.8463,
                "amino_acid": "L"
            },
            "223": {
                "conservation": 0.8356,
                "amino_acid": "D"
            },
            "227": {
                "conservation": 0.8168,
                "amino_acid": "E"
            },
            "230": {
                "conservation": 0.8642,
                "amino_acid": "Y"
            },
            "231": {
                "conservation": 0.8046,
                "amino_acid": "K"
            },
            "232": {
                "conservation": 0.8844,
                "amino_acid": "D"
            },
            "233": {
                "conservation": 0.8355,
                "amino_acid": "S"
            },
            "234": {
                "conservation": 0.832,
                "amino_acid": "T"
            },
            "235": {
                "conservation": 0.8149,
                "amino_acid": "L"
            },
            "236": {
                "conservation": 0.9163,
                "amino_acid": "I"
            },
            "237": {
                "conservation": 0.8082,
                "amino_acid": "M"
            },
            "238": {
                "conservation": 0.8983,
                "amino_acid": "Q"
            },
            "239": {
                "conservation": 0.9404,
                "amino_acid": "L"
            },
            "240": {
                "conservation": 0.9178,
                "amino_acid": "L"
            },
            "241": {
                "conservation": 0.8887,
                "amino_acid": "R"
            },
            "242": {
                "conservation": 0.9231,
                "amino_acid": "D"
            },
            "243": {
                "conservation": 0.9685,
                "amino_acid": "N"
            },
            "244": {
                "conservation": 0.9009,
                "amino_acid": "L"
            },
            "245": {
                "conservation": 0.8181,
                "amino_acid": "T"
            },
            "246": {
                "conservation": 0.904,
                "amino_acid": "L"
            },
            "247": {
                "conservation": 0.9993,
                "amino_acid": "W"
            },
            "248": {
                "conservation": 0.8651,
                "amino_acid": "T"
            }
    }
}

@pytest.fixture
def mapping_content_Q9NU22_and_H0YB80_domains():
    """Fixture for mapping content of Q9NU22 and H0YB80 domains."""
    return "InterPro_ID\tPfam_ID\nIPR011704\tPF07728\nIPR023410\tPF00244\nIPR003959\tPF00004"

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
def multi_logger():
    """Create a mock multi_logger that matches the get_multi_logger behavior"""
    def multi_log(level: str, message: str, *args):
        # This simulates the behavior of the actual multi_logger
        # which formats the message with args and logs to multiple loggers
        return getattr(logging.getLogger(), level)(message, *args)

    mock = MagicMock(side_effect=multi_log)
    return mock

@pytest.fixture
def good_eco_codes_all():
    return ["ECO:0000255", "ECO:0000269", "ECO:0000313", "ECO:0007669"]

@pytest.fixture
def target_sequence_Q9NU22():
    return ".........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................FQFFAT-----.........--------------rrllscgg...................."

@pytest.fixture
def target_sequence_continuous_Q9NU22(target_sequence_Q9NU22):
    target_sequence_to_process = target_sequence_Q9NU22
    target_sequence_continuous_Q9NU22 = "".join([char.upper() for char in target_sequence_to_process if char.isalpha()])
    return target_sequence_continuous_Q9NU22

@pytest.fixture
def annot_sequence_Q9NU22():
    return ".........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYCPN............GVGFRRKDGIFYNFCQqak...........eqpeKKYIFIIDEINRANlS........K....VFG....E......VMM.LMEHD.K...RGENWSVPLTYSENDE.eRFYVp..........................enVYIIGLMNTAD.........RSLAVVDYALRRRF............................"


@pytest.fixture
def base_annotations():
    """Base structure with BINDING annotations only"""
    return {
        "201": [
            {
                "type": "BINDING",
                "description": "Interacts with GTP",
                "ligand_id": "ChEBI:CHEBI:37565",
                "evidence": "ECO:0000255",
                "entry": "P15005",
                "aminoacid": "G"
            }
        ],
        "202": [
            {
                "type": "BINDING",
                "description": "Interacts with GTP",
                "ligand_id": "ChEBI:CHEBI:37565",
                "evidence": "ECO:0000255",
                "entry": "P15005",
                "aminoacid": "P"
            }
        ]
    }

@pytest.fixture
def entry_annotations_binding_only(base_annotations):
    """Just BINDING annotations"""
    return base_annotations

@pytest.fixture
def entry_annotations_disulfid_pair():
    """Just DISULFID pair annotations"""
    return {
        # ADDED FOR TESTING PAIRED ANNOTATIONS, MODIFIED ORIGINAL SEQUENCE
        "205": [
            {
                "type": "DISULFID",
                "description": "Intrachain (with C-246); in linked form",
                "evidence": "ECO:0000269|PubMed:12345678",
                "entry": "P15005",
                "aminoacid": "C", # NIILQGPPGC
                "paired_position": "246"
            }],
        "246": [
            {
                "type": "DISULFID",
                "description": "Intrachain (with C-205); in linked form",
                "evidence": "ECO:0000269|PubMed:12345678",
                "entry": "P15005",
                "aminoacid": "C", # SYEDFIQGYRCN
                "paired_position": "205"
            }
        ]
    }

@pytest.fixture
def transfer_dict():
    return {}

@pytest.fixture
def transfer_dict_initialized_structure_Q9NU22():
    """Fixture for transfer_dict with structure initialized but no annotation data"""
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
                            "annotations": {"positions": {}, "indices": {"matches": set(), "misses": set()}},
                            "conservations": {"positions": {}, "indices": {"matches": set(), "misses": set()}},
                            "position_conversion": {
                                "target_to_aln": {},
                                "aln_to_target": {}
                            },
                            "annotation_ranges": {},
                            "conservation_ranges": {}
                        }
                    }
                }
            }
        }
    }

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
                            },
                            "conservation_ranges": {}
                        }
                    }
                }
            }
        }
    }

@pytest.fixture
def transfer_dict_populated_pre_process_Q9NU22(transfer_dict_populated_disulfid_Q9NU22):
    """Pre-processed version with modified structure"""
    original = copy.deepcopy(transfer_dict_populated_disulfid_Q9NU22)
    content = original["DOMAIN"]["sequence_id"]

    return {
        "domain": {
            "PF07728": {
                "sequence_id": content
            }
        }
    }

@pytest.fixture
def transfer_dict_populated_disulfid_list_Q9NU22(transfer_dict_populated_pre_process_Q9NU22):
    """Same dictionary but with lists instead of sets"""
    return convert_sets_and_tuples_to_lists(copy.deepcopy(transfer_dict_populated_pre_process_Q9NU22))

@pytest.fixture
def transfer_dict_populated_disulfid_post_conservation_Q9NU22(transfer_dict_populated_pre_process_Q9NU22):
    """Creates an extended version of transfer_dict_populated_disulfid_Q9NU22 with GO and conservation data."""
    extended_dict = copy.deepcopy(transfer_dict_populated_pre_process_Q9NU22)
    interval_path = extended_dict["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]

    # Update positions and scores
    interval_path["conservations"]["positions"].update({
        "329": {"conservation": 0.9853, "hit": True},
        "332": {"conservation": 0.8806, "hit": True},
        "334": {"conservation": 0.9084, "hit": True},
        "335": {"conservation": 0.9562, "hit": True}
    })
    interval_path["conservations"]["indices"]["matches"].update(["329", "332", "334", "335"])

    # Add conservation positions to existing position mappings
    conservation_mappings = {
        "target_to_aln": {
            "329": "14", "332": "17", "334": "19", "335": "20"
        },
        "aln_to_target": {
            "14": "329", "17": "332", "19": "334", "20": "335"
        }
    }

    # Merge conservation mappings with existing annotation mappings
    interval_path["position_conversion"]["target_to_aln"].update(conservation_mappings["target_to_aln"])
    interval_path["position_conversion"]["aln_to_target"].update(conservation_mappings["aln_to_target"])

    # Update conservation ranges
    interval_path["conservation_ranges"] = {
        "conserved_positions": {
            "positions": {329, 332, 334, 335},
            "ranges": [(329, 329), (332, 332), (334, 335)]
        }
    }

    return extended_dict

@pytest.fixture
def transfer_dict_populated_disulfid_post_gos_Q9NU22(transfer_dict_populated_disulfid_post_conservation_Q9NU22):
    """Adds GO terms data to the transfer dict"""
    extended_dict = copy.deepcopy(transfer_dict_populated_disulfid_post_conservation_Q9NU22)

    # Path to annotations
    interval_path = extended_dict["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]

    # Add GO data to both positions
    for pos in ["333", "373"]:
        anno_key = list(interval_path["annotations"]["positions"][pos].keys())[0]
        interval_path["annotations"]["positions"][pos][anno_key]["GO"] = {
            "MCRB_ECOLI": {
                "jaccard_index": 0.1667,
                "terms": {
                    "GO:0005524": "ATP binding",
                    "GO:0016887": "ATP hydrolysis activity"
                }
            }
        }

    return extended_dict

@pytest.fixture
def transfer_dict_populated_disulfid_post_gos_list_Q9NU22(transfer_dict_populated_disulfid_post_gos_Q9NU22):
    """Post-processed version with lists and GO terms matching expected output structure"""
    transfer_dict_converted = convert_sets_and_tuples_to_lists(copy.deepcopy(transfer_dict_populated_disulfid_post_gos_Q9NU22))

    hit_intervals_content = transfer_dict_converted["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]

    final_dict = {
        "sequence_id": "sp|Q9NU22|MDN1_HUMAN",
        "domain": {
            "PF07728": {
                "hit_intervals": hit_intervals_content
            }
        }
    }

    return final_dict

@pytest.fixture
def transfer_dict_success_binding_Q9NU22():
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
                            "329": {
                                "BINDING | Interacts with GTP": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with GTP",
                                        "count": 1,
                                        "annot_amino_acid": "G",
                                        "target_amino_acid": "G"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000255": {
                                            "rep_primary_accession": "P15005",
                                            "rep_mnemo_name": "MCRB_ECOLI",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_id": {
                                            "ChEBI:CHEBI:37565": {
                                                "rep_primary_accession": "P15005",
                                                "rep_mnemo_name": "MCRB_ECOLI",
                                                "count": 1
                                            }
                                        }
                                    },
                                    "GO": {
                                        "MCRB_ECOLI": {
                                            "terms": {
                                                "GO:0005524": "ATP binding",
                                                "GO:0016887": "ATP hydrolysis activity"
                                            },
                                            "jaccard_index": 0.1667
                                        }
                                    }
                                }
                            },
                            "330": {
                                "BINDING | Interacts with GTP": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with GTP",
                                        "count": 1,
                                        "annot_amino_acid": "P",
                                        "target_amino_acid": "P"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000255": {
                                            "rep_primary_accession": "P15005",
                                            "rep_mnemo_name": "MCRB_ECOLI",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_id": {
                                            "ChEBI:CHEBI:37565": {
                                                "rep_primary_accession": "P15005",
                                                "rep_mnemo_name": "MCRB_ECOLI",
                                                "count": 1
                                            }
                                        }
                                    },
                                    "GO": {
                                        "MCRB_ECOLI": {
                                            "terms": {
                                                "GO:0005524": "ATP binding",
                                                "GO:0016887": "ATP hydrolysis activity"
                                            },
                                            "jaccard_index": 0.1667
                                        }
                                    }
                                }
                            }
                        },
                        "indices": {
                            "matches": [
                                "329",
                                "330"
                            ],
                            "misses": []
                        }
                    },
                    "conservation_ranges": {
                        "conserved_positions": {
                            "positions": [329, 332, 334, 335],
                            "ranges": [
                                [329, 329],
                                [332, 332],
                                [334, 335],
                            ],
                        }
                    },
                    "conservations": {
                        "positions": {
                            "329": {
                                "conservation": 0.9853,
                                "hit": True
                            },
                            "332": {
                                "conservation": 0.8806,
                                "hit": True
                            },
                            "334": {
                                "conservation": 0.9084,
                                "hit": True
                            },
                            "335": {
                                "conservation": 0.9562,
                                "hit": True
                            }
                        },
                        "indices": {
                            "matches": [
                                "329",
                                "332",
                                "334",
                                "335"
                            ],
                            "misses": []
                        }
                    },
                    "position_conversion": {
                        "target_to_aln": {
                            "329": "14",
                            "330": "15",
                            "332": "17",
                            "334": "19",
                            "335": "20"
                        },
                        "aln_to_target": {
                            "14": "329",
                            "15": "330",
                            "17": "332",
                            "19": "334",
                            "20": "335"
                        }
                    },
                    "annotation_ranges": {
                        "BINDING | Interacts with GTP": {
                            "positions": [
                                329,
                                330
                            ],
                            "ranges": [
                                [
                                    329,
                                    330
                                ]
                            ]
                        }
                    }
                }
            }
        }
    }
}

@pytest.fixture
def transfer_dict_success_all_types_H0YB80():
    return {
    "sequence_id": "tr|H0YB80|H0YB80_HUMAN",
    "domain": {
        "PF00244": {
            "hit_intervals": {
                "1-114": {
                    "sequence": "VFYLKMKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLSEESYKDSTLIMQLLRDNLTLWT",
                    "length": 114,
                    "hit_start": 1,
                    "hit_end": 114,
                    "annotations": {
                        "positions": {
                            "12": {
                                "BINDING | Interacts with peptide": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with peptide",
                                        "count": 1,
                                        "annot_amino_acid": "R",
                                        "target_amino_acid": "R"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:21853016": {
                                            "rep_primary_accession": "Q5CUW0",
                                            "rep_mnemo_name": "Q5CUW0_CRYPI",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_ccd_id": {
                                            "peptide": {
                                                "rep_primary_accession": "Q5CUW0",
                                                "rep_mnemo_name": "Q5CUW0_CRYPI",
                                                "count": 1
                                            }
                                        }
                                    },
                                    "GO": {
                                        "Q5CUW0_CRYPI": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                },
                                "BINDING | Interacts with O-phospho-L-serine": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with O-phospho-L-serine",
                                        "count": 1,
                                        "annot_amino_acid": "R",
                                        "target_amino_acid": "R"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:26551337": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "paired_position": {
                                        "13": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_id": {
                                            "ChEBI:CHEBI:57524": {
                                                "rep_primary_accession": "E2RU97",
                                                "rep_mnemo_name": "1433_GIAIC",
                                                "count": 1
                                            }
                                        }
                                    },
                                    "GO": {
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "13": {
                                "BINDING | Interacts with peptide": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with peptide",
                                        "count": 1,
                                        "annot_amino_acid": "Y",
                                        "target_amino_acid": "Y"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:21853016": {
                                            "rep_primary_accession": "Q5CUW0",
                                            "rep_mnemo_name": "Q5CUW0_CRYPI",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_ccd_id": {
                                            "peptide": {
                                                "rep_primary_accession": "Q5CUW0",
                                                "rep_mnemo_name": "Q5CUW0_CRYPI",
                                                "count": 1
                                            }
                                        }
                                    },
                                    "GO": {
                                        "Q5CUW0_CRYPI": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                },
                                "BINDING | Interacts with O-phospho-L-serine": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with O-phospho-L-serine",
                                        "count": 1,
                                        "annot_amino_acid": "Y",
                                        "target_amino_acid": "Y"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:26551337": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "paired_position": {
                                        "12": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_id": {
                                            "ChEBI:CHEBI:57524": {
                                                "rep_primary_accession": "E2RU97",
                                                "rep_mnemo_name": "1433_GIAIC",
                                                "count": 1
                                            }
                                        }
                                    },
                                    "GO": {
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "57": {
                                "BINDING | Interacts with peptide": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with peptide",
                                        "count": 2,
                                        "annot_amino_acid": "L",
                                        "target_amino_acid": "L"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:21853016": {
                                            "rep_primary_accession": "Q5CUW0",
                                            "rep_mnemo_name": "Q5CUW0_CRYPI",
                                            "count": 1
                                        },
                                        "ECO:0000269|PubMed:26551337": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_ccd_id": {
                                            "peptide": {
                                                "rep_primary_accession": "Q5CUW0",
                                                "rep_mnemo_name": "Q5CUW0_CRYPI",
                                                "count": 2
                                            }
                                        }
                                    },
                                    "GO": {
                                        "Q5CUW0_CRYPI": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        },
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "58": {
                                "BINDING | Interacts with peptide": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with peptide",
                                        "count": 2,
                                        "annot_amino_acid": "N",
                                        "target_amino_acid": "N"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:21853016": {
                                            "rep_primary_accession": "Q5CUW0",
                                            "rep_mnemo_name": "Q5CUW0_CRYPI",
                                            "count": 1
                                        },
                                        "ECO:0000269|PubMed:26551337": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_ccd_id": {
                                            "peptide": {
                                                "rep_primary_accession": "Q5CUW0",
                                                "rep_mnemo_name": "Q5CUW0_CRYPI",
                                                "count": 2
                                            }
                                        }
                                    },
                                    "GO": {
                                        "Q5CUW0_CRYPI": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        },
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "65": {
                                "BINDING | Interacts with peptide": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with peptide",
                                        "count": 1,
                                        "annot_amino_acid": "E",
                                        "target_amino_acid": "E"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:21853016": {
                                            "rep_primary_accession": "Q5CUW0",
                                            "rep_mnemo_name": "Q5CUW0_CRYPI",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_ccd_id": {
                                            "peptide": {
                                                "rep_primary_accession": "Q5CUW0",
                                                "rep_mnemo_name": "Q5CUW0_CRYPI",
                                                "count": 1
                                            }
                                        }
                                    },
                                    "GO": {
                                        "Q5CUW0_CRYPI": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "109": {
                                "BINDING | Interacts with peptide": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with peptide",
                                        "count": 2,
                                        "annot_amino_acid": "N",
                                        "target_amino_acid": "N"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:21853016": {
                                            "rep_primary_accession": "Q5CUW0",
                                            "rep_mnemo_name": "Q5CUW0_CRYPI",
                                            "count": 1
                                        },
                                        "ECO:0000269|PubMed:26551337": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_ccd_id": {
                                            "peptide": {
                                                "rep_primary_accession": "Q5CUW0",
                                                "rep_mnemo_name": "Q5CUW0_CRYPI",
                                                "count": 2
                                            }
                                        }
                                    },
                                    "GO": {
                                        "Q5CUW0_CRYPI": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        },
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "113": {
                                "BINDING | Interacts with peptide": {
                                    "essentials": {
                                        "type": "BINDING",
                                        "description": "Interacts with peptide",
                                        "count": 1,
                                        "annot_amino_acid": "W",
                                        "target_amino_acid": "W"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:21853016": {
                                            "rep_primary_accession": "Q5CUW0",
                                            "rep_mnemo_name": "Q5CUW0_CRYPI",
                                            "count": 1
                                        }
                                    },
                                    "additional_keys": {
                                        "ligand_ccd_id": {
                                            "peptide": {
                                                "rep_primary_accession": "Q5CUW0",
                                                "rep_mnemo_name": "Q5CUW0_CRYPI",
                                                "count": 1
                                            }
                                        }
                                    },
                                    "GO": {
                                        "Q5CUW0_CRYPI": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "61": {
                                "MUTAGEN | V->D: Loss of binding to difopein.": {
                                    "essentials": {
                                        "type": "MUTAGEN",
                                        "description": "V->D: Loss of binding to difopein.",
                                        "count": 1,
                                        "annot_amino_acid": "V",
                                        "target_amino_acid": "V"
                                    },
                                    "hit": True,
                                    "evidence": {
                                        "ECO:0000269|PubMed:19733174": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "GO": {
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "78": {
                                "MUTAGEN | R->K: Increased oligomerization.": {
                                    "essentials": {
                                        "type": "MUTAGEN",
                                        "description": "R->K: Increased oligomerization.",
                                        "count": 1,
                                        "annot_amino_acid": "R",
                                        "target_amino_acid": "K"
                                    },
                                    "hit": False,
                                    "evidence": {
                                        "ECO:0000269|PubMed:24658679": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "GO": {
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "86": {
                                "MUTAGEN | T->A: Slightly decreased oligomerization.": {
                                    "essentials": {
                                        "type": "MUTAGEN",
                                        "description": "T->A: Slightly decreased oligomerization.",
                                        "count": 1,
                                        "annot_amino_acid": "T",
                                        "target_amino_acid": "A"
                                    },
                                    "hit": False,
                                    "evidence": {
                                        "ECO:0000269|PubMed:24658679": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "GO": {
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            },
                            "92": {
                                "MOD_RES | Phosphothreonine": {
                                    "essentials": {
                                        "type": "MOD_RES",
                                        "description": "Phosphothreonine",
                                        "count": 1,
                                        "annot_amino_acid": "T",
                                        "target_amino_acid": "S"
                                    },
                                    "hit": False,
                                    "evidence": {
                                        "ECO:0000269|PubMed:16368691, ECO:0000269|PubMed:19733174, ECO:0000269|PubMed:24147113, ECO:0000269|PubMed:24658679": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "GO": {
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                },
                                "MUTAGEN | T->A: Loss of phosphorylation by a protein kinase. No effect on subcellular localization. Dramatic decrease in the number of encysting parasites and cysts, but a large increase in the number of trophozoites. In encysting cells of 12 hours, significantly slower cyst conversion rate compared to the wild-type. No effect on binding to difopein. Decreased binding to a number of synthetic phosphopeptides.": {
                                    "essentials": {
                                        "type": "MUTAGEN",
                                        "description": "T->A: Loss of phosphorylation by a protein kinase. No effect on subcellular localization. Dramatic decrease in the number of encysting parasites and cysts, but a large increase in the number of trophozoites. In encysting cells of 12 hours, significantly slower cyst conversion rate compared to the wild-type. No effect on binding to difopein. Decreased binding to a number of synthetic phosphopeptides.",
                                        "count": 1,
                                        "annot_amino_acid": "T",
                                        "target_amino_acid": "S"
                                    },
                                    "hit": False,
                                    "evidence": {
                                        "ECO:0000269|PubMed:16368691, ECO:0000269|PubMed:19733174": {
                                            "rep_primary_accession": "E2RU97",
                                            "rep_mnemo_name": "1433_GIAIC",
                                            "count": 1
                                        }
                                    },
                                    "GO": {
                                        "1433_GIAIC": {
                                            "terms": {},
                                            "jaccard_index": 0.0
                                        }
                                    }
                                }
                            }
                        },
                        "indices": {
                            "matches": [
                                "109",
                                "113",
                                "12",
                                "13",
                                "57",
                                "58",
                                "61",
                                "65"
                            ],
                            "misses": [
                                "78",
                                "86",
                                "92"
                            ]
                        }
                    },
                    "conservations": {
                        "positions": {
                            "1": {
                                "conservation": 0.8855,
                                "hit": True
                            },
                            "2": {
                                "conservation": 0.9213,
                                "hit": True
                            },
                            "3": {
                                "conservation": 0.8509,
                                "hit": True
                            },
                            "5": {
                                "conservation": 0.943,
                                "hit": True
                            },
                            "6": {
                                "conservation": 0.8856,
                                "hit": True
                            },
                            "7": {
                                "conservation": 0.8162,
                                "hit": True
                            },
                            "8": {
                                "conservation": 0.8751,
                                "hit": True
                            },
                            "9": {
                                "conservation": 0.986,
                                "hit": True
                            },
                            "10": {
                                "conservation": 0.9316,
                                "hit": True
                            },
                            "12": {
                                "conservation": 0.9776,
                                "hit": True
                            },
                            "13": {
                                "conservation": 0.9815,
                                "hit": True
                            },
                            "15": {
                                "conservation": 0.8356,
                                "hit": True
                            },
                            "16": {
                                "conservation": 0.9682,
                                "hit": True
                            },
                            "34": {
                                "conservation": 0.9736,
                                "hit": True

                            },
                            "37": {
                                "conservation": 0.9657,
                                "hit": True

                            },
                            "41": {
                                "conservation": 0.812,
                                "hit": False
                            },
                            "48": {
                                "conservation": 0.8392,
                                "hit": True

                            },
                            "50": {
                                "conservation": 0.9312,
                                "hit": True

                            },
                            "52": {
                                "conservation": 0.9043,
                                "hit": True

                            },
                            "53": {
                                "conservation": 0.9743,
                                "hit": True

                            },
                            "54": {
                                "conservation": 0.9191,
                                "hit": True

                            },
                            "55": {
                                "conservation": 0.9228,
                                "hit": True

                            },
                            "56": {
                                "conservation": 0.8745,
                                "hit": True

                            },
                            "57": {
                                "conservation": 0.978,
                                "hit": True

                            },
                            "58": {
                                "conservation": 0.9648,
                                "hit": True
                            },
                            "60": {
                                "conservation": 0.9051,
                                "hit": True
                            },
                            "61": {
                                "conservation": 0.9377,
                                "hit": True
                            },
                            "62": {
                                "conservation": 0.9346,
                                "hit": True
                            },
                            "64": {
                                "conservation": 0.8826,
                                "hit": True
                            },
                            "65": {
                                "conservation": 0.9215,
                                "hit": True
                            },
                            "66": {
                                "conservation": 0.8339,
                                "hit": True
                            },
                            "73": {
                                "conservation": 0.9722,
                                "hit": True
                            },
                            "74": {
                                "conservation": 0.8229,
                                "hit": True
                            },
                            "77": {
                                "conservation": 0.9375,
                                "hit": True
                            },
                            "80": {
                                "conservation": 0.8853,
                                "hit": True
                            },
                            "81": {
                                "conservation": 0.87,
                                "hit": True
                            },
                            "82": {
                                "conservation": 0.8009,
                                "hit": True
                            },
                            "84": {
                                "conservation": 0.909,
                                "hit": True
                            },
                            "88": {
                                "conservation": 0.8463,
                                "hit": True
                            },
                            "89": {
                                "conservation": 0.8356,
                                "hit": True
                            },
                            "93": {
                                "conservation": 0.8168,
                                "hit": True
                            },
                            "96": {
                                "conservation": 0.8642,
                                "hit": True
                            },
                            "97": {
                                "conservation": 0.8046,
                                "hit": True
                            },
                            "98": {
                                "conservation": 0.8844,
                                "hit": True
                            },
                            "99": {
                                "conservation": 0.8355,
                                "hit": True
                            },
                            "100": {
                                "conservation": 0.832,
                                "hit": True
                            },
                            "101": {
                                "conservation": 0.8149,
                                "hit": True
                            },
                            "102": {
                                "conservation": 0.9163,
                                "hit": True
                            },
                            "103": {
                                "conservation": 0.8082,
                                "hit": True
                            },
                            "104": {
                                "conservation": 0.8983,
                                "hit": True
                            },
                            "105": {
                                "conservation": 0.9404,
                                "hit": True
                            },
                            "106": {
                                "conservation": 0.9178,
                                "hit": True
                            },
                            "107": {
                                "conservation": 0.8887,
                                "hit": True
                            },
                            "108": {
                                "conservation": 0.9231,
                                "hit": True
                            },
                            "109": {
                                "conservation": 0.9685,
                                "hit": True
                            },
                            "110": {
                                "conservation": 0.9009,
                                "hit": True
                            },
                            "111": {
                                "conservation": 0.8181,
                                "hit": True
                            },
                            "112": {
                                "conservation": 0.904,
                                "hit": True
                            },
                            "113": {
                                "conservation": 0.9993,
                                "hit": True
                            },
                            "114": {
                                "conservation": 0.8651,
                                "hit": True
                            }
                        },
                        "indices": {
                            "matches": [
                                "1",
                                "10",
                                "100",
                                "101",
                                "102",
                                "103",
                                "104",
                                "105",
                                "106",
                                "107",
                                "108",
                                "109",
                                "110",
                                "111",
                                "112",
                                "113",
                                "114",
                                "12",
                                "13",
                                "15",
                                "16",
                                "2",
                                "3",
                                "34",
                                "37",
                                "48",
                                "5",
                                "50",
                                "52",
                                "53",
                                "54",
                                "55",
                                "56",
                                "57",
                                "58",
                                "6",
                                "60",
                                "61",
                                "62",
                                "64",
                                "65",
                                "66",
                                "7",
                                "73",
                                "74",
                                "77",
                                "8",
                                "80",
                                "81",
                                "82",
                                "84",
                                "88",
                                "89",
                                "9",
                                "93",
                                "96",
                                "97",
                                "98",
                                "99"
                            ],
                            "misses": [
                                "41"
                            ]
                        }
                    },
                    "position_conversion": {
                        "target_to_aln": {
                            "12": "164",
                            "13": "165",
                            "57": "227",
                            "58": "228",
                            "65": "235",
                            "109": "294",
                            "113": "298",
                            "61": "231",
                            "78": "251",
                            "86": "262",
                            "92": "274",
                            "1": "151",
                            "2": "153",
                            "3": "154",
                            "5": "157",
                            "6": "158",
                            "7": "159",
                            "8": "160",
                            "9": "161",
                            "10": "162",
                            "15": "167",
                            "16": "168",
                            "34": "199",
                            "37": "202",
                            "41": "206",
                            "48": "218",
                            "50": "220",
                            "52": "222",
                            "53": "223",
                            "54": "224",
                            "55": "225",
                            "56": "226",
                            "60": "230",
                            "62": "232",
                            "64": "234",
                            "66": "236",
                            "73": "246",
                            "74": "247",
                            "77": "250",
                            "80": "253",
                            "81": "257",
                            "82": "258",
                            "84": "260",
                            "88": "266",
                            "89": "270",
                            "93": "276",
                            "96": "280",
                            "97": "281",
                            "98": "282",
                            "99": "283",
                            "100": "284",
                            "101": "285",
                            "102": "287",
                            "103": "288",
                            "104": "289",
                            "105": "290",
                            "106": "291",
                            "107": "292",
                            "108": "293",
                            "110": "295",
                            "111": "296",
                            "112": "297",
                            "114": "299"
                        },
                        "aln_to_target": {
                            "164": "12",
                            "165": "13",
                            "227": "57",
                            "228": "58",
                            "235": "65",
                            "294": "109",
                            "298": "113",
                            "231": "61",
                            "251": "78",
                            "262": "86",
                            "274": "92",
                            "151": "1",
                            "153": "2",
                            "154": "3",
                            "157": "5",
                            "158": "6",
                            "159": "7",
                            "160": "8",
                            "161": "9",
                            "162": "10",
                            "167": "15",
                            "168": "16",
                            "199": "34",
                            "202": "37",
                            "206": "41",
                            "218": "48",
                            "220": "50",
                            "222": "52",
                            "223": "53",
                            "224": "54",
                            "225": "55",
                            "226": "56",
                            "230": "60",
                            "232": "62",
                            "234": "64",
                            "236": "66",
                            "246": "73",
                            "247": "74",
                            "250": "77",
                            "253": "80",
                            "257": "81",
                            "258": "82",
                            "260": "84",
                            "266": "88",
                            "270": "89",
                            "276": "93",
                            "280": "96",
                            "281": "97",
                            "282": "98",
                            "283": "99",
                            "284": "100",
                            "285": "101",
                            "287": "102",
                            "288": "103",
                            "289": "104",
                            "290": "105",
                            "291": "106",
                            "292": "107",
                            "293": "108",
                            "295": "110",
                            "296": "111",
                            "297": "112",
                            "299": "114"
                        }
                    },
                    "annotation_ranges": {
                        "BINDING | Interacts with peptide": {
                            "positions": [
                                12,
                                13,
                                57,
                                58,
                                65,
                                109,
                                113
                            ],
                            "ranges": [
                                [
                                    12,
                                    13
                                ],
                                [
                                    57,
                                    58
                                ],
                                [
                                    65,
                                    65
                                ],
                                [
                                    109,
                                    109
                                ],
                                [
                                    113,
                                    113
                                ]
                            ]
                        },
                        "BINDING | Interacts with O-phospho-L-serine": {
                            "positions": [
                                12,
                                13
                            ],
                            "ranges": [
                                [
                                    12,
                                    13
                                ]
                            ]
                        },
                        "MUTAGEN | V->D: Loss of binding to difopein.": {
                            "positions": [
                                61
                            ],
                            "ranges": [
                                [
                                    61,
                                    61
                                ]
                            ]
                        },
                        "MUTAGEN | R->K: Increased oligomerization.": {
                            "positions": [
                                78
                            ],
                            "ranges": [
                                [
                                    78,
                                    78
                                ]
                            ]
                        },
                        "MUTAGEN | T->A: Slightly decreased oligomerization.": {
                            "positions": [
                                86
                            ],
                            "ranges": [
                                [
                                    86,
                                    86
                                ]
                            ]
                        },
                        "MOD_RES | Phosphothreonine": {
                            "positions": [
                                92
                            ],
                            "ranges": [
                                [
                                    92,
                                    92
                                ]
                            ]
                        },
                        "MUTAGEN | T->A: Loss of phosphorylation by a protein kinase. No effect on subcellular localization. Dramatic decrease in the number of encysting parasites and cysts, but a large increase in the number of trophozoites. In encysting cells of 12 hours, significantly slower cyst conversion rate compared to the wild-type. No effect on binding to difopein. Decreased binding to a number of synthetic phosphopeptides.": {
                            "positions": [
                                92
                            ],
                            "ranges": [
                                [
                                    92,
                                    92
                                ]
                            ]
                        }
                    },
                    "conservation_ranges": {
                        "conserved_positions": {
                            "positions": [
                                1,
                                2,
                                3,
                                5,
                                6,
                                7,
                                8,
                                9,
                                10,
                                12,
                                13,
                                15,
                                16,
                                34,
                                37,
                                41,
                                48,
                                50,
                                52,
                                53,
                                54,
                                55,
                                56,
                                57,
                                58,
                                60,
                                61,
                                62,
                                64,
                                65,
                                66,
                                73,
                                74,
                                77,
                                80,
                                81,
                                82,
                                84,
                                88,
                                89,
                                93,
                                96,
                                97,
                                98,
                                99,
                                100,
                                101,
                                102,
                                103,
                                104,
                                105,
                                106,
                                107,
                                108,
                                109,
                                110,
                                111,
                                112,
                                113,
                                114
                            ],
                            "ranges": [
                                [
                                    1,
                                    3
                                ],
                                [
                                    5,
                                    10
                                ],
                                [
                                    12,
                                    13
                                ],
                                [
                                    15,
                                    16
                                ],
                                [
                                    34,
                                    34
                                ],
                                [
                                    37,
                                    37
                                ],
                                [
                                    41,
                                    41
                                ],
                                [
                                    48,
                                    48
                                ],
                                [
                                    50,
                                    50
                                ],
                                [
                                    52,
                                    58
                                ],
                                [
                                    60,
                                    62
                                ],
                                [
                                    64,
                                    66
                                ],
                                [
                                    73,
                                    74
                                ],
                                [
                                    77,
                                    77
                                ],
                                [
                                    80,
                                    82
                                ],
                                [
                                    84,
                                    84
                                ],
                                [
                                    88,
                                    89
                                ],
                                [
                                    93,
                                    93
                                ],
                                [
                                    96,
                                    114
                                ]
                            ]
                        }
                    }
                }
            }
        }
    }
}

## Added in Process_annotation
@pytest.fixture
def get_annotation_dict():
    """Helper fixture that provides individual annotation dictionaries"""
    def get_annotation(fixture, pos, index=0):
        """Extract specific annotation dict from fixture
        Args:
            fixture: The source fixture (for now only
            annotations_content_binding_fixture_Q9NU22_PF07728 or annotations_content_disulfid_fixture_Q9NU22_PF07728)
            pos: Position key (e.g. "201")
            index: Index in annotation list (default 0)
        """
        return fixture["MCRB_ECOLI"][pos][index]
    return get_annotation

@pytest.fixture
def annotation_dict_205_Q9NU22(annotations_content_disulfid_fixture_Q9NU22_PF07728, get_annotation_dict):
    return get_annotation_dict(annotations_content_disulfid_fixture_Q9NU22_PF07728, "205")

@pytest.fixture
def annotation_dict_246_Q9NU22(annotations_content_disulfid_fixture_Q9NU22_PF07728, get_annotation_dict):
    return get_annotation_dict(annotations_content_disulfid_fixture_Q9NU22_PF07728, "246")

@pytest.fixture
def mock_make_anno_total_disulfid_return_205_P15005(annotation_dict_205_Q9NU22):
    return {
        "annotation": annotation_dict_205_Q9NU22,
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-246); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-246); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "333",
            "annot_position": "205",
            "index_position": "18",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
            "paired_target_position": "373"
        },
        "paired_annot_pos_str": "246"
    }

@pytest.fixture
def mock_make_anno_total_disulfid_return_246_P15005(annotation_dict_205_Q9NU22):
   return {
        "annotation": annotation_dict_205_Q9NU22,
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-205); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "373",
            "annot_position": "246",
            "index_position": "72",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
            "paired_target_position": "333"
        },
        "paired_annot_pos_str": "205"
    }

@pytest.fixture
def mock_make_anno_total_disulfid_return_205_Q9NU22():
    return {
        "annotation": {
                    "type": "DISULFID",
                    "description": "Intrachain (with C-246); in linked form",
                    "evidence": "ECO:0000250", # Changed from 269 to 250, testing different ECOs
                    "entry": "P00750", # Random Mock
                    "aminoacid": "C",
                    "paired_position": "246"
                },
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-246); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-246); in linked form",
            "count": 1,
            "evidence": "ECO:0000250", # Changed from 269 to 250, testing different ECOs
            "target_position": "333",
            "annot_position": "205",
            "index_position": "18",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
            "paired_target_position": "373"
        },
        "paired_annot_pos_str": "246"
    }

@pytest.fixture
def mock_make_anno_total_disulfid_return_246_Q9NU22():
    return {
        "annotation": {
                    "type": "DISULFID",
                    "description": "Intrachain (with C-205); in linked form",
                    "evidence": "ECO:0000250", # Changed from 269 to 250, testing different ECOs
                    "entry": "P00750", # Random Mock
                    "aminoacid": "C",
                    "paired_position": "205"
                },
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-205); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "count": 1,
            "evidence": "ECO:0000250", # Changed from 269 to 250, testing different ECOs
            "target_position": "373",
            "annot_position": "246",
            "index_position": "72",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
            "paired_target_position": "333"
        },
        "paired_annot_position": "205"
    }


@pytest.fixture
def mock_map_filter_disulfid_return(annotation_dict_246_Q9NU22):
    return (True, {
        "annotation": annotation_dict_246_Q9NU22,
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-205); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "373",
            "annot_position": "246",
            "index_position": "72",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
            "paired_target_position": "333"
        },
        "paired_annot_position": "205"
    })

# add_transfer_dict

@pytest.fixture
def anno_total_disulfid_MCRB_ECOLI_Q9NU22_205():
    """Fixture for mock annotation data."""
    return {
        "type": "DISULFID",
        "description": "Intrachain (with C-246); in linked form",
        "count": 1,
        "evidence": "ECO:0000269|PubMed:12345678",
        "target_position": "333",
        "annot_position": "205",
        "paired_target_position": "373",
        "index_position": "18",
        "annot_amino_acid": "C",
        "target_amino_acid": "C"
    }

@pytest.fixture
def anno_total_disulfid_MCRB_ECOLI_Q9NU22_246():
    """Fixture for mock paired annotation data."""
    return {
        "type": "DISULFID",
        "description": "Intrachain (with C-205); in linked form",
        "count": 1,
        "evidence": "ECO:0000269|PubMed:12345678",
        "target_position": "373",
        "annot_position": "246",
        "paired_target_position": "333",
        "index_position": "72",
        "annot_amino_acid": "C",
        "target_amino_acid": "C"
    }

@pytest.fixture
def minimal_hmmalign_lines_fixture_Q9NU22():
    return """# STOCKHOLM 1.0
Q7US48_RHOBA/138-284                                    .........GLMMVGEPGTAKSMLGE.LLAVAI.............sGTGSLAVQGTA.GT..TEDQIKYGWNYAmll......drgPVHEALVPSPVMTAMR..................DGKVVRFEEITRCL.P........E....VQD....A......LIS.ILSER.R...MMIPEMDQAGDDNANS.vFASP............................gFSIIATANLRD.........RGVSEMSAALKRRF............................
MCRB_ECOLI/196-350                                      .........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYCCN............GVGFRRKDGIFYNFCQqak...........eqpeKKYIFIIDEINRANlS........K....VFG....E......VMM.LMEHD.K...RGENWSVPLTYSENDE.eRFYVp..........................enVYIIGLMNTAD.........RSLAVVDYALRRRF............................
sp|Q9NU22|MDN1_HUMANtarget//325-451                     .........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................-FQFFAT-----.........--------------rrllscgg...................
//
"""

@pytest.fixture
def minimal_hmmalign_lines_fixture_H0YB80():
    return """# STOCKHOLM 1.0
C1E2K1_MICCC/19-248                         KAKLAEQ.AERYDEMMEYMSKVAESgsggaDAEL.SVEERNLLSVAYKNVIGAR.RASWRIIS........S.....I.ETKEESKYt...kESDVVR......LIQKYKSNVEKELSDICDRILKLLRDHLE.ET..S...SAG...ESKV.FYK.KMKGDYYRYLAE.FK.......GG..E...ARKNAAEETLLAYKEAENIA.....SNELAPTHPIRLGLALNFSVFYYEILNAPER...ACDMAKKA...FDEAIAE..L...DT.LG.EES.YKDSTL.IMQLLRDNLTLWT
14338_ARATH/13-238                          MAKLAEQ.AERYEEMVQFMEQLVSGa..tpAGEL.TVEERNLLSVAYKNVIGSL.RAAWRIVS........S.....I.EQKEESR-....kNEEHVS......LVKDYRSKVETELSSICSGILRLLDSHLI.PS..A...TAS...ESKV.FYL.KMKGDYHRYLAE.FK.......SG..D...ERKTAAEDTMIAYKAAQDVA.....VADLAPTHPIRLGLALNFSVFYYEILNSSEK...ACSMAKQA...FEEAIAE..L...DT.LG.EES.YKDSTL.IMQLLRDNLTLWT
1433_GIAIC/13-236                           MAQLNEN.AERYDEMVETMRKISGM.....EGEL.SDKERNLLSVAYKNVIGPR.RAAWRIVS........S.....I.EAKEKGRQk...pNAKRIE......QIRVYRQKIEKELSDICNDILKLLQEQFV.PR..S...TNA...DAKV.FYY.KMQGDYYRYLAE.YS.......SG..E...DKEKIAGSALNAYNSAFEIS.....-QQLPPTHPIRLGLALNFSVFYYEILASPDR...ACELARKA...FDAAITD..L...DK.LT.EES.YKDSTL.IMQLLRDNLNLWV
Q5CUW0_CRYPI/33-257                         MAKLAEQ.AERYDEMAKYMKDVVEAr..qeSEEL.TVEERNLLSVAYKNAVGSR.RSSWRIIS........S.....V.EQKEHSRN.....AEDASK......MCGKYRSKVEAELTDICNDILTMLDKHLI.PT..A...TSP...DSKV.FYF.KMKGDYHRYISE.FS.......TG..D...SKQSSAEDALKAYKDATVVA.....-KDLEPTHPIRLGLALNFSVFHYEILNEPRA...AIDMAKEA...FEMAIEQ..L...DK.LS.EDC.YKDSTL.IMQLLRDNLTLWT
tr|H0YB80|H0YB80_HUMANtarget//1-114         -------.-----------------.....----.-------------------.--------........-.....-.--------.....------......-----------------------------.--..-...---...---V.FYL.KMKGDYYRYLAE.VA.......AG..D...DKKGIVDQSQQAYQEAFEIS.....KKEMQPTHPIRLGLALNFSVFYYEILNSPEK...ACSLAKTA...FDEAIAE..L...DT.LS.EES.YKDSTL.IMQLLRDNLTLWT
#=GR tr|H0YB80|H0YB80_HUMANtarget//1-114 PP .......................................................................................................................................................8.***.************.**.......**..*...********************.....*******************************...********...*******..*...**.**.***.******.************8
//
"""

### Fixtures added in cleanup_improve_transfer_dict and helper functions

def mock_populate_conservation(transfer_dict, source_dict):
    """Helper to populate conservation data from source"""
    target_path = transfer_dict["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
    source_path = source_dict["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
    target_path["conservations"] = copy.deepcopy(source_path["conservations"])
    target_path["position_conversion"] = copy.deepcopy(source_path["position_conversion"])
    target_path["conservation_ranges"] = copy.deepcopy(source_path["conservation_ranges"])

def mock_populate_go_data(transfer_dict, source_dict):
    """Helper to populate GO data from source"""
    target_path = transfer_dict["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
    source_path = source_dict["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
    for pos in ["333", "373"]:
        anno_key = list(target_path["annotations"]["positions"][pos].keys())[0]
        target_path["annotations"]["positions"][pos][anno_key]["GO"] = \
            copy.deepcopy(source_path["annotations"]["positions"][pos][anno_key]["GO"])

@pytest.fixture
def target_gos_Q9NU22_MCRB_ECOLI():
    return {"GO:0005524", "GO:0016887"}

def format_tsv_line(*args):
    """Helper to format TSV lines consistently"""
    return "\t".join(str(arg).strip() for arg in args)

@pytest.fixture
def mock_iprscan_data():
    """Real test data from iprscan_content_Q9NU22_PF07728"""
    return [
        ["sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF07728",
         "AAA domain", "325", "448", "2.7E-13", "T", "09-01-2025", "IPR011704", "ATPase domain",
         "GO:0005524(InterPro)|GO:0016887(InterPro)", "-"],
        ["sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Gene3D",
         "G3DSA:3.40.50.300", "-", "309", "486", "2.7E-28", "T", "09-01-2025", "IPR027417",
         "P-loop domain", "GO:0042623(InterPro)", "-"],
        ["sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "SMART",
         "SM00382", "AAA_5", "321", "634", "19.0", "T", "09-01-2025", "IPR003593",
         "ATPase domain", "GO:0016887(InterPro)", "-"]
    ]

@pytest.fixture
def iprscan_df_Q9NU22_PF07728(mock_iprscan_data):
    """Creates DataFrame with proper columns"""
    return pd.DataFrame(mock_iprscan_data, columns=[
        "Protein Accession", "Sequence MD5 digest", "Sequence Length", "Analysis",
        "Signature Accession", "Signature Description", "Start Location", "Stop Location",
        "Score", "Status", "Date", "InterPro Accession", "InterPro Description",
        "GO Annotations", "Pathways Annotations"
    ])

@pytest.fixture
def mock_empty_df():
    """Creates empty DataFrame with correct columns"""
    return pd.DataFrame(columns=[
        "Protein Accession", "Sequence MD5 digest", "Sequence Length", "Analysis",
        "Signature Accession", "Signature Description", "Start Location", "Stop Location",
        "Score", "Status", "Date", "InterPro Accession", "InterPro Description",
        "GO Annotations", "Pathways Annotations"
    ])

@pytest.fixture
def iprscan_content_Q9NU22():
    """Create mock InterProScan TSV content string for test."""
    lines = [
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Gene3D", "G3DSA:3.40.50.410", "von Willebrand factor, type A domain", "5379", "5545", "3.8E-8", "T", "09-01-2025", "IPR036465", "von Willebrand factor A-like domain superfamily", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF17865", "Midasin AAA lid domain", "915", "1015", "2.7E-29", "T", "09-01-2025", "IPR041190", "Midasin AAA lid domain 5", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Gene3D", "G3DSA:3.40.50.300", "-", "309", "486", "2.7E-28", "T", "09-01-2025", "IPR027417", "P-loop containing nucleoside triphosphate hydrolase", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Gene3D", "G3DSA:3.40.50.300", "-", "1732", "1903", "6.5E-39", "T", "09-01-2025", "IPR027417", "P-loop containing nucleoside triphosphate hydrolase", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Gene3D", "G3DSA:3.40.50.300", "-", "1063", "1233", "3.0E-36", "T", "09-01-2025", "IPR027417", "P-loop containing nucleoside triphosphate hydrolase", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "SUPERFAMILY", "SSF52540", "P-loop containing nucleoside triphosphate hydrolases", "1048", "1316", "7.47E-43", "T", "09-01-2025", "IPR027417", "P-loop containing nucleoside triphosphate hydrolase", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF17867", "Midasin AAA lid domain", "1900", "2000", "4.2E-18", "T", "09-01-2025", "IPR040848", "Midasin, AAA lid domain 7", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF17867", "Midasin AAA lid domain", "1227", "1328", "6.5E-17", "T", "09-01-2025", "IPR040848", "Midasin, AAA lid domain 7", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF17867", "Midasin AAA lid domain", "480", "604", "2.2E-13", "T", "09-01-2025", "IPR040848", "Midasin, AAA lid domain 7", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF07728", "AAA domain (dynein-related subfamily)", "325", "448", "2.7E-13", "T", "09-01-2025", "IPR011704", "ATPase, dynein-related, AAA domain", "GO:0005524(InterPro)|GO:0016887(InterPro)", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF07728", "AAA domain (dynein-related subfamily)", "1749", "1888", "7.0E-12", "T", "09-01-2025", "IPR011704", "ATPase, dynein-related, AAA domain", "GO:0005524(InterPro)|GO:0016887(InterPro)", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF07728", "AAA domain (dynein-related subfamily)", "672", "756", "9.7E-4", "T", "09-01-2025", "IPR011704", "ATPase, dynein-related, AAA domain", "GO:0005524(InterPro)|GO:0016887(InterPro)", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "SMART", "SM00382", "AAA_5", "321", "634", "19.0", "T", "09-01-2025", "IPR003593", "AAA+ ATPase domain", "GO:0016887(InterPro)", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "PANTHER", "PTHR48103", "MIDASIN-RELATED", "29", "5594", "0.0", "T", "09-01-2025", "-", "-", "GO:0000027(PANTHER)|GO:0000055(PANTHER)|GO:0005634(PANTHER)|GO:0030687(PANTHER)", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "SUPERFAMILY", "SSF52540", "P-loop containing nucleoside triphosphate hydrolases", "1362", "1616", "1.54E-32", "T", "09-01-2025", "IPR027417", "P-loop containing nucleoside triphosphate hydrolase", "-", "-")
    ]
    return "\n".join(lines)

@pytest.fixture
def iprscan_content_H0YB80():
    """Create mock InterProScan TSV content string for test."""
    lines = [
        format_tsv_line("tr|H0YB80|H0YB80_HUMAN", "1a8174c6dd031b64374951bc37250399", "130", "Pfam", "PF00244", "14-3-3 protein", "1", "114", "8.4E-58", "T", "25-12-2024", "IPR023410", "14-3-3 domain", "-", "-"),
        format_tsv_line("tr|H0YB80|H0YB80_HUMAN", "1a8174c6dd031b64374951bc37250399", "130", "PANTHER", "PTHR18860", "14-3-3 PROTEIN", "1", "118", "1.6E-62", "T", "25-12-2024", "IPR000308", "14-3-3 protein", "-", "-"),
        format_tsv_line("tr|H0YB80|H0YB80_HUMAN", "1a8174c6dd031b64374951bc37250399", "130", "SUPERFAMILY", "SSF48445", "14-3-3 protein", "1", "121", "1.44E-59", "T", "25-12-2024", "IPR036815", "14-3-3 domain superfamily", "-", "-"),
        format_tsv_line("tr|H0YB80|H0YB80_HUMAN", "1a8174c6dd031b64374951bc37250399", "130", "Gene3D", "G3DSA:1.20.190.20", "-", "1", "126", "1.5E-71", "T", "25-12-2024", "IPR036815", "14-3-3 domain superfamily", "-", "-"),
        format_tsv_line("tr|H0YB80|H0YB80_HUMAN", "1a8174c6dd031b64374951bc37250399", "130", "SMART", "SM00101", "1433_4", "1", "127", "9.4E-36", "T", "25-12-2024", "IPR023410", "14-3-3 domain", "-", "-")
    ]
    return "\n".join(lines)

@pytest.fixture(scope="session")
def go_terms_dir_mock(tmp_path_factory):
    """Creates a mock directory with both full GO and no-GO data files"""
    base_dir = tmp_path_factory.mktemp("mock_data")

    # Set up full GO terms file
    target_dir = base_dir / "target_name"
    target_dir.mkdir(parents=True)
    iprscan_file = target_dir / "iprscan.tsv"

    lines = [
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Gene3D", "G3DSA:3.40.50.410", "von Willebrand factor, type A domain", "5379", "5545", "3.8E-8", "T", "09-01-2025", "IPR036465", "von Willebrand factor A-like domain superfamily", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF17865", "Midasin AAA lid domain", "915", "1015", "2.7E-29", "T", "09-01-2025", "IPR041190", "Midasin AAA lid domain 5", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Gene3D", "G3DSA:3.40.50.300", "-", "309", "486", "2.7E-28", "T", "09-01-2025", "IPR027417", "P-loop containing nucleoside triphosphate hydrolase", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF07728", "AAA domain (dynein-related subfamily)", "325", "448", "2.7E-13", "T", "09-01-2025", "IPR011704", "ATPase, dynein-related, AAA domain", "GO:0005524(InterPro)|GO:0016887(InterPro)", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "SMART", "SM00382", "AAA_5", "321", "634", "19.0", "T", "09-01-2025", "IPR003593", "AAA+ ATPase domain", "GO:0016887(InterPro)", "-")
    ]
    iprscan_file.write_text("\n".join(lines))

    # Set up no-GO terms file
    target_dir_no_go = base_dir / "target_name_no_go"
    target_dir_no_go.mkdir(parents=True)
    iprscan_file_no_go = target_dir_no_go / "iprscan.tsv"

    lines_no_go = [
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Gene3D", "G3DSA:3.40.50.410", "von Willebrand factor, type A domain", "5379", "5545", "3.8E-8", "T", "09-01-2025", "IPR036465", "von Willebrand factor A-like domain superfamily", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Pfam", "PF17865", "Midasin AAA lid domain", "915", "1015", "2.7E-29", "T", "09-01-2025", "IPR041190", "Midasin AAA lid domain 5", "-", "-"),
        format_tsv_line("sp|Q9NU22|MDN1_HUMAN", "7b57368205d60b6e7b538d2181ad7f2b", "5596", "Gene3D", "G3DSA:3.40.50.300", "-", "309", "486", "2.7E-28", "T", "09-01-2025", "IPR027417", "P-loop containing nucleoside triphosphate hydrolase", "-", "-")
    ]
    iprscan_file_no_go.write_text("\n".join(lines_no_go))

    return str(base_dir)

### Populate Conservation & Populate Go Data For Annotations

@pytest.fixture
def transfer_dict_populated_disulfid_inside_cleanup_Q9NU22(transfer_dict_populated_disulfid_Q9NU22):
    """Transfer dict with DOMAIN replaced by PF07728."""
    cleanup_transfer_dict_ver = copy.deepcopy(transfer_dict_populated_disulfid_Q9NU22)
    cleanup_transfer_dict_ver["PF07728"] = cleanup_transfer_dict_ver.pop("DOMAIN")
    return cleanup_transfer_dict_ver

@pytest.fixture
def transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22(transfer_dict_populated_disulfid_post_conservation_Q9NU22):
    """Creates version of post-processed dict with PF07728 as top key."""
    post_processed_cleanup = copy.deepcopy(transfer_dict_populated_disulfid_post_conservation_Q9NU22)
    content = post_processed_cleanup["domain"]["PF07728"]
    return {"PF07728": content}

@pytest.fixture
def annotations_content_Q9NU22_PF07728():
    """Annotations dictionary with GO terms info for reference annotated sequence in PF07728 annotations JSON."""
    return {
        "MCRB_ECOLI": {
            "0": {  # go_terms_annot_key
                "Molecular Function": {
                    "GO:0005524": "ATP binding",
                    "GO:0016887": "ATP hydrolysis activity",
                    "GO:0003677": "DNA binding",
                    "GO:0010385": "double-stranded methylated DNA binding",
                    "GO:0004519": "endonuclease activity",
                    "GO:0005525": "GTP binding",
                    "GO:0003924": "GTPase activity",
                    "GO:0044729": "hemi-methylated DNA-binding",
                    "GO:0042802": "identical protein binding",
                    "GO:0015666": "restriction endodeoxyribonuclease activity"
                },
                "Biological Process": {
                    "GO:0006308": "DNA catabolic process",
                    "GO:0009307": "DNA restriction-modification system"
            }
        }
    }
}

@pytest.fixture
def target_id_plus_seq_Q9NU22(minimal_hmmalign_lines_fixture_Q9NU22):
    """Returns target ID and sequence as tuple."""
    lines = minimal_hmmalign_lines_fixture_Q9NU22.splitlines()
    try:
        line = next(line for line in lines if line.startswith("sp|Q9NU22|MDN1_HUMANtarget//325-451"))
    except StopIteration:
        raise ValueError("Target ID 'sp|Q9NU22|MDN1_HUMANtarget//325-451' not found in hmmalign_lines.")
    id_seq = line.split()[:2]  # Returns [id, seq]
    return id_seq

@pytest.fixture
def conservation_id_plus_seq_Q9NU22_PF07728(minimal_hmmalign_lines_fixture_Q9NU22):
    """Returns conservation ID and sequence as tuple."""
    lines = minimal_hmmalign_lines_fixture_Q9NU22.splitlines()
    try:
        line = next(line for line in lines if line.startswith("Q7US48_RHOBA/138-284"))
    except StopIteration:
        raise ValueError("Conservation ID 'MCRB_ECOLI/196-350' not found in hmmalign_lines.")
    id_seq = line.split()[:2]  # Returns [id, seq]
    return id_seq


@pytest.fixture()
def mock_json_filepaths(tmp_path):
    """Create mock JSON files for testing"""
    cons_file = tmp_path / "conservations.json"
    annot_file = tmp_path / "annotations.json"

    return cons_file, annot_file

### Added for testing miscellaneous/utility functions (convert_*, read_*)
## Sent transfer_dict_populated_disulfid_list_original_structure_Q9NU22 to test_utils.py
## since only convert_lists_to_original_types made use of it and it went to utils, not transfer_annotations.

@pytest.fixture
def mock_open_files():
    def _mock_open_files(mock_files):
        def _opener(file, *args, **kwargs):
            if file in mock_files:
                return mock_open(read_data=mock_files[file])()
            raise FileNotFoundError(f"File not found: {file}")
        return _opener
    return _mock_open_files

@pytest.fixture
def mock_target_info():
    return {
        "sp|Q9NU22|MDN1_HUMAN": {
            "325-451": [(325, 451, ".........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................FQFFAT-----.........--------------rrllscgg....................")]
        }
    }

### Added for testing add_to_transfer_dict and helper functions
@pytest.fixture
def interval_dict_base():
    """Base interval dict structure used by other fixtures"""
    return {
        "sequence": "TESTSEQUENCE",
        "length": 10,
        "hit_start": 325,
        "hit_end": 451,
        "annotations": {
            "positions": {
                "333": {
                    "DISULFID | Test": {
                        "essentials": {"type": "DISULFID"},
                        "hit": True
                    }
                }
            },
            "indices": {"matches": {"333"}, "misses": set()}
        },
        "position_conversion": {
            "target_to_aln": {"333": "18"},
            "aln_to_target": {"18": "333"}
        },
        "annotation_ranges": {}
    }

@pytest.fixture
def interval_dict_empty(interval_dict_base):
    """Full interval dict structure with empty ranges"""
    return interval_dict_base

@pytest.fixture
def interval_dict_with_range(interval_dict_base):
    """Interval dict with existing range"""
    result = dict(interval_dict_base)
    result["annotation_ranges"] = {
        "DISULFID | Test": {
            "ranges": [(333, 333)],
            "positions": {333}
        }
    }
    return result

@pytest.fixture
def interval_dict_with_multiple_ranges(interval_dict_base):
    """Interval dict with multiple non-adjacent ranges"""
    result = dict(interval_dict_base)
    result["annotation_ranges"] = {
        "DISULFID | Test": {
            "ranges": [(333, 334), (336, 337)],
            "positions": {333, 334, 336, 337}
        }
    }
    return result

@pytest.fixture
def interval_dict_with_position_in_range(interval_dict_base):
    """Interval dict with range containing position to be added"""
    result = dict(interval_dict_base)
    result["annotation_ranges"] = {
        "DISULFID | Test": {
            "ranges": [(333, 335)],
            "positions": {333, 334, 335}
        }
    }
    return result


###### Configuration Constants
### Added in Map_and_filter_annot_pos
target_name_mock_Q9NU22 = "sp|Q9NU22|MDN1_HUMAN"
target_hit_start_mock_Q9NU22 = 325
target_hit_end_mock_Q9NU22 = 451
offset_start_mock_Q9NU22 = 196
offset_end_mock_Q9NU22 = 350
entry_mnemo_name_mock_Q9NU22 = "MCRB_ECOLI"

## Extra variables for paired position testing
paired_annot_pos_str_mock_Q9NU22 = "246"
caller_target_pos_str_mock_Q9NU22 = "333"
# - Removed, no use as of now, kept as documentation only
# caller_annot_pos_str_mock = "205"
# annotation_dict = use the get_annotation_dict fixture

### Added Defaults for use in validate_annotations and validate_paired_annotations
### (values change in downstream functions from these cited)
counter_target_pos_mock = None
counter_annot_pos_mock = None

### Added/Modified in Process_annotation
# res_hit = Define in the test
# annotation_dict = use the get_annotation_dict fixture
# counter_target_pos = Define in the test
# counter_annot_pos_str = Define in the test

### Added/Modified in Make_anno_total_dict
# annotation_dict = use the get_annotation_dict fixture
# counter_target_pos = Define in the test
# counter_annot_pos_str = Define in the test

## If called from validate_paired_annotations:
# caller_target_pos = Define in the test
# logger = use the logger fixture
# entry_annotations = use the entry_annotations_* fixture, according to the test

### Added/Modified in Add_to_transfer_dict
entry_primary_accesion_mock_Q9NU22 = "P15005"

# Repeated Anno ID test - Using _return_348 and _return_725 fixtures
entry_primary_accession_mock_repeated_Q9NU22 = "P00720"
entry_mnemo_name_mock_repeated_Q9NU22 = "TPA_HUMAN"

# anno_id = Define in the test, new
# anno_total = Define in the test, new

## If called from a paireable type in process_annotation, also note that:
# paired_annot_pos = Define in the test
# caller_target_pos = Define in the test
# annotation_dict = use the get_annotation_dict fixture

### Added/Modified in _add_single_annotation
# paired_additional_keys/additional_keys = Define in the test
# paired_position_res_hit/hit = Define in the test
# paired_anno_id = Define in the test
# paired_anno_total = Define in the test

### Added in clenup_improve_transfer_dict and helper functions
pfam_id_mock_Q9NU22 = "PF07728"


###T parse_arguments

def test_parse_arguments_required():
    test_args = [
        "-iA", hmmalign_result_mock,
        "-r", resource_dir_mock,
        "-d", domain_accession_mock,
        "-o", output_dir_mock
    ]

    with pytest.raises(SystemExit):
        parse_arguments()

    with pytest.MonkeyPatch().context() as mp:
        mp.setattr("sys.argv", ["script_name.py"] + test_args)
        args = parse_arguments()

    expected = Namespace(
        dom_align="/home/user/results/human/PF07728_hmmalign.sth",
        resource_dir="/home/user/resources/",
        domain_accession="PF07728",
        output_dir="/home/user/results/human/PF07728/",
        eco_codes=[],
        log="logs/transfer_annotations.log"
    )
    assert vars(args) == vars(expected)

def test_parse_arguments_optional():
    test_args = [
        "-iA", hmmalign_result_mock,
        "-r", resource_dir_mock,
        "-d", domain_accession_mock,
        "-o", output_dir_mock,
        "-e", *good_eco_codes_mock, # Unpacking list of ECO codes to simupaired space-separated input
        "-l", log_filepath_mock
    ]

    with pytest.MonkeyPatch().context() as mp:
        mp.setattr("sys.argv", ["script_name.py"] + test_args)
        args = parse_arguments()

    assert args.log == log_filepath_mock
    assert args.eco_codes == good_eco_codes_mock

###T get_pfam_id_from_hmmalign_result

def test_get_pfam_id_from_hmmalign_result():
    assert get_pfam_id_from_hmmalign_result(hmmalign_result_mock) == "PF07728"
    hmmalign_result_mock_2 = "/home/user/results/human/PF00001_hmmalign.sth"
    assert get_pfam_id_from_hmmalign_result(hmmalign_result_mock_2) == "PF00001"

###T get_annotation_filepath

def test_get_annotation_filepath():
    assert get_annotation_filepath(resource_dir_mock, "PF07728") == ("/home/user/resources/PF07728/annotations.json", "/home/user/resources/PF07728/conservations.json")

###T read_files

def test_read_files_basic(annotations_content_binding_fixture_Q9NU22_PF07728, mock_open_files):
    annotations_content = json.dumps(annotations_content_binding_fixture_Q9NU22_PF07728)
    mock_files = {
        "dummy_hmmalign": hmmalign_result_content_mock,
        "dummy_annotations": annotations_content
    }

    with patch("builtins.open", new=mock_open_files(mock_files)):
        hmmalign_lines, annotations = read_files("dummy_hmmalign", "dummy_annotations")
        assert hmmalign_lines == hmmalign_result_content_mock.splitlines()
        assert annotations == annotations_content_binding_fixture_Q9NU22_PF07728

def test_read_files_missing_annotations(mock_open_files):
    mock_files = {
        "dummy_hmmalign": hmmalign_result_content_mock
    }

    with patch("builtins.open", new=mock_open_files(mock_files)):
        hmmalign_lines, annotations = read_files("dummy_hmmalign", "nonexistent_annotations")
        assert hmmalign_lines == hmmalign_result_content_mock.splitlines()
        assert annotations == {"sequence_id": {}}

###T iterate_aligned_sequences

@pytest.fixture
def basic_seq_pair():
    return "ABC", "XYZ", 1, 10, 3, 12

@pytest.fixture
def gapped_seq_pair():
    return "A-BC", "X-YZ", 1, 10, 3, 12

def test_iterate_aligned_sequences_basic(basic_seq_pair):
    source_sequence, target_sequence, source_start, target_start, source_end, target_end = basic_seq_pair
    result = list(iterate_aligned_sequences(source_sequence, target_sequence, source_start, target_start, source_end, target_end))

    expected = [
        (0, 1, 10, "A", "X"),
        (1, 2, 11, "B", "Y"),
        (2, 3, 12, "C", "Z")
    ]

    assert result == expected

def test_iterate_aligned_sequences_with_gaps(gapped_seq_pair):
    source_sequence, target_sequence, source_start, target_start, source_end, target_end = gapped_seq_pair
    result = list(iterate_aligned_sequences(source_sequence, target_sequence, source_start, target_start, source_end, target_end))

    expected = [
        (0, 1, 10, "A", "X"),
        (1, 1, None, "-", "-"), # Gap in target will be None
        (2, 2, 11, "B", "Y"),
        (3, 3, 12, "C", "Z")
    ]

    assert result == expected

def test_iterate_aligned_sequences_early_stop():
    result = list(iterate_aligned_sequences(
        "ABCDEF", "XYZUVW",
        source_start=1, target_start=10,
        source_end=3, target_end=12  # Should stop after "C"/"Z"
    ))

    expected = [
        (0, 1, 10, "A", "X"),
        (1, 2, 11, "B", "Y"),
        (2, 3, 12, "C", "Z")
    ]

    assert result == expected

def test_iterate_aligned_sequences_empty():
    result = list(iterate_aligned_sequences("", "", 0, 0, 1, 1))
    assert result == []

def test_iterate_aligned_sequences_real_Q9NU22():
    # A - MCRB_ECOLI ; B - sp|Q9NU22|MDN1_HUMAN
    source_sequence = ".........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYCCN"
    target_sequence = ".........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCT"
    source_start, target_start = 196, 325
    source_end, target_end = 350, 451
    result = list(iterate_aligned_sequences(source_sequence, target_sequence, source_start, target_start, source_end, target_end))

    expected_indices = {
        14: (14, 201, 329, "G", "G"),
        15: (15, 202, 330, "P", "P"),
        18: (18, 205, 333, "C", "C"),
        72: (72, 246, 373, "C", "C"),
    }
    for idx, expected_tuple in expected_indices.items():
        assert result[idx] == expected_tuple

###T extract_target_info_from_hmmalign

def test_extract_target_info_empty_lines(logger, multi_logger):
    """Test empty hmmalign lines."""
    mock_lines = ["# STOCKHOLM 1.0"]
    result = extract_target_info_from_hmmalign(logger, multi_logger, mock_lines)
    assert result == {}
    multi_logger.assert_not_called()

def test_extract_target_info_single_target(logger, multi_logger):
    """Test single target line extraction."""
    mock_lines = [
        "# STOCKHOLM 1.0",
        "sp|Q9NU22|MDN1_HUMANtarget//325-451  VLLEGPIGCGKTSLVE"
    ]
    expected = {
        "sp|Q9NU22|MDN1_HUMAN": {
            "325-451": [(325, 451, "VLLEGPIGCGKTSLVE")]
        }
    }
    result = extract_target_info_from_hmmalign(logger, multi_logger, mock_lines)
    assert result == expected
    multi_logger.assert_not_called()

def test_extract_target_info_multiple_targets(logger, multi_logger):
    """Test multiple target lines extraction."""
    mock_lines = [
        "# STOCKHOLM 1.0",
        "sp|Q9NU22|MDN1_HUMANtarget//325-451  VLLEGPIGCGKTSLVE",
        "sp|Q9NU22|MDN1_HUMANtarget//452-580  DVFLIGPPGPLRRSIAM"
    ]
    expected = {
        "sp|Q9NU22|MDN1_HUMAN": {
            "325-451": [(325, 451, "VLLEGPIGCGKTSLVE")],
            "452-580": [(452, 580, "DVFLIGPPGPLRRSIAM")]
        }
    }
    result = extract_target_info_from_hmmalign(logger, multi_logger, mock_lines)
    assert result == expected
    multi_logger.assert_not_called()

def test_extract_target_info_invalid_format(logger, multi_logger):
    """Test handling of invalid format lines."""
    mock_lines = [
        "# STOCKHOLM 1.0",
        "sp|Q9NU22|MDN1_HUMANtarget/invalid_format  SEQUENCE"
    ]
    result = extract_target_info_from_hmmalign(logger, multi_logger, mock_lines)
    assert result == {}
    multi_logger.assert_called_once_with(
        "warning",
        "TRANSFER_ANNOTS --- EXTRACT_TARGET_INFO --- Invalid format in line: %s. Error: %s",
        "sp|Q9NU22|MDN1_HUMANtarget/invalid_format  SEQUENCE",
        ANY  # The exact error message may vary
    )


###T setup_for_conservations_only

def test_setup_for_conservations_only_basic(logger, multi_logger):
    """Test basic case with single target"""
    mock_lines = [
        "# STOCKHOLM 1.0",
        "sp|Q9NU22|MDN1_HUMANtarget//325-451  VLLEGPIGCGKTSLVE"
    ]
    expected = {
        "PF07728": {
            "sequence_id": {
                "sp|Q9NU22|MDN1_HUMAN": {
                    "hit_intervals": {
                        "325-451": {
                            "sequence": "VLLEGPIGCGKTSLVE",
                            "length": 16,
                            "hit_start": 325,
                            "hit_end": 451,
                            "annotations": {},
                            "conservations": {
                                "positions": {},
                                "indices": {"matches": set(), "misses": set()}
                            },
                            "position_conversion": {
                                "target_to_aln": {},
                                "aln_to_target": {}
                            },
                            "annotation_ranges": {},
                            "conservation_ranges": {}
                        }
                    }
                }
            }
        }
    }

    result = setup_for_conservations_only(logger, multi_logger, mock_lines, "PF07728")
    assert result == expected

def test_setup_for_conservations_only_empty_lines(logger, multi_logger):
    """Test empty hmmalign lines case"""
    mock_lines = ["# STOCKHOLM 1.0"]
    result = setup_for_conservations_only(logger, multi_logger, mock_lines, "PF07728")
    assert result == {}
    multi_logger.assert_called_once_with(
    "error",
    "TRANSFER_ANNOTS --- SETUP_CONS_ONLY --- No target sequences found in hmmalign lines!"
    )

def test_setup_for_conservations_only_invalid_format(logger, multi_logger):
    """Test invalid format handling"""
    mock_lines = [
        "# STOCKHOLM 1.0",
        "sp|Q9NU22|MDN1_HUMANtarget/invalid_format  SEQUENCE"
    ]
    result = setup_for_conservations_only(logger, multi_logger, mock_lines, "PF07728")
    assert result == {}
    multi_logger.assert_called_with(
        "error",
        "TRANSFER_ANNOTS --- SETUP_CONS_ONLY --- No target sequences found in hmmalign lines!"
    )


###T find_and_map_annots

def test_find_and_map_annots_basic_case(annotations_content_binding_fixture_Q9NU22_PF07728, logger, multi_logger, mock_target_info, target_sequence_Q9NU22, annot_sequence_Q9NU22):
    """Test basic case with single binding annotation"""
    mock_lines = [
        "# STOCKHOLM 1.0\n",
        "MCRB_ECOLI/196-350      .........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYCPN............GVGFRRKDGIFYNFCQqak...........eqpeKKYIFIIDEINRANlS........K....VFG....E......VMM.LMEHD.K...RGENWSVPLTYSENDE.eRFYVp..........................enVYIIGLMNTAD.........RSLAVVDYALRRRF............................\n",
        "sp|Q9NU22|MDN1_HUMANtarget//325-451                     .........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................FQFFAT-----.........--------------rrllscgg....................\n"
    ]

    with patch("transfer_annotations.extract_target_info_from_hmmalign", return_value=mock_target_info) as mock_extract, \
         patch("transfer_annotations.map_and_filter_annot_pos") as mock_map_and_filter:

        find_and_map_annots(
            logger=logger,
            multi_logger=multi_logger,
            hmmalign_lines=mock_lines,
            annotations=annotations_content_binding_fixture_Q9NU22_PF07728,
            good_eco_codes=[]
        )

        mock_extract.assert_called_once_with(logger, multi_logger, mock_lines)

        mock_map_and_filter.assert_called_with(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=[],
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=annotations_content_binding_fixture_Q9NU22_PF07728["MCRB_ECOLI"],
            transfer_dict={},
            processed_annotations=set(),
        )

def test_find_and_map_annots_no_target_sequence(annotations_content_binding_fixture_Q9NU22_PF07728, logger, multi_logger):
    """Test case with no matching annotations"""
    mock_lines = [
        "# STOCKHOLM 1.0\n",
        "VWA8_HUMAN/105-261    .........DVFLIGPPGPLRRSIAM.QYLELT..............KREVEYIALSR.DT..TETDLKQRREIR\n"
    ]
    with patch("transfer_annotations.extract_target_info_from_hmmalign", return_value={}):
        transfer_dict = find_and_map_annots(
            logger=logger,
            multi_logger=multi_logger,
            hmmalign_lines=mock_lines,
            annotations=annotations_content_binding_fixture_Q9NU22_PF07728,
            good_eco_codes=[]
        )

        assert not transfer_dict

        # Verify the logger captured the correct error message
        multi_logger.assert_called_once_with(
            "error",
            "TRANSFER_ANNOTS --- FIND_AND_MAP --- No target sequences found in hmmalign lines!"
        )

###T read_conservations_and_annotations

def test_read_conservations_and_annotations_success(
    tmp_path, mock_json_filepaths, conservations_content_Q9NU22_PF07728, annotations_content_disulfid_fixture_Q9NU22_PF07728):
    # Create temporary files
    cons_file, annot_file = mock_json_filepaths

    with open(cons_file, "w") as f:
        json.dump(conservations_content_Q9NU22_PF07728, f)
    with open(annot_file, "w") as f:
        json.dump(annotations_content_disulfid_fixture_Q9NU22_PF07728, f)

    # Test function
    conservations, annotations = read_conservations_and_annotations(str(cons_file), str(annot_file))

    assert conservations == conservations_content_Q9NU22_PF07728
    assert annotations == annotations_content_disulfid_fixture_Q9NU22_PF07728

def test_read_conservations_and_annotations_conservations_only(tmp_path):
    """Test handling when only conservations file exists"""
    cons_file = tmp_path / "conservations.json"
    mock_conservations = {"sequence_id/1-100": {"1": 0.5, "2": 0.8}}

    with open(cons_file, "w") as f:
        json.dump(mock_conservations, f)

    conservations, annotations = read_conservations_and_annotations(
        str(cons_file),
        str(tmp_path / "nonexistent_annotations.json")
    )

    assert conservations == mock_conservations
    assert annotations == {"sequence_id": {}}

def test_read_conservations_and_annotations_missing_files(tmp_path):
    """Test handling of missing files - should return empty dicts."""
    conservations, annotations = read_conservations_and_annotations(
        str(tmp_path / "nonexistent_cons.json"),
        str(tmp_path / "nonexistent_annot.json")
    )

    assert conservations == {"sequence_id/range": {}}
    assert annotations == {"sequence_id": {}}

def test_read_conservations_and_annotations_one_file_missing(tmp_path, mock_json_filepaths):
    """Test handling when only one file exists."""
    cons_file, _ = mock_json_filepaths

    # Create only conservations file
    with open(cons_file, "w") as f:
        json.dump({"some": "data"}, f)

    conservations, annotations = read_conservations_and_annotations(
        str(cons_file),
        str(tmp_path / "nonexistent_annot.json")
    )

    assert conservations == {"some": "data"}
    assert annotations == {"sequence_id": {}}

def test_read_conservations_and_annotations_invalid_json(tmp_path, mock_json_filepaths):
    # Create file with invalid JSON
    cons_file, annot_file = mock_json_filepaths

    with open(cons_file, "w") as f:
        f.write("invalid json")
    with open(annot_file, "w") as f:
        json.dump({}, f)

    with pytest.raises(json.JSONDecodeError):
        read_conservations_and_annotations(str(cons_file), str(annot_file))

def test_read_conservations_and_annotations_empty_files(mock_json_filepaths):
    """Test handling of empty but valid JSON files"""
    cons_file, annot_file = mock_json_filepaths

    with open(cons_file, "w", encoding="utf-8") as f:
        json.dump({}, f)
    with open(annot_file, "w", encoding="utf-8") as f:
        json.dump({}, f)

    conservations, annotations = read_conservations_and_annotations(str(cons_file), str(annot_file))
    assert conservations == {"sequence_id/range": {}}
    assert annotations == {"sequence_id": {}}

###T parse_go_annotations

def test_parse_go_annotations_multiple_terms():
    """Test parsing GO terms with multiple terms and sources."""
    go_column = "GO:0000027(PANTHER)|GO:0000055(InterPro)|GO:0005634()|GO:0030687(PANTHER)"
    result = parse_go_annotations(go_column)
    assert result == ["GO:0000027", "GO:0000055", "GO:0005634", "GO:0030687"]

def test_parse_go_annotations_single_term():
    """Test parsing single GO term."""
    go_column = "GO:0016887(InterPro)"
    result = parse_go_annotations(go_column)
    assert result == ["GO:0016887"]

def test_parse_go_annotations_empty():
    """Test parsing empty GO annotations."""
    assert parse_go_annotations("") == []
    assert parse_go_annotations("-") == []
    assert parse_go_annotations("  ") == []
    assert parse_go_annotations(None) == []


###T check_interval_overlap

def test_check_interval_overlap_exact_match(multi_logger):
    """Test when intervals match exactly."""
    assert check_interval_overlap(multi_logger, 100, 200, 100, 200) == True

def test_check_interval_overlap_within_margins(multi_logger):
    """Test when interval is within allowed margins."""
    # Hit interval 100-200 (length 101)
    # Default 10% margin = 10 positions
    # Allowed range: [90-210]

    # Exactly at margin boundaries
    assert check_interval_overlap(multi_logger, 90, 210, 100, 200) == True

    # One end within margin, other end inside interval
    assert check_interval_overlap(multi_logger, 90, 150, 100, 200) == True  # Start at margin
    assert check_interval_overlap(multi_logger, 150, 210, 100, 200) == True  # End at margin

    # Just outside margins
    assert check_interval_overlap(multi_logger, 89, 200, 100, 200) == False  # Start beyond margin
    assert check_interval_overlap(multi_logger, 100, 211, 100, 200) == False  # End beyond margin

def test_check_interval_overlap_outside_margins(multi_logger):
    """Test when interval is outside allowed margins."""
    assert check_interval_overlap(multi_logger, 50, 150, 200, 300) == False
    assert check_interval_overlap(multi_logger, 250, 350, 100, 200) == False

def test_check_interval_overlap_invalid_inputs(multi_logger):
    """Test handling of invalid inputs."""
    check_interval_overlap(multi_logger, -1, 100, 100, 200)
    multi_logger.assert_called_once_with(
        "error",
        "TRANSFER_ANNOTS --- CHECK_INTERVAL --- All positions must be non-negative numbers"
    )
    check_interval_overlap(multi_logger, 100, 200, 100, 200, 1.5)
    multi_logger.assert_called_with(
        "error",
        "TRANSFER_ANNOTS --- CHECK_INTERVAL --- margin_percent must be a float between 0 and 1 (exclusive, inclusive)"
    )

def test_check_interval_overlap_edge_margins(multi_logger):
    """Test edge cases for margin percentages."""
    check_interval_overlap(multi_logger, 100, 200, 100, 200, 0.0)
    multi_logger.assert_called_once_with(
        "error",
        "TRANSFER_ANNOTS --- CHECK_INTERVAL --- margin_percent must be a float between 0 and 1 (exclusive, inclusive)"
    )
    assert check_interval_overlap(multi_logger, 90, 210, 100, 200, 1.0) == True

###T gather_go_terms_for_target

def test_gather_go_terms_for_target_pfam_match(multi_logger, iprscan_df_Q9NU22_PF07728):
    """Test gathering GO terms when matching PFAM ID"""
    with patch("os.path.exists", return_value=True), \
         patch("builtins.open", mock_open()), \
         patch("pandas.read_csv", return_value=iprscan_df_Q9NU22_PF07728), \
         patch("transfer_annotations.check_interval_overlap", return_value=False), \
         patch("transfer_annotations.parse_go_annotations") as mock_parse:

        mock_parse.return_value = ["GO:0005524", "GO:0016887"]

        result = gather_go_terms_for_target(
            multi_logger=multi_logger,
            target_name="sp|Q9NU22|MDN1_HUMAN",
            pfam_id="PF07728",
            go_terms_dir="/mock/path",
            interpro_conv_id="IPR011704",
            hit_start=325,
            hit_end=451
        )

        assert result == {"GO:0005524", "GO:0016887"}
        mock_parse.assert_called_once_with("GO:0005524(InterPro)|GO:0016887(InterPro)")

def test_gather_go_terms_for_target_interpro_match(multi_logger, iprscan_df_Q9NU22_PF07728):
    """Test gathering GO terms when matching InterPro ID"""
    with patch("os.path.exists", return_value=True), \
         patch("builtins.open", mock_open()), \
         patch("pandas.read_csv", return_value=iprscan_df_Q9NU22_PF07728), \
         patch("transfer_annotations.check_interval_overlap", return_value=False), \
         patch("transfer_annotations.parse_go_annotations") as mock_parse:

        mock_parse.return_value = ["GO:0042623"]

        result = gather_go_terms_for_target(
            multi_logger=multi_logger,
            target_name="sp|Q9NU22|MDN1_HUMAN",
            pfam_id="PF00000",
            go_terms_dir="/mock/path",
            interpro_conv_id="IPR027417",
            hit_start=300,
            hit_end=500
        )

        assert result == {"GO:0042623"}
        mock_parse.assert_called_once_with("GO:0042623(InterPro)")

def test_gather_go_terms_for_target_interval_only(multi_logger, iprscan_df_Q9NU22_PF07728):
    """Test gathering GO terms when only interval matches"""
    with patch("os.path.exists", return_value=True), \
         patch("builtins.open", mock_open()), \
         patch("pandas.read_csv", return_value=iprscan_df_Q9NU22_PF07728), \
         patch("transfer_annotations.check_interval_overlap", return_value=True), \
         patch("transfer_annotations.parse_go_annotations") as mock_parse:

        mock_parse.return_value = ["GO:0016887"]

        result = gather_go_terms_for_target(
            multi_logger=multi_logger,
            target_name="sp|Q9NU22|MDN1_HUMAN",
            pfam_id="PF00000",
            go_terms_dir="/mock/path",
            interpro_conv_id="IPR000000",
            hit_start=320,
            hit_end=635
        )

        assert result == {"GO:0016887"}

def test_gather_go_terms_for_target_no_matches(multi_logger, iprscan_df_Q9NU22_PF07728):
    """Test when no entries match any criteria"""
    with patch("os.path.exists", return_value=True), \
         patch("builtins.open", mock_open()), \
         patch("pandas.read_csv", return_value=iprscan_df_Q9NU22_PF07728), \
         patch("transfer_annotations.check_interval_overlap", return_value=False):

        result = gather_go_terms_for_target(
            multi_logger=multi_logger,
            target_name="sp|Q9NU22|MDN1_HUMAN",
            pfam_id="PF00000",
            go_terms_dir="/mock/path",
            interpro_conv_id="IPR000000",
            hit_start=1,
            hit_end=100
        )

        assert result == set()
        multi_logger.assert_called_once_with(
            "warning",
            "TRANSFER_ANNOTS --- GO_TERMS_TARGET --- No usable InterProScan lines in iprscan.tsv for %s-%s",
            "PF00000", "sp-Q9NU22-MDN1_HUMAN"
        )

def test_gather_go_terms_for_target_file_not_found(multi_logger):
    """Test handling of missing iprscan.tsv file"""
    with patch("os.path.exists", return_value=False):
        result = gather_go_terms_for_target(
            multi_logger=multi_logger,
            target_name="nonexistent",
            pfam_id="PF00000",
            go_terms_dir="/mock/path",
            interpro_conv_id="IPR000000",
            hit_start=1,
            hit_end=100
        )

        assert result == set()
        multi_logger.assert_called_once_with(
            "warning",
            "TRANSFER_ANNOTS --- GO_TERMS_TARGET --- Missing iprscan.tsv file for Pfam ID %s and Target Name %s at %s",
            "PF00000", "nonexistent", "/mock/path/nonexistent/iprscan.tsv"
        )

def test_gather_go_terms_for_target_empty_file(multi_logger, mock_empty_df):
    """Test handling of empty but valid TSV file"""
    with patch("os.path.exists", return_value=True), \
         patch("builtins.open", mock_open()), \
         patch("pandas.read_csv", return_value=mock_empty_df), \
         patch("transfer_annotations.check_interval_overlap", return_value=False):

        result = gather_go_terms_for_target(
            multi_logger=multi_logger,
            target_name="sp|Q9NU22|MDN1_HUMAN",
            pfam_id="PF07728",
            go_terms_dir="/mock/path",
            interpro_conv_id="IPR011704",
            hit_start=325,
            hit_end=451
        )

        assert result == set()
        multi_logger.assert_called_once_with(
            "warning",
            "TRANSFER_ANNOTS --- GO_TERMS_TARGET --- No usable InterProScan lines in iprscan.tsv for %s-%s",
            "PF07728", "sp-Q9NU22-MDN1_HUMAN")

def test_gather_go_terms_for_target_multiple_matches(multi_logger, iprscan_df_Q9NU22_PF07728):
    """Test gathering GO terms when both accession and interval match"""
    with patch("os.path.exists", return_value=True), \
         patch("builtins.open", mock_open()), \
         patch("pandas.read_csv", return_value=iprscan_df_Q9NU22_PF07728), \
         patch("transfer_annotations.check_interval_overlap", return_value=True), \
         patch("transfer_annotations.parse_go_annotations") as mock_parse:

        # Need to provide enough return values for all potential matches in iprscan_df_Q9NU22_PF07728
        mock_parse.side_effect = [
            ["GO:0005524", "GO:0016887"],  # PF07728 match
            ["GO:0042623"],                 # Gene3D match
            ["GO:0016887"],                 # SMART match
        ]

        result = gather_go_terms_for_target(
            multi_logger=multi_logger,
            target_name="sp|Q9NU22|MDN1_HUMAN",
            pfam_id="PF07728",
            go_terms_dir="/mock/path",
            interpro_conv_id="IPR011704",
            hit_start=325,
            hit_end=451
        )

        # All GO terms from all matching rows should be combined
        assert result == {"GO:0005524", "GO:0016887", "GO:0042623"}
        assert mock_parse.call_count == 3

def test_gather_go_terms_for_target_dash_go_terms(multi_logger, iprscan_df_Q9NU22_PF07728):
    """Test handling of "-" GO annotations"""
    mock_df_dash = iprscan_df_Q9NU22_PF07728.copy()
    mock_df_dash.loc[0, "GO Annotations"] = "-"

    with patch("os.path.exists", return_value=True), \
         patch("builtins.open", mock_open()), \
         patch("pandas.read_csv", return_value=mock_df_dash), \
         patch("transfer_annotations.check_interval_overlap", return_value=False), \
         patch("transfer_annotations.parse_go_annotations", wraps=parse_go_annotations):

        result = gather_go_terms_for_target(
            multi_logger=multi_logger,
            target_name="sp|Q9NU22|MDN1_HUMAN",
            pfam_id="PF07728",
            go_terms_dir="/mock/path",
            interpro_conv_id="IPR011704",
            hit_start=325,
            hit_end=451
        )

        assert result == set()

## Integration tests so to speak
def test_gather_go_terms_for_target(multi_logger, go_terms_dir_mock):
    go_terms = gather_go_terms_for_target(
        multi_logger=multi_logger,
        target_name="target_name",
        pfam_id="PF07728",
        go_terms_dir=go_terms_dir_mock,
        interpro_conv_id="IPR011704",
        hit_start=325,
        hit_end=451
    )
    assert go_terms == {"GO:0005524", "GO:0016887"}

def test_gather_go_terms_for_target_no_go(multi_logger, go_terms_dir_mock):
    go_terms = gather_go_terms_for_target(
        multi_logger=multi_logger,
        target_name="target_name_no_go",
        pfam_id="PF07728",
        go_terms_dir=go_terms_dir_mock,
        interpro_conv_id="IPR011704",
        hit_start=325,
        hit_end=451
    )
    assert go_terms == set()

###T get_alignment_sequences

def test_get_alignment_sequences(minimal_hmmalign_lines_fixture_Q9NU22, target_id_plus_seq_Q9NU22, conservation_id_plus_seq_Q9NU22_PF07728, logger, multi_logger):
    hmmalign_lines_list = minimal_hmmalign_lines_fixture_Q9NU22.splitlines()
    target_seq, conservation_seq = get_alignment_sequences(
        hmmalign_lines_list,
        "sp|Q9NU22|MDN1_HUMANtarget//325-451",
        "Q7US48_RHOBA/138-284",
        logger,
        multi_logger
    )

    _, target_seq_expected = target_id_plus_seq_Q9NU22
    _, conservation_seq_expected = conservation_id_plus_seq_Q9NU22_PF07728

    assert target_seq == target_seq_expected
    assert conservation_seq == conservation_seq_expected

def test_get_alignment_sequence_no_target(minimal_hmmalign_lines_fixture_Q9NU22, logger, multi_logger):
    hmmalign_lines_list = minimal_hmmalign_lines_fixture_Q9NU22.splitlines()

    # Remove the line with the target ID
    hmmalign_lines_list = [line for line in hmmalign_lines_list if not line.startswith("sp|Q9NU22|MDN1_HUMANtarget//325-451")]

    get_alignment_sequences(
        hmmalign_lines_list,
        "sp|Q9NU22|MDN1_HUMANtarget//325-451",
        "Q7US48_RHOBA/138-284",
        logger,
        multi_logger
    )

    multi_logger.assert_called_once_with(
        "error",
        "TRANSFER_ANNOTS --- GET_ALIGN_SEQS --- Target ID '%s' not found in hmmalign_lines.",
        "sp|Q9NU22|MDN1_HUMANtarget//325-451"
    )


def test_get_alignment_sequence_no_conservation(minimal_hmmalign_lines_fixture_Q9NU22, logger, multi_logger):
    hmmalign_lines_list = minimal_hmmalign_lines_fixture_Q9NU22.splitlines()

    # Remove the line with the conservation ID
    hmmalign_lines_list = [line for line in hmmalign_lines_list if not line.startswith("Q7US48_RHOBA/138-284")]

    get_alignment_sequences(
        hmmalign_lines_list,
        "sp|Q9NU22|MDN1_HUMANtarget//325-451",
        "Q7US48_RHOBA/138-284",
        logger,
        multi_logger
    )

    multi_logger.assert_called_once_with(
        "error",
        "TRANSFER_ANNOTS --- GET_ALIGN_SEQS --- Conservation ID '%s' not found in hmmalign_lines.",
        "Q7US48_RHOBA/138-284"
    )

###T populate_conservation

def test_populate_conservation_base(transfer_dict_populated_disulfid_inside_cleanup_Q9NU22, conservations_content_Q9NU22_PF07728, conservation_id_plus_seq_Q9NU22_PF07728, target_id_plus_seq_Q9NU22, logger):
    """Test conservation data population."""
    # Test data
    conservation_key, conservation_seq = conservation_id_plus_seq_Q9NU22_PF07728
    conserved_positions = ["143", "146", "149", "267", "283", "284"]
    interval_key = "325-451"
    target_seq = target_id_plus_seq_Q9NU22[1]

    populate_conservation(
        transfer_dict=transfer_dict_populated_disulfid_inside_cleanup_Q9NU22,
        pfam_id="PF07728",
        target_name="sp|Q9NU22|MDN1_HUMAN",
        target_seq=target_seq,
        conservation_seq=conservation_seq,
        conservation_key=conservation_key,
        conservations=conservations_content_Q9NU22_PF07728,
        conserved_positions=conserved_positions,
        target_hit_start=325,
        target_hit_end=451,
        interval_key=interval_key,
        conservation_start=138,
        conservation_end=284,
        logger=logger
    )

    # Verify conservation data was populated correctly
    interval_data = transfer_dict_populated_disulfid_inside_cleanup_Q9NU22["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"][interval_key]

    # Test positions data
    positions_data = interval_data["conservations"]["positions"]
    assert "329" in positions_data
    assert positions_data["329"]["conservation"] == 0.9853
    assert positions_data["329"]["hit"] is True

    assert "332" in positions_data
    assert positions_data["332"]["conservation"] == 0.8806
    assert positions_data["332"]["hit"] is True

    assert "335" in positions_data
    assert positions_data["335"]["conservation"] == 0.9562
    assert positions_data["335"]["hit"] is True

    # Test indices data
    indices_data = interval_data["conservations"]["indices"]
    assert indices_data["matches"] == {"329", "332", "335"}
    assert indices_data["misses"] == set()


def test_populate_conservation_gap_position(transfer_dict_populated_disulfid_inside_cleanup_Q9NU22, conservations_content_Q9NU22_PF07728, conservation_id_plus_seq_Q9NU22_PF07728, target_id_plus_seq_Q9NU22, logger):
    """Test conservation data when target sequence has gap at conserved position."""
    original_target_seq = target_id_plus_seq_Q9NU22[1]
    conservation_key, conservation_seq = conservation_id_plus_seq_Q9NU22_PF07728
    conserved_positions = ["143", "146", "149", "267", "283", "284"]
    interval_key = "325-451"

    # Insert gap at index position 20 (cons_pos 149 -> would be target_pos 335)
    offset_for_gap = 20  # Index position for cons_pos 149
    mod_target_seq = original_target_seq[:offset_for_gap] + "--" + original_target_seq[offset_for_gap:]

    populate_conservation(
        transfer_dict=transfer_dict_populated_disulfid_inside_cleanup_Q9NU22,
        pfam_id="PF07728",
        target_name="sp|Q9NU22|MDN1_HUMAN",
        target_seq=mod_target_seq,  # Use modified sequence with gap
        conservation_seq=conservation_seq,
        conservation_key=conservation_key,
        conservations=conservations_content_Q9NU22_PF07728,
        conserved_positions=conserved_positions,
        target_hit_start=325,
        target_hit_end=451,
        interval_key=interval_key,
        conservation_start=138,
        conservation_end=284,
        logger=logger
    )
    # Verify debug msg about missing pos
    assert logger.debug.called
    # Get the format string and arguments separately
    debug_msg = logger.debug.call_args[0][0]  # Format string
    debug_args = logger.debug.call_args[0][1:]  # Arguments passed to format string
    assert "TRANSFER_ANNOTS --- POP_CONS --- Pfam ID: %s - Target Name: %s - Conserv-Rep Name: %s - Missing conserved positions: %s" == debug_msg

    # The missing positions will be in the last argument
    missing_positions_str = str(debug_args[3])  # Get the actual missing positions set

    # Parse the missing positions from the set string
    # The string will be something like "{'149', '283', '284'}"
    actual_missing = set(pos.strip("'\" {}") for pos in missing_positions_str.split(","))
    expected_missing = {"149", "283", "284"}

    assert expected_missing == actual_missing

    interval_data = transfer_dict_populated_disulfid_inside_cleanup_Q9NU22["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"][interval_key]
    positions_data = interval_data["conservations"]["positions"]
    # Verify cons_pos 143 -> target_pos 329 hits
    assert positions_data["329"]["hit"] is True
    assert positions_data["329"]["conservation"] == 0.9853
    # Verify cons_pos 149 -> target_pos 335 misses due to inserted gap
    assert "335" not in positions_data

def test_populate_conservation_mismatch_position(
    transfer_dict_populated_disulfid_inside_cleanup_Q9NU22,
    conservations_content_Q9NU22_PF07728,
    conservation_id_plus_seq_Q9NU22_PF07728,
    target_id_plus_seq_Q9NU22,
    logger
):
    """Test conservation data when target sequence has different residue than conserved sequence."""
    original_target_seq = target_id_plus_seq_Q9NU22[1]
    conservation_key, conservation_seq = conservation_id_plus_seq_Q9NU22_PF07728
    conserved_positions = ["143", "146", "149"]
    interval_key = "325-451"

    # Change residue at index position 20 (cons_pos 149 -> target_pos 335)
    mod_target_seq_list = list(original_target_seq)
    mod_target_seq_list[20] = "R"  # Different amino acid
    mod_target_seq = "".join(mod_target_seq_list)

    populate_conservation(
        transfer_dict=transfer_dict_populated_disulfid_inside_cleanup_Q9NU22,
        pfam_id="PF07728",
        target_name="sp|Q9NU22|MDN1_HUMAN",
        target_seq=mod_target_seq,  # Use modified sequence with gap
        conservation_seq=conservation_seq,
        conservation_key=conservation_key,
        conservations=conservations_content_Q9NU22_PF07728,
        conserved_positions=conserved_positions,
        target_hit_start=325,
        target_hit_end=451,
        interval_key=interval_key,
        conservation_start=138,
        conservation_end=284,
        logger=logger
    )

    interval_data = transfer_dict_populated_disulfid_inside_cleanup_Q9NU22["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"][interval_key]
    positions_data = interval_data["conservations"]["positions"]

    # Position should be present but marked as miss
    assert "335" in positions_data
    assert positions_data["335"]["conservation"] == 0.9562
    assert positions_data["335"]["hit"] is False
    assert "335" in interval_data["conservations"]["indices"]["misses"]


###T populate_go_data_for_annotations

def test_populate_go_data_basic_match(
    transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22,
    annotations_content_Q9NU22_PF07728
):
    """Test GO term population with matching terms."""
    # Target GOs from test_gather_go_terms_for_target_multiple_matches
    target_gos = {"GO:0005524", "GO:0016887", "GO:0042623"}

    populate_go_data_for_annotations(
        transfer_dict=transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22,
        pfam_id="PF07728",
        target_name="sp|Q9NU22|MDN1_HUMAN",
        annotations=annotations_content_Q9NU22_PF07728,
        go_terms_annot_key="0",
        target_gos=target_gos
    )

    anno_data = transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]["annotations"]["positions"]["333"]["DISULFID | Intrachain (with C-246); in linked form"]

    assert "GO" in anno_data
    assert "MCRB_ECOLI" in anno_data["GO"]
    assert len(anno_data["GO"]["MCRB_ECOLI"]["terms"]) == 2 # Missing GO:0042623, not in annotations
    assert anno_data["GO"]["MCRB_ECOLI"]["jaccard_index"] == 0.1538


def test_populate_go_data_no_matches(
    transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22,
    annotations_content_Q9NU22_PF07728
):
    """Test GO term population with no matching terms."""
    target_gos = {"GO:0000000"}  # Non-existent term

    populate_go_data_for_annotations(
        transfer_dict=transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22,
        pfam_id="PF07728",
        target_name="sp|Q9NU22|MDN1_HUMAN",
        annotations=annotations_content_Q9NU22_PF07728,
        go_terms_annot_key="0",
        target_gos=target_gos
    )

    interval_data = transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
    anno_data = interval_data["annotations"]["positions"]["333"]["DISULFID | Intrachain (with C-246); in linked form"]

    assert "GO" in anno_data
    assert "MCRB_ECOLI" in anno_data["GO"]
    assert anno_data["GO"]["MCRB_ECOLI"]["terms"] == {}
    assert anno_data["GO"]["MCRB_ECOLI"]["jaccard_index"] == 0.0

def test_populate_go_data_empty_target_gos(
    transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22,
    annotations_content_Q9NU22_PF07728
):
    """Test GO term population with empty target GOs."""
    target_gos = set()

    populate_go_data_for_annotations(
        transfer_dict=transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22,
        pfam_id="PF07728",
        target_name="sp|Q9NU22|MDN1_HUMAN",
        annotations=annotations_content_Q9NU22_PF07728,
        go_terms_annot_key="0",
        target_gos=target_gos
    )

    interval_data = transfer_dict_populated_disulfid_post_conservation_inside_cleanup_Q9NU22["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
    anno_data = interval_data["annotations"]["positions"]["333"]["DISULFID | Intrachain (with C-246); in linked form"]

    assert "GO" in anno_data
    assert "MCRB_ECOLI" in anno_data["GO"]
    assert anno_data["GO"]["MCRB_ECOLI"]["terms"] == {}
    assert anno_data["GO"]["MCRB_ECOLI"]["jaccard_index"] == 0.0



###T cleanup_improve_transfer_dict

def test_cleanup_improve_transfer_dict_basic(
    logger,
    multi_logger,
    transfer_dict_populated_disulfid_Q9NU22,
    conservations_content_Q9NU22_PF07728,
    annotations_content_disulfid_fixture_Q9NU22_PF07728,
    mapping_content_Q9NU22_and_H0YB80_domains,
    target_gos_Q9NU22_MCRB_ECOLI,
    target_id_plus_seq_Q9NU22,
    conservation_id_plus_seq_Q9NU22_PF07728,
    transfer_dict_populated_disulfid_post_conservation_Q9NU22,
    transfer_dict_populated_disulfid_post_gos_Q9NU22
):
    # Define mocks that directly apply fixture data
    def mock_populate_conservation(transfer_dict, **kwargs):
        target_path = transfer_dict["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
        source_path = transfer_dict_populated_disulfid_post_conservation_Q9NU22["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]

        # Apply complete conservation and position tracking state from fixture
        target_path["conservations"] = copy.deepcopy(source_path["conservations"])
        target_path["position_conversion"] = copy.deepcopy(source_path["position_conversion"])
        target_path["conservation_ranges"] = copy.deepcopy(source_path["conservation_ranges"])

    def mock_populate_go_data(transfer_dict, **kwargs):
        # Get paths
        target_path = transfer_dict["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
        source_path = transfer_dict_populated_disulfid_post_gos_Q9NU22["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]

        # Apply complete annotation state
        for pos in ["333", "373"]:
            anno_key = list(target_path["annotations"]["positions"][pos].keys())[0]
            target_path["annotations"]["positions"][pos][anno_key]["GO"] = \
                copy.deepcopy(source_path["annotations"]["positions"][pos][anno_key]["GO"])

    mapping_df = pd.read_csv(StringIO(mapping_content_Q9NU22_and_H0YB80_domains), sep="\t")

    with patch(
        "transfer_annotations.read_conservations_and_annotations",
        return_value=(conservations_content_Q9NU22_PF07728, annotations_content_disulfid_fixture_Q9NU22_PF07728)
    ), patch(
        "pandas.read_csv",
        return_value=mapping_df
    ), patch(
        "transfer_annotations.gather_go_terms_for_target",
        return_value=target_gos_Q9NU22_MCRB_ECOLI
    ), patch(
        "transfer_annotations.get_alignment_sequences",
        return_value=(target_id_plus_seq_Q9NU22[1], conservation_id_plus_seq_Q9NU22_PF07728[1])
    ), patch(
        "transfer_annotations.populate_conservation",
        side_effect=mock_populate_conservation
    ), patch(
        "transfer_annotations.populate_go_data_for_annotations",
        side_effect=mock_populate_go_data
    ):
        result = cleanup_improve_transfer_dict(
            logger=logger,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict_populated_disulfid_Q9NU22,
            pfam_id="PF07728",
            hmmalign_lines=[],
            conservations_filepath="fake.json",
            annotations_filepath="fake.json",
            output_dir="fake_dir",
            pfam_interpro_map_filepath="fake.tsv"
        )
        assert result == transfer_dict_populated_disulfid_post_gos_Q9NU22

def test_cleanup_improve_transfer_dict_conservations_only(
    logger,
    multi_logger,
    transfer_dict_populated_disulfid_Q9NU22,
    conservations_content_Q9NU22_PF07728,
    mapping_content_Q9NU22_and_H0YB80_domains,
    target_id_plus_seq_Q9NU22,
    conservation_id_plus_seq_Q9NU22_PF07728,
    transfer_dict_populated_disulfid_post_conservation_Q9NU22
):
    """Test handling when only conservations data exists"""
    mapping_df = pd.read_csv(StringIO(mapping_content_Q9NU22_and_H0YB80_domains), sep="\t")

    with patch(
        "transfer_annotations.read_conservations_and_annotations",
        return_value=(conservations_content_Q9NU22_PF07728, {"sequence_id": {}})
    ), patch(
        "pandas.read_csv",
        return_value=mapping_df
    ), patch(
        "transfer_annotations.get_alignment_sequences",
        return_value=(target_id_plus_seq_Q9NU22[1], conservation_id_plus_seq_Q9NU22_PF07728[1])
    ), patch(
        "transfer_annotations.populate_conservation",
        side_effect=lambda transfer_dict, **kwargs: mock_populate_conservation(
            transfer_dict, transfer_dict_populated_disulfid_post_conservation_Q9NU22
        )
    ):
        result = cleanup_improve_transfer_dict(
            logger=logger,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict_populated_disulfid_Q9NU22,
            pfam_id="PF07728",
            hmmalign_lines=[],
            conservations_filepath="fake.json",
            annotations_filepath="nonexistent.json",
            output_dir="fake_dir",
            pfam_interpro_map_filepath="fake.tsv"
        )

        # Verify only conservation data was populated
        interval_path = result["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
        assert interval_path["conservations"] == transfer_dict_populated_disulfid_post_conservation_Q9NU22["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]["conservations"]
        assert "GO" not in interval_path["annotations"]["positions"]["333"][list(interval_path["annotations"]["positions"]["333"].keys())[0]]

def test_cleanup_improve_transfer_dict_handles_empty_conservations(
    logger,
    multi_logger,
    transfer_dict_populated_disulfid_Q9NU22,
    mapping_content_Q9NU22_and_H0YB80_domains,
    target_gos_Q9NU22_MCRB_ECOLI,
    target_id_plus_seq_Q9NU22,
    conservation_id_plus_seq_Q9NU22_PF07728,
):
    """Test handling of empty conservations data"""
    mapping_df = pd.read_csv(StringIO(mapping_content_Q9NU22_and_H0YB80_domains), sep="\t")

    with patch(
        "transfer_annotations.read_conservations_and_annotations",
        return_value=({}, {})  # Empty conservations and annotations
    ), patch(
        "pandas.read_csv",
        return_value=mapping_df
    ), patch(
        "transfer_annotations.gather_go_terms_for_target",
        return_value=target_gos_Q9NU22_MCRB_ECOLI
    ), patch(
        "transfer_annotations.get_alignment_sequences",
        return_value=(target_id_plus_seq_Q9NU22[1], conservation_id_plus_seq_Q9NU22_PF07728[1])
    ):
        result = cleanup_improve_transfer_dict(
            logger=logger,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict_populated_disulfid_Q9NU22,
            pfam_id="PF07728",
            hmmalign_lines=[],
            conservations_filepath="fake.json",
            annotations_filepath="fake.json",
            output_dir="fake_dir",
            pfam_interpro_map_filepath="fake.tsv"
        )

        # Verify structure is maintained but without conservation data
        interval = result["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
        assert interval["conservations"]["positions"] == {}
        assert interval["conservations"]["indices"]["matches"] == set()
        assert "GO" not in interval["annotations"]["positions"]["333"][list(interval["annotations"]["positions"]["333"].keys())[0]]

def test_cleanup_improve_transfer_dict_only_conservations_valid(
    logger,
    multi_logger,
    transfer_dict_populated_disulfid_Q9NU22,
    conservations_content_Q9NU22_PF07728,
    mapping_content_Q9NU22_and_H0YB80_domains,
    target_id_plus_seq_Q9NU22,
    conservation_id_plus_seq_Q9NU22_PF07728,
    transfer_dict_populated_disulfid_post_conservation_Q9NU22
):
    """Test handling when only conservations data is valid"""
    mapping_df = pd.read_csv(StringIO(mapping_content_Q9NU22_and_H0YB80_domains), sep="\t")
    expected_dict = copy.deepcopy(transfer_dict_populated_disulfid_post_conservation_Q9NU22)

    with patch(
        "transfer_annotations.read_conservations_and_annotations",
        return_value=(conservations_content_Q9NU22_PF07728, {"sequence_id": {}})
    ), patch(
        "pandas.read_csv",
        return_value=mapping_df
    ), patch(
        "transfer_annotations.get_alignment_sequences",
        return_value=(target_id_plus_seq_Q9NU22[1], conservation_id_plus_seq_Q9NU22_PF07728[1])
    ), patch(
        "transfer_annotations.populate_conservation",
        side_effect=lambda transfer_dict, **kwargs: mock_populate_conservation(
            transfer_dict, transfer_dict_populated_disulfid_post_conservation_Q9NU22
        )
    ):
        result = cleanup_improve_transfer_dict(
            logger=logger,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict_populated_disulfid_Q9NU22,
            pfam_id="PF07728",
            hmmalign_lines=[],
            conservations_filepath="fake.json",
            annotations_filepath="fake.json",
            output_dir="fake_dir",
            pfam_interpro_map_filepath="fake.tsv"
        )

        # Verify only conservation data was added
        interval_path = result["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
        assert interval_path["conservations"] == expected_dict["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]["conservations"]
        assert "GO" not in interval_path["annotations"]["positions"]["333"][list(interval_path["annotations"]["positions"]["333"].keys())[0]]

def test_cleanup_improve_transfer_dict_only_annotations_valid(
    logger,
    multi_logger,
    transfer_dict_populated_disulfid_Q9NU22,
    annotations_content_disulfid_fixture_Q9NU22_PF07728,
    mapping_content_Q9NU22_and_H0YB80_domains,
    target_gos_Q9NU22_MCRB_ECOLI,
    transfer_dict_populated_disulfid_post_gos_Q9NU22
):
    """Test handling when only annotations data is valid"""
    mapping_df = pd.read_csv(StringIO(mapping_content_Q9NU22_and_H0YB80_domains), sep="\t")
    expected_dict = copy.deepcopy(transfer_dict_populated_disulfid_post_gos_Q9NU22)

    with patch(
        "transfer_annotations.read_conservations_and_annotations",
        return_value=({"sequence_id/range": {}}, annotations_content_disulfid_fixture_Q9NU22_PF07728)
    ), patch(
        "pandas.read_csv",
        return_value=mapping_df
    ), patch(
        "transfer_annotations.gather_go_terms_for_target",
        return_value=target_gos_Q9NU22_MCRB_ECOLI
    ), patch(
        "transfer_annotations.populate_go_data_for_annotations",
        side_effect=lambda transfer_dict, **kwargs: mock_populate_go_data(
            transfer_dict, transfer_dict_populated_disulfid_post_gos_Q9NU22
        )
    ):
        result = cleanup_improve_transfer_dict(
            logger=logger,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict_populated_disulfid_Q9NU22,
            pfam_id="PF07728",
            hmmalign_lines=[],
            conservations_filepath="fake.json",
            annotations_filepath="fake.json",
            output_dir="fake_dir",
            pfam_interpro_map_filepath="fake.tsv"
        )

        # Verify only GO data was added
        interval_path = result["domain"]["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]
        assert interval_path["conservations"]["positions"] == {}
        assert interval_path["conservations"]["indices"]["matches"] == set()
        assert "GO" in interval_path["annotations"]["positions"]["333"][list(interval_path["annotations"]["positions"]["333"].keys())[0]]


###T convert_sets_and_tuples_to_lists

def test_convert_sets_and_tuples_to_lists_real_data(transfer_dict_populated_disulfid_Q9NU22):
    """Test conversion of real transfer dict data"""
    result = convert_sets_and_tuples_to_lists(transfer_dict_populated_disulfid_Q9NU22)

    # Check conversion of set in indices
    matches = result["DOMAIN"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]["annotations"]["indices"]["matches"]
    assert isinstance(matches, list)
    assert sorted(matches) == ["333", "373"]

    # Check conversion of tuples in annotation_ranges
    ranges = result["DOMAIN"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["hit_intervals"]["325-451"]["annotation_ranges"]["DISULFID | Intrachain (with C-246); in linked form"]["ranges"]
    assert isinstance(ranges[0], list)
    assert ranges == [[333, 333]]

def test_convert_sets_and_tuples_to_lists_simple():
    """Test conversion of simple set and tuple data"""
    test_data = {
        "set": {"a", "b", "c"},
        "tuple": (1, 2, 3)
    }

    result = convert_sets_and_tuples_to_lists(test_data)

    assert isinstance(result["set"], list)
    assert sorted(result["set"]) == ["a", "b", "c"]
    assert isinstance(result["tuple"], list)
    assert result["tuple"] == [1, 2, 3]

def test_convert_sets_and_tuples_to_lists_nested():
    """Test conversion of nested structures with sets and tuples"""
    test_data = {
        "nested": {
            "set": {"x", "y"},
            "data": [
                {"inner_set": {"1", "2"}},
                ("a", "b")
            ]
        }
    }

    result = convert_sets_and_tuples_to_lists(test_data)

    assert isinstance(result["nested"]["set"], list)
    assert sorted(result["nested"]["set"]) == ["x", "y"]
    assert isinstance(result["nested"]["data"][0]["inner_set"], list)
    assert sorted(result["nested"]["data"][0]["inner_set"]) == ["1", "2"]
    assert isinstance(result["nested"]["data"][1], list)
    assert result["nested"]["data"][1] == ["a", "b"]

def test_convert_sets_and_tuples_to_lists_annotation_ranges():
    """Test special handling of annotation_ranges structure"""
    test_data = {
        "positions": {1, 2, 3},
        "ranges": [(1, 2), (3, 4)]
    }

    result = convert_sets_and_tuples_to_lists(test_data)

    assert isinstance(result["positions"], list)
    assert result["positions"] == [1, 2, 3]
    assert isinstance(result["ranges"], list)
    assert isinstance(result["ranges"][0], list)
    assert isinstance(result["ranges"][1], list)
    assert result["ranges"] == [[1, 2], [3, 4]]

###T write_reports

def test_write_reports_with_data(
    tmp_path, logger, multi_logger,
    transfer_dict_populated_disulfid_post_conservation_Q9NU22
):
    """Test write_reports with actual annotation data"""
    mock_dict = transfer_dict_populated_disulfid_post_conservation_Q9NU22

    with patch("builtins.open", mock_open()) as mock_file:
        write_reports(logger, multi_logger, mock_dict, str(tmp_path))

        # Verify correct files were created
        calls = mock_file.call_args_list
        assert len(calls) == 2

        # Verify paths
        expected_paths = [
            os.path.join(str(tmp_path), "PF07728", "PF07728_report.json"),
            os.path.join(str(tmp_path), "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json")
        ]
        actual_paths = [call[0][0] for call in calls]
        assert sorted(actual_paths) == sorted(expected_paths)

        # Verify debug logs were called
        assert logger.debug.call_count == 3

def test_write_reports_empty_data(tmp_path, logger, multi_logger):
    """Test write_reports with empty transfer dict"""
    transfer_dict = {
        "domain": {
            "PF07728": {
                "sequence_id": {
                    "sp|Q9NU22|MDN1_HUMAN": {}
                }
            }
        }
    }

    with patch("builtins.open", mock_open()) as mock_file:
        write_reports(logger, multi_logger, transfer_dict, str(tmp_path))

        calls = mock_file.call_args_list
        assert len(calls) == 2  # Main report + target report
        assert logger.debug.call_count == 3

def test_write_reports_multiple_targets(tmp_path, logger, multi_logger):
    """Test write_reports with multiple target sequences"""
    transfer_dict = {
        "domain": {
            "PF07728": {
                "sequence_id": {
                    "TEST1": {"some_data": "value1"},
                    "TEST2": {"some_data": "value2"}
                }
            }
        }
    }

    with patch("builtins.open", mock_open()) as mock_file:
        write_reports(logger, multi_logger, transfer_dict, str(tmp_path))

        calls = mock_file.call_args_list
        assert len(calls) == 3  # Main report + 2 target reports
        assert logger.debug.call_count == 5  # 1 main + 2 data writes + 2 confirmations

def test_write_reports_pipe_in_target_name(tmp_path, logger, multi_logger):
    """Test write_reports handles target names with pipes"""
    transfer_dict = {
        "domain": {
            "PF07728": {
                "sequence_id": {
                    "TEST|A": {"some_data": "value"}
                }
            }
        }
    }

    with patch("builtins.open", mock_open()) as mock_file:
        write_reports(logger, multi_logger, transfer_dict, str(tmp_path))

        calls = mock_file.call_args_list
        assert len(calls) == 2  # Main report + target report

        # Check paths
        paths = [call[0][0] for call in calls]
        assert any("TEST-A" in path for path in paths)
        assert logger.debug.call_count == 3

def test_write_reports_completely_empty_dict(tmp_path, logger, multi_logger):
    """Test write_reports with completely empty transfer dict"""
    transfer_dict = {"domain" : {}}

    write_reports(logger, multi_logger, transfer_dict, str(tmp_path))
    multi_logger.assert_called_once_with(
        "warning",
        "TRANSFER_ANNOTS --- WRITE_REPORT --- Transfer dictionary is empty - no value assigned to 'domain' - skipping report writing")

def test_write_reports_empty_sequence_id(tmp_path, logger, multi_logger):
    """Test write_reports with empty sequence_id dict"""
    transfer_dict = {
        "domain": {
            "PF07728": {
                "sequence_id": {}
            }
        }
    }

    with patch("builtins.open", mock_open()) as mock_file:
        write_reports(logger, multi_logger, transfer_dict, str(tmp_path))

        calls = mock_file.call_args_list
        assert len(calls) == 1  # Only main report
        assert logger.debug.call_count == 1

def test_write_reports_none_target_data(tmp_path, logger, multi_logger):
    """Test write_reports with None as target data"""
    transfer_dict = {
        "domain": {
            "PF07728": {
                "sequence_id": {
                    "TEST": None
                }
            }
        }
    }

    with patch("json.dump") as mock_json_dump:
        write_reports(logger, multi_logger, transfer_dict, str(tmp_path))

        # Verify number of calls
        assert mock_json_dump.call_count == 2
        assert logger.debug.call_count == 3

        # Get the data passed to json.dump calls
        dump_calls = mock_json_dump.call_args_list

        # Check main report structure (first call)
        main_report_data = dump_calls[0][0][0]  # First arg of first call
        expected_main_report = {
            "domain_id": "PF07728",
            "sequences": {"TEST": None}
        }
        assert main_report_data == expected_main_report

        # Check target report structure (second call)
        target_report_data = dump_calls[1][0][0]  # First arg of second call
        expected_target_report = {
            "sequence_id": "TEST",
            "domain": {
                "PF07728": {"annotations": "None"}
            }
        }
        assert target_report_data == expected_target_report


###T map_and_filter_annot_pos

def test_map_and_filter_with_paired_positions(
    logger, multi_logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_disulfid_pair, transfer_dict,
    get_annotation_dict, annotations_content_disulfid_fixture_Q9NU22_PF07728
):
    """Unit test for map_and_filter_annot_pos with paired positions."""
    annotation_dict = get_annotation_dict(annotations_content_disulfid_fixture_Q9NU22_PF07728, "246")

    with patch("transfer_annotations.validate_paired_annotations") as mock_validate_paired_annotations:
        # Configure the mock to return expected values
        mock_validate_paired_annotations.return_value = (True, {"dummy": "result"})

        # Call the function under test
        result = map_and_filter_annot_pos(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22,
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
            paired_annotation_dict=annotation_dict
        )

        # Assert that validate_paired_annotations was called with expected arguments
        mock_validate_paired_annotations.assert_called_once_with(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            paired_annotation_dict=annotation_dict,
            entry_annotations=entry_annotations_disulfid_pair,
            counter_target_pos=ANY,  # Since this is set inside the function
            counter_annot_pos=ANY,  # Since this is set inside the function
            paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22,
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
            processed_annotations=set()
        )

        # Verify the result
        assert result == (True, {"dummy": "result"})

def test_map_and_filter_without_paired_positions(
    logger, multi_logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, transfer_dict
):
    """Unit test for map_and_filter_annot_pos without paired positions."""

    with patch("transfer_annotations.validate_annotations") as mock_validate_annotations:
        # Call the function under test
        result = map_and_filter_annot_pos(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            # No paired position parameters
        )

        # Assert that validate_annotations was called with expected arguments
        mock_validate_annotations.assert_called_once_with(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            counter_target_pos=ANY,  # These are set within the function
            counter_annot_pos=ANY,
        )

        # Since the function doesn't return anything in this case
        assert result is None

###T add_to_transfer_dict

def test_add_to_transfer_dict_initializes_structure(
    multi_logger, transfer_dict, target_sequence_continuous_Q9NU22,
    anno_total_disulfid_MCRB_ECOLI_Q9NU22_205
):
    """Test dict initialization with minimal data"""
    anno_id_mock = anno_total_disulfid_MCRB_ECOLI_Q9NU22_205["type"] + " | " + anno_total_disulfid_MCRB_ECOLI_Q9NU22_205["description"]
    with patch("transfer_annotations._add_single_annotation") as mock_add_single:
        add_to_transfer_dict(
            hit=True,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_continuous_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id=anno_id_mock,
            anno_total=anno_total_disulfid_MCRB_ECOLI_Q9NU22_205,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22
        )

        # Check structure initialization
        assert "DOMAIN" in transfer_dict
        assert "sequence_id" in transfer_dict["DOMAIN"]
        assert target_name_mock_Q9NU22 in transfer_dict["DOMAIN"]["sequence_id"]

        interval_dict = transfer_dict["DOMAIN"]["sequence_id"][target_name_mock_Q9NU22]["hit_intervals"]["325-451"]
        assert interval_dict["sequence"] == target_sequence_continuous_Q9NU22
        assert interval_dict["length"] == len(target_sequence_continuous_Q9NU22)
        assert interval_dict["hit_start"] == target_hit_start_mock_Q9NU22
        assert interval_dict["hit_end"] == target_hit_end_mock_Q9NU22
        assert interval_dict["annotations"] == {"positions": {}, "indices": {"matches": set(), "misses": set()}}
        assert interval_dict["conservations"] == {"positions": {}, "indices": {"matches": set(), "misses": set()}}
        assert interval_dict["position_conversion"] == {"target_to_aln": {}, "aln_to_target": {}}
        assert interval_dict["annotation_ranges"] == {}

        # Verify expected call to helper
        mock_add_single.assert_called_once()

def test_add_to_transfer_dict_handles_paired_annotations(
    multi_logger, transfer_dict, target_sequence_continuous_Q9NU22,
    anno_total_disulfid_MCRB_ECOLI_Q9NU22_205, anno_total_disulfid_MCRB_ECOLI_Q9NU22_246
):
    """Test proper handling of paired annotations"""
    with patch("transfer_annotations._add_single_annotation") as mock_add_single:
        add_to_transfer_dict(
            hit=True,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_continuous_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id="DISULFID | Intrachain (with C-246); in linked form",
            anno_total=anno_total_disulfid_MCRB_ECOLI_Q9NU22_205,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            # Paired annotation parameters
            paired_position_res_hit=True,
            paired_anno_id="DISULFID | Intrachain (with C-205); in linked form",
            paired_anno_total=anno_total_disulfid_MCRB_ECOLI_Q9NU22_246
        )

        # Verify both annotations processed
        assert mock_add_single.call_count == 2
        assert mock_add_single.call_args_list[0][1]["anno_id"] == "DISULFID | Intrachain (with C-246); in linked form"
        assert mock_add_single.call_args_list[1][1]["anno_id"] == "DISULFID | Intrachain (with C-205); in linked form"

def test_add_to_transfer_dict_paired_repeated_anno_id(
    multi_logger, transfer_dict_populated_disulfid_Q9NU22,
    target_sequence_continuous_Q9NU22,
    mock_make_anno_total_disulfid_return_205_Q9NU22,
    mock_make_anno_total_disulfid_return_246_Q9NU22,
):
    """Test paired disulfide annotation processing"""
    anno_id = mock_make_anno_total_disulfid_return_205_Q9NU22["anno_id"]
    anno_total = mock_make_anno_total_disulfid_return_205_Q9NU22["anno_total"]
    paired_anno_id = mock_make_anno_total_disulfid_return_246_Q9NU22["anno_id"]
    paired_anno_total = mock_make_anno_total_disulfid_return_246_Q9NU22["anno_total"]

    interval_key = f"{target_hit_start_mock_Q9NU22}-{target_hit_end_mock_Q9NU22}"
    positions_dict = transfer_dict_populated_disulfid_Q9NU22["DOMAIN"]["sequence_id"][target_name_mock_Q9NU22]["hit_intervals"][interval_key]["annotations"]["positions"]

    add_to_transfer_dict(
        hit=True,
        multi_logger=multi_logger,
        transfer_dict=transfer_dict_populated_disulfid_Q9NU22,
        target_name=target_name_mock_Q9NU22,
        target_sequence_continuous=target_sequence_continuous_Q9NU22,
        target_hit_start=target_hit_start_mock_Q9NU22,
        target_hit_end=target_hit_end_mock_Q9NU22,
        anno_id=anno_id,
        anno_total=anno_total,
        entry_mnemo_name=entry_mnemo_name_mock_repeated_Q9NU22,
        entry_primary_accession=entry_primary_accession_mock_repeated_Q9NU22,
        paired_position_res_hit=True,
        paired_anno_id=paired_anno_id,
        paired_anno_total=paired_anno_total
    )

    # Helper function to access nested dictionary paths
    def get_annotation(target_position_str, positions_dict):
        return positions_dict[target_position_str]

    # Test data
    positions = {
        "333": {"partner": "373", "annot_pos": "205"},
        "373": {"partner": "333", "annot_pos": "246"}
    }

    for target_position_str, data in positions.items():
        anno = get_annotation(target_position_str, positions_dict)["DISULFID | Intrachain (with C-246); in linked form" if target_position_str == "333" else "DISULFID | Intrachain (with C-205); in linked form"]

        # Count assertions
        assert anno["essentials"]["count"] == 2
        assert anno["paired_position"][data["partner"]]["count"] == 2
        assert anno["additional_keys"]["annot_position"][data["annot_pos"]]["count"] == 2

        # Evidence assertions
        for eco_code in ["ECO:0000250", "ECO:0000269|PubMed:12345678"]:
            assert eco_code in anno["evidence"]
            if eco_code == "ECO:0000250":
                evidence = anno["evidence"][eco_code]
                assert evidence["rep_mnemo_name"] == "TPA_HUMAN"
                assert evidence["rep_primary_accession"] == "P00720"
                assert evidence["count"] == 1

    # Verify no misses
    # print(json.dumps(transfer_dict_populated_disulfid_Q9NU22, indent=4))
    assert positions_dict["333"][anno_id]["hit"] is True
    assert positions_dict["373"][paired_anno_id]["hit"] is True

def test_add_to_transfer_dict_hit_miss_pair(
    multi_logger, transfer_dict,
    target_sequence_continuous_Q9NU22,
    mock_make_anno_total_disulfid_return_205_P15005,
    mock_make_anno_total_disulfid_return_246_P15005,
    anno_total_disulfid_MCRB_ECOLI_Q9NU22_205,
    anno_total_disulfid_MCRB_ECOLI_Q9NU22_246
):
    """Test disulfide pair where one residue hits and the other misses"""
    anno_id = mock_make_anno_total_disulfid_return_205_P15005["anno_id"]
    anno_total = anno_total_disulfid_MCRB_ECOLI_Q9NU22_205.copy()
    paired_anno_total = anno_total_disulfid_MCRB_ECOLI_Q9NU22_246.copy()
    paired_anno_id = mock_make_anno_total_disulfid_return_246_P15005["anno_id"]
    hit_interval = f"{target_hit_start_mock_Q9NU22}-{target_hit_end_mock_Q9NU22}"

    add_to_transfer_dict(
        hit=True,
        multi_logger=multi_logger,
        transfer_dict=transfer_dict,
        target_name=target_name_mock_Q9NU22,
        target_sequence_continuous=target_sequence_continuous_Q9NU22,
        target_hit_start=target_hit_start_mock_Q9NU22,
        target_hit_end=target_hit_end_mock_Q9NU22,
        anno_id=anno_id,
        anno_total=anno_total,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
        paired_position_res_hit=False,
        paired_anno_id=paired_anno_id,
        paired_anno_total=paired_anno_total
    )

    positions_dict = transfer_dict["DOMAIN"]["sequence_id"][target_name_mock_Q9NU22]["hit_intervals"][hit_interval]["annotations"]["positions"]

    # Assert hit position
    assert "333" in positions_dict
    assert positions_dict["333"][anno_id]["hit"] is True

    # Assert missed position
    assert "373" in positions_dict
    assert positions_dict["373"][paired_anno_id]["hit"] is False


def test_add_to_transfer_dict_disulfid_single_mocked_helper(
    multi_logger, transfer_dict,
    target_sequence_continuous_Q9NU22,
    mock_make_anno_total_disulfid_return_205_P15005,
    anno_total_disulfid_MCRB_ECOLI_Q9NU22_205
):
    """Test add_to_transfer_dict with successful annotation addition"""

    anno_id = mock_make_anno_total_disulfid_return_205_P15005["anno_id"]
    anno_total = anno_total_disulfid_MCRB_ECOLI_Q9NU22_205
    interval_key = f"{target_hit_start_mock_Q9NU22}-{target_hit_end_mock_Q9NU22}"

    # Call the function
    with patch("transfer_annotations._add_single_annotation") as mock_add_single:
        add_to_transfer_dict(
            hit=True,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_continuous_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id=anno_id,
            anno_total=anno_total,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
        )

        # Assertions to verify correct addition to transfer_dict
        mock_add_single.assert_called_once_with(
            hit=True,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            interval_key=interval_key,
            anno_id=anno_id,
            anno_total=anno_total,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            additional_keys={}
        )

def test_add_to_transfer_dict_disulfid_paired_mocked_helper(
    multi_logger, transfer_dict,
    target_sequence_continuous_Q9NU22,
    mock_make_anno_total_disulfid_return_205_P15005,
    mock_make_anno_total_disulfid_return_246_P15005,
    anno_total_disulfid_MCRB_ECOLI_Q9NU22_205,
    anno_total_disulfid_MCRB_ECOLI_Q9NU22_246
):
    """
    Test add_to_transfer_dict with successful annotation pair
    addition. Intent on seeing that expected calls went through.
    """

    anno_id = mock_make_anno_total_disulfid_return_205_P15005["anno_id"]
    anno_total = anno_total_disulfid_MCRB_ECOLI_Q9NU22_205.copy()
    paired_anno_total = anno_total_disulfid_MCRB_ECOLI_Q9NU22_246.copy()
    paired_anno_id = mock_make_anno_total_disulfid_return_246_P15005["anno_id"]
    interval_key = f"{target_hit_start_mock_Q9NU22}-{target_hit_end_mock_Q9NU22}"

    # Call the function
    with patch("transfer_annotations._add_single_annotation") as mock_add_single:
        add_to_transfer_dict(
            hit=True,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_continuous_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id=anno_id,
            anno_total=anno_total,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            paired_position_res_hit=True,
            paired_anno_id=paired_anno_id,
            paired_anno_total=paired_anno_total
        )

        expected_calls = [
            call(
                hit=True,
                transfer_dict=transfer_dict,
                target_name=target_name_mock_Q9NU22,
                interval_key=interval_key,
                anno_id=anno_id,
                anno_total=anno_total,
                entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
                entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
                additional_keys={}
            ),
            call(
                hit=True,
                transfer_dict=transfer_dict,
                target_name=target_name_mock_Q9NU22,
                interval_key=interval_key,
                anno_id=paired_anno_id,
                anno_total=paired_anno_total,
                entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
                entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
                additional_keys={}
            )
        ]

        mock_add_single.assert_has_calls(expected_calls, any_order=False)
        assert mock_add_single.call_count == 2

###T _add_single_annotation

def test_add_single_annotation_complete(
    logger, transfer_dict_initialized_structure_Q9NU22, anno_total_disulfid_MCRB_ECOLI_Q9NU22_205,
):
    """Test all aspects of _add_single_annotation data population"""
    with patch("transfer_annotations.update_position_ranges"):
        interval_dict = transfer_dict_initialized_structure_Q9NU22["DOMAIN"]["sequence_id"][target_name_mock_Q9NU22]["hit_intervals"]["325-451"]

        _add_single_annotation(
            hit=True,
            transfer_dict=transfer_dict_initialized_structure_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            interval_key="325-451",
            anno_id="DISULFID | Intrachain (with C-246); in linked form",
            anno_total=anno_total_disulfid_MCRB_ECOLI_Q9NU22_205,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            additional_keys={"annot_position": "205"}
        )

        position_data = interval_dict["annotations"]["positions"]["333"]["DISULFID | Intrachain (with C-246); in linked form"]

        # Verify essentials
        assert position_data["essentials"]["type"] == "DISULFID"
        assert position_data["essentials"]["count"] == 1
        assert position_data["hit"] is True

        # Verify evidence
        assert "ECO:0000269|PubMed:12345678" in position_data["evidence"]
        evidence = position_data["evidence"]["ECO:0000269|PubMed:12345678"]
        assert evidence["rep_primary_accession"] == entry_primary_accesion_mock_Q9NU22
        assert evidence["count"] == 1

        # Verify paired position
        assert "373" in position_data["paired_position"]
        paired = position_data["paired_position"]["373"]
        assert paired["rep_primary_accession"] == entry_primary_accesion_mock_Q9NU22
        assert paired["count"] == 1
        assert "gapped_paired" not in paired
        assert "insert_column_paired" not in paired

        # Verify additional keys
        assert "annot_position" in position_data["additional_keys"]
        add_key = position_data["additional_keys"]["annot_position"]["205"]
        assert add_key["rep_primary_accession"] == entry_primary_accesion_mock_Q9NU22
        assert add_key["count"] == 1

        # Verify indices tracking
        assert "333" in interval_dict["annotations"]["indices"]["matches"]
        assert interval_dict["position_conversion"]["target_to_aln"]["333"] == "18"
        assert interval_dict["position_conversion"]["aln_to_target"]["18"] == "333"

def test_add_single_annotation_gapped_paired(
    logger, transfer_dict_initialized_structure_Q9NU22
):
    """Test _add_single_annotation handling of gapped paired position"""
    with patch("transfer_annotations.update_position_ranges"):
        # Modify anno_total to include gapped_paired flag
        anno_total = {
            "type": "DISULFID",
            "description": "Intrachain (with C-246); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "333",
            "annot_position": "205",
            "paired_target_position": "373",
            "index_position": "18",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
            "gapped_paired": True  # Add the flag
        }

        _add_single_annotation(
            hit=True,
            transfer_dict=transfer_dict_initialized_structure_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            interval_key="325-451",
            anno_id="DISULFID | Intrachain (with C-246); in linked form",
            anno_total=anno_total,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            additional_keys={}
        )

        position_data = transfer_dict_initialized_structure_Q9NU22["DOMAIN"]["sequence_id"][target_name_mock_Q9NU22]["hit_intervals"]["325-451"]["annotations"]["positions"]["333"]["DISULFID | Intrachain (with C-246); in linked form"]

        # Verify paired position has gapped flag
        assert "373" in position_data["paired_position"]
        paired = position_data["paired_position"]["373"]
        assert paired["gapped_target"] is True
        assert "insert_column_paired" not in paired

def test_add_single_annotation_insert_column_paired(
    logger, transfer_dict_initialized_structure_Q9NU22
):
    """Test _add_single_annotation handling of insert column paired position"""
    with patch("transfer_annotations.update_position_ranges"):
        # Modify anno_total to include insert_column_paired flag
        anno_total = {
            "type": "DISULFID",
            "description": "Intrachain (with C-246); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "333",
            "annot_position": "205",
            "paired_target_position": "373",
            "index_position": "18",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
            "insert_column_paired": True  # Add the flag
        }

        _add_single_annotation(
            hit=True,
            transfer_dict=transfer_dict_initialized_structure_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            interval_key="325-451",
            anno_id="DISULFID | Intrachain (with C-246); in linked form",
            anno_total=anno_total,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            additional_keys={}
        )

        position_data = transfer_dict_initialized_structure_Q9NU22["DOMAIN"]["sequence_id"][target_name_mock_Q9NU22]["hit_intervals"]["325-451"]["annotations"]["positions"]["333"]["DISULFID | Intrachain (with C-246); in linked form"]

        # Verify paired position has insert column flag
        assert "373" in position_data["paired_position"]
        paired = position_data["paired_position"]["373"]
        assert paired["insert_column"] is True
        assert "gapped_target" not in paired


###T update_position_range

def test_update_position_ranges_new_anno(interval_dict_empty):
    """Test creation of new annotation range"""
    with patch("transfer_annotations.merge_adjacent_ranges") as mock_merge:
        update_position_ranges(
            interval_dict=interval_dict_empty,
            target_position_str="333",
            ranges_dict_key="annotation_ranges",
            range_id="DISULFID | Intrachain (with C-246); in linked form"
            )

        ranges = interval_dict_empty["annotation_ranges"]["DISULFID | Intrachain (with C-246); in linked form"]
        assert ranges["ranges"] == [(333, 333)]
        assert ranges["positions"] == {333}
        mock_merge.assert_not_called()

def test_update_position_ranges_extend_integration(interval_dict_with_range):
    """Test extending existing range with integration"""
    update_position_ranges(
        interval_dict=interval_dict_with_range,
        target_position_str="334",
        ranges_dict_key="annotation_ranges",
        range_id="DISULFID | Test"
        )

    ranges = interval_dict_with_range["annotation_ranges"]["DISULFID | Test"]
    assert ranges["ranges"] == [(333, 334)]
    assert ranges["positions"] == {333, 334}

def test_update_position_ranges_position_in_range(interval_dict_with_position_in_range):
    """Test adding position that falls within existing range"""
    update_position_ranges(
        interval_dict=interval_dict_with_position_in_range,
        target_position_str="334",
        ranges_dict_key="annotation_ranges",
        range_id="DISULFID | Test"
        )

    ranges = interval_dict_with_position_in_range["annotation_ranges"]["DISULFID | Test"]
    assert ranges["ranges"] == [(333, 335)]  # Range shouldn't change
    assert ranges["positions"] == {333, 334, 335}  # Position already in set

def test_update_position_ranges_merge_multiple(interval_dict_with_multiple_ranges):
    """Test merging of multiple ranges when adding connecting position"""
    update_position_ranges(
        interval_dict=interval_dict_with_multiple_ranges,
        target_position_str="335",
        ranges_dict_key="annotation_ranges",
        range_id="DISULFID | Test"
        )

    ranges = interval_dict_with_multiple_ranges["annotation_ranges"]["DISULFID | Test"]
    assert ranges["ranges"] == [(333, 337)]  # All ranges merged
    assert ranges["positions"] == {333, 334, 335, 336, 337}

def test_update_position_ranges_non_adjacent(interval_dict_with_multiple_ranges):
    """Test adding non-adjacent position"""
    update_position_ranges(
        interval_dict=interval_dict_with_multiple_ranges,
        target_position_str="340",
        ranges_dict_key="annotation_ranges",
        range_id="DISULFID | Test"
        )

    ranges = interval_dict_with_multiple_ranges["annotation_ranges"]["DISULFID | Test"]
    assert ranges["ranges"] == [(333, 334), (336, 337), (340, 340)]
    assert ranges["positions"] == {333, 334, 336, 337, 340}

def test_update_position_ranges_duplicate_pos(interval_dict_with_range):
    """Test adding duplicate position"""
    with patch("transfer_annotations.merge_adjacent_ranges") as mock_merge:
        update_position_ranges(
            interval_dict=interval_dict_with_range,
            target_position_str="333",
            ranges_dict_key="annotation_ranges",
            range_id="DISULFID | Test"
            )

        ranges = interval_dict_with_range["annotation_ranges"]["DISULFID | Test"]
        assert ranges["ranges"] == [(333, 333)]
        assert ranges["positions"] == {333}
        mock_merge.assert_not_called()


###T merge_adjacent_ranges

def test_merge_adjacent_ranges():
    """Test merging of adjacent ranges"""
    ranges = [(333, 334), (335, 337), (338, 340)]
    merge_adjacent_ranges(ranges)
    assert ranges == [(333, 340)]

def test_merge_overlapping_ranges():
    """Test merging of overlapping ranges"""
    ranges = [(333, 335), (334, 337), (336, 340)]
    merge_adjacent_ranges(ranges)
    assert ranges == [(333, 340)]

def test_merge_non_adjacent_ranges():
    """Test preservation of non-adjacent ranges"""
    ranges = [(333, 334), (336, 337), (340, 341)]
    merge_adjacent_ranges(ranges)
    assert ranges == [(333, 334), (336, 337), (340, 341)]

def test_merge_adjacent_ranges_empty():
    """Test merging with empty ranges list"""
    ranges = []
    merge_adjacent_ranges(ranges)
    assert ranges == []

def test_merge_adjacent_ranges_single():
    """Test merging with single range"""
    ranges = [(333, 334)]
    merge_adjacent_ranges(ranges)
    assert ranges == [(333, 334)]

###T make_anno_total_dict

def test_make_anno_total_dict_basic(
    logger, good_eco_codes_all, entry_annotations_disulfid_pair,
    annotation_dict_205_Q9NU22
):
    """Test basic case with valid DISULFID annotation"""
    result = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_dict_205_Q9NU22,
        counter_target_pos_str="333",
        counter_annot_pos_str="205",
        index=18,
        target_amino="C",
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair,
    )

    expected = {
        "annotation": {
            **annotation_dict_205_Q9NU22,
        },
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-246); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-246); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "333",
            "annot_position": "205",
            "index_position": "18",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
        },
        "paired_annot_pos_str": "246"
    }

    assert result == expected

def test_make_anno_total_dict_multiple_evidence(
    logger, good_eco_codes_all, entry_annotations_disulfid_pair,
    annotation_dict_205_Q9NU22
):
    """Test case with multiple evidences in the string"""
    annotation_dict_205_Q9NU22["evidence"] = ["ECO:0000269|PubMed:12345678, ECO:0007744|PDB:1ABC"]

    result_single_pass = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_dict_205_Q9NU22,
        counter_target_pos_str="333",
        counter_annot_pos_str="205",
        index=18,
        target_amino="C",
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair
    )

    expected_single_pass = {
        "annotation": {
            **annotation_dict_205_Q9NU22,
        },
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-246); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-246); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "333",
            "annot_position": "205",
            "index_position": "18",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
        },
        "paired_annot_pos_str": "246"
    }

    expected_multiple_pass = copy.deepcopy(expected_single_pass)
    expected_multiple_pass["anno_total"]["evidence"] = "ECO:0000269|PubMed:12345678, ECO:0007744|PDB:1ABC"

    assert result_single_pass == expected_single_pass

    good_eco_codes_all_modded = good_eco_codes_all + ["ECO:0007744"]

    result_multiple_pass = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all_modded,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_dict_205_Q9NU22,
        counter_target_pos_str="333",
        counter_annot_pos_str="205",
        index=18,
        target_amino="C",
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair
    )

    assert result_multiple_pass == expected_multiple_pass

def test_make_anno_total_dict_eco_filtering(
    logger, entry_annotations_binding_only,
    get_annotation_dict, annotations_content_binding_fixture_Q9NU22_PF07728
):
    """Test ECO code filtering with BINDING annotation"""
    annotation_dict = get_annotation_dict(annotations_content_binding_fixture_Q9NU22_PF07728, "201")

    result = make_anno_total_dict(
        good_eco_codes=["ECO:0000269"],  # Exclude all other codes, including current: ECO:0000255
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_dict,
        counter_target_pos_str="329",
        counter_annot_pos_str="201",
        index=14,
        target_amino="G",
        logger=logger,
        entry_annotations=entry_annotations_binding_only
    )

    assert result["anno_total"] is None

def test_make_anno_total_dict_binding_type(
    logger, good_eco_codes_all, entry_annotations_binding_only,
    get_annotation_dict, annotations_content_binding_fixture_Q9NU22_PF07728
):
    """Test BINDING type annotation with additional keys"""
    annotation_dict = get_annotation_dict(annotations_content_binding_fixture_Q9NU22_PF07728, "201")

    result = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_dict,
        counter_target_pos_str="329",
        counter_annot_pos_str="201",
        index=14,
        target_amino="G",
        logger=logger,
        entry_annotations=entry_annotations_binding_only
    )

    expected = {
        "annotation": {
            **annotation_dict,
        },
        "anno_type": "BINDING",
        "anno_id": "BINDING | Interacts with GTP",
        "anno_total": {
            "type": "BINDING",
            "description": "Interacts with GTP",
            "count": 1,
            "evidence": "ECO:0000255",
            "ligand_id": "ChEBI:CHEBI:37565",
            "target_position": "329",
            "annot_position": "201",
            "index_position": "14",
            "annot_amino_acid": "G",
            "target_amino_acid": "G"
        },
        "paired_annot_pos_str": None
    }

    assert result["annotation"] is annotation_dict
    assert result["anno_type"] == "BINDING"
    assert result["anno_id"] == "BINDING | Interacts with GTP"
    assert result["paired_annot_pos_str"] is None
    assert result == expected

def test_make_anno_total_dict_with_caller_target_pos(
    logger, good_eco_codes_all, entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
):
    """Test with caller_target_pos provided (paired annotation case)"""
    result = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_dict_246_Q9NU22,
        counter_target_pos_str="373",
        counter_annot_pos_str="246",
        caller_target_pos_str="333",
        index=18,
        target_amino="C",
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair
    )

    expected = {
        "annotation": {
            **annotation_dict_246_Q9NU22,
        },
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-205); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "373",
            "annot_position": "246",
            "index_position": "18",
            "annot_amino_acid": "C",
            "target_amino_acid": "C",
            "paired_target_position": "333",
        },
        "paired_annot_pos_str": "205"
    }

    assert result["annotation"] is annotation_dict_246_Q9NU22
    assert result["anno_type"] == "DISULFID"
    assert result["anno_id"] == "DISULFID | Intrachain (with C-205); in linked form"
    assert result["paired_annot_pos_str"] == "205"
    assert result == expected

def test_make_anno_total_dict_no_evidence(
    logger, good_eco_codes_all
):
    """Test case where annotation has no evidence"""
    annotation_without_evidence = {
        "type": "BINDING",
        "description": "Test binding",
        "entry": "P12345"
    }

    result = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_without_evidence,
        counter_target_pos_str="100",
        counter_annot_pos_str="201",
        index=14,
        target_amino="G",
        logger=logger,
        entry_annotations={"201": [annotation_without_evidence]}
    )

    assert result["annotation"] is annotation_without_evidence
    assert result["anno_type"] is "BINDING"
    assert result["anno_id"] == "BINDING | Test binding"
    assert result["paired_annot_pos_str"] is None
    assert result["anno_total"] is None

###T process_annotation

def test_process_annotation_basic_single(
    logger, multi_logger, good_eco_codes_all, transfer_dict,
    target_sequence_Q9NU22, target_sequence_continuous_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, get_annotation_dict,
    annotations_content_binding_fixture_Q9NU22_PF07728
):
    """Test processing of a single BINDING annotation"""
    annotation_dict = get_annotation_dict(annotations_content_binding_fixture_Q9NU22_PF07728, "201")
    mock_make_anno_total = {
        "annotation": annotation_dict,
        "anno_type": "BINDING",
        "anno_id": "BINDING | Interacts with GTP",
        "anno_total": {
            "type": "BINDING",
            "description": "Interacts with GTP",
            "count": 1,
            "evidence": "ECO:0000255",
            "ligand_id": "ChEBI:CHEBI:37565",
            "target_position": "329",
            "annot_position": "201",
            "index_position": "14",
            "annot_amino_acid": "G",
            "target_amino_acid": "G",
        },
        "paired_annot_pos_str": None
    }
    processed_annotations = set()

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_make_anno_total) as mock_make_total, \
         patch("transfer_annotations.add_to_transfer_dict") as mock_add_dict:
        process_annotation(
            res_hit=True,
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            annotation_dict=annotation_dict,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            counter_target_pos_str="329",
            counter_annot_pos_str="201",
            index=14,
            target_amino="G",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "201", "329", "BINDING")
        )

        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict,
            counter_target_pos_str="329",
            counter_annot_pos_str="201",
            index=14,
            target_amino="G",
            logger=logger,
            entry_annotations=entry_annotations_binding_only
        )

        mock_add_dict.assert_called_once_with(
            hit=True,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_continuous_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id="BINDING | Interacts with GTP",
            anno_total=mock_make_anno_total["anno_total"],
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
        )

        assert ("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "201", "329", "BINDING") in processed_annotations

def test_process_annotation_paired_failure(
    logger, multi_logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair,
    target_sequence_Q9NU22, target_sequence_continuous_Q9NU22, annot_sequence_Q9NU22,
    mock_make_anno_total_disulfid_return_205_P15005,
    mock_map_filter_disulfid_return,
    annotation_dict_205_Q9NU22, annotation_dict_246_Q9NU22
):
    """
    Test processing of paired DISULFID annotations
    where paired pair fails to match, different residues
    """

    # Modify target sequence to ensure failure at paired position
    target_sequence_mod_to_fail_73 = list(target_sequence_Q9NU22)
    # Change C to R in MYCCT substring, note 0-based index
    target_sequence_mod_to_fail_73[72] = "R"
    target_sequence_mod_to_fail_73 = "".join(target_sequence_mod_to_fail_73)

    target_sequence_mod_to_fail_73_continuous = list(target_sequence_continuous_Q9NU22)
    # Index accomodates removal of gaps and insertions
    target_sequence_mod_to_fail_73_continuous[48] = "R"
    target_sequence_mod_to_fail_73_continuous = "".join(target_sequence_mod_to_fail_73_continuous)

    _, mock_paired_result_dict  = mock_map_filter_disulfid_return
    mock_map_filter_disulfid_fail_return = (False, mock_paired_result_dict)

    processed_annotations = set()

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch("transfer_annotations.map_and_filter_annot_pos", return_value=mock_map_filter_disulfid_fail_return) as mock_map_filter, \
         patch("transfer_annotations.add_to_transfer_dict") as mock_add_transfer:
        process_annotation(
            res_hit=True,
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence_mod_to_fail_73,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID", "match_column")
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair
        )

        # Verify map_and_filter_annot_pos was called correctly
        mock_map_filter.assert_called_once_with(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_mod_to_fail_73,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            processed_annotations=processed_annotations,
            paired_annot_pos_str="246",
            caller_target_pos_str="333",
            paired_annotation_dict=annotation_dict_246_Q9NU22
        )

        # Verify remove_failed_annotations was called correctly
        mock_add_transfer.assert_called_once_with(
            hit=True,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_mod_to_fail_73_continuous,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id="DISULFID | Intrachain (with C-246); in linked form",
            anno_total=mock_make_anno_total_disulfid_return_205_P15005["anno_total"],
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            paired_position_res_hit=False,
            paired_anno_id="DISULFID | Intrachain (with C-205); in linked form",
            paired_anno_total=mock_map_filter_disulfid_return[1]["anno_total"]
        )

        assert ("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID", "match_column") in processed_annotations
        assert ("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "246", "373", "DISULFID", "match_column") in processed_annotations


def test_process_annotation_paired_success(
    logger, multi_logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair, annotation_dict_205_Q9NU22,
    annotation_dict_246_Q9NU22, mock_make_anno_total_disulfid_return_205_P15005,
    mock_map_filter_disulfid_return,
    target_sequence_Q9NU22, target_sequence_continuous_Q9NU22, annot_sequence_Q9NU22
):
    """Test processing of paired DISULFID annotations where paired pair succeeds"""

    processed_annotations = set()

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch("transfer_annotations.map_and_filter_annot_pos", return_value=mock_map_filter_disulfid_return) as mock_map_filter, \
         patch("transfer_annotations.add_to_transfer_dict") as mock_add_transfer:
        process_annotation(
            res_hit=True,
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID", "match_column")
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair
        )

        # Verify map_and_filter_annot_pos was called correctly
        mock_map_filter.assert_called_once_with(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            processed_annotations=processed_annotations,
            paired_annot_pos_str="246",
            caller_target_pos_str="333",
            paired_annotation_dict=annotation_dict_246_Q9NU22
        )

        # Verify remove_failed_annotations was called correctly
        mock_add_transfer.assert_called_once_with(
            hit=True,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_continuous_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id="DISULFID | Intrachain (with C-246); in linked form",
            anno_total=mock_make_anno_total_disulfid_return_205_P15005["anno_total"],
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            paired_position_res_hit=True,
            paired_anno_id="DISULFID | Intrachain (with C-205); in linked form",
            paired_anno_total=mock_map_filter_disulfid_return[1]["anno_total"]
        )

        assert ("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID", "match_column") in processed_annotations
        assert ("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "246", "373", "DISULFID", "match_column") in processed_annotations


def test_process_annotation_paired_no_paired_dict(
    logger, multi_logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair, annotation_dict_205_Q9NU22,
    mock_make_anno_total_disulfid_return_205_P15005,
    target_sequence_Q9NU22, target_sequence_continuous_Q9NU22, annot_sequence_Q9NU22
):
    """Test processing of paired DISULFID annotations where no paired annotation dict is found"""

    # Create a modified entry_annotations without the paired annotation
    modified_entry_annotations = {
        "205": entry_annotations_disulfid_pair["205"],
        "246": []  # Empty list to trigger the else branch
    }

    processed_annotations = set()

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch("transfer_annotations.add_to_transfer_dict") as mock_add_transfer:
        process_annotation(
            res_hit=True,
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            entry_annotations=modified_entry_annotations,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID", "match_column")
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            logger=logger,
            entry_annotations=modified_entry_annotations
        )

        mock_add_transfer.assert_called_once_with(
            hit=True,
            multi_logger=multi_logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_continuous_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id="DISULFID | Intrachain (with C-246); in linked form",
            anno_total=mock_make_anno_total_disulfid_return_205_P15005["anno_total"],
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            paired_position_res_hit=False,
            paired_anno_id=None,
            paired_anno_total=None
        )

        assert ("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID", "match_column") in processed_annotations

def test_process_annotation_paired_gapped_paired(
    logger, multi_logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair, annotation_dict_205_Q9NU22,
    mock_make_anno_total_disulfid_return_205_P15005,
    target_sequence_Q9NU22, target_sequence_continuous_Q9NU22, annot_sequence_Q9NU22
):
    """Test processing paired annotation where paired position is a gap"""
    mock_paired_result_dict = {
        "gapped_paired": True,
        "insert_column_paired": False,
    }
    mock_map_filter_gapped_return = (False, mock_paired_result_dict)
    processed_annotations = set()

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch("transfer_annotations.map_and_filter_annot_pos", return_value=mock_map_filter_gapped_return) as mock_map_filter, \
         patch("transfer_annotations.add_to_transfer_dict") as mock_add_transfer:

        process_annotation(
            res_hit=True,
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID", "match_column")
        )

        expected_paired_key = (
            entry_mnemo_name_mock_Q9NU22,
            target_name_mock_Q9NU22,
            "246",  # paired_annot_pos_str
            "0",    # paired_target_position_str for gap
            "DISULFID",
            "gapped_target_position"
        )

        assert expected_paired_key in processed_annotations
        assert mock_add_transfer.call_args[1]["paired_anno_total"] is None
        assert mock_add_transfer.call_args[1]["paired_anno_id"] is None

def test_process_annotation_paired_insert_column(
    logger, multi_logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair, annotation_dict_205_Q9NU22,
    mock_make_anno_total_disulfid_return_205_P15005,
    target_sequence_Q9NU22, target_sequence_continuous_Q9NU22, annot_sequence_Q9NU22
):
    """Test processing paired annotation where paired position is in insert column"""
    mock_paired_result_dict = {
        "gapped_paired": False,
        "insert_column_paired": True,
        "anno_total": None,
        "anno_id": None
    }
    mock_map_filter_insert_return = (False, mock_paired_result_dict)
    processed_annotations = set()

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch("transfer_annotations.map_and_filter_annot_pos", return_value=mock_map_filter_insert_return) as mock_map_filter, \
         patch("transfer_annotations.add_to_transfer_dict") as mock_add_transfer:

        process_annotation(
            res_hit=True,
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID", "match_column")
        )

        expected_paired_key = (
            entry_mnemo_name_mock_Q9NU22,
            target_name_mock_Q9NU22,
            "246",  # paired_annot_pos_str
            "0",    # paired_target_position_str for insert
            "DISULFID",
            "insert_column"
        )

        assert expected_paired_key in processed_annotations
        assert mock_add_transfer.call_args[1]["paired_anno_total"] is None
        assert mock_add_transfer.call_args[1]["paired_anno_id"] is None

def test_process_annotation_paired_match_column(
    logger, multi_logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair, annotation_dict_205_Q9NU22,
    mock_make_anno_total_disulfid_return_205_P15005,
    target_sequence_Q9NU22, target_sequence_continuous_Q9NU22, annot_sequence_Q9NU22
):
    """Test processing paired annotation with successful match"""
    mock_paired_result_dict = {
        "gapped_paired": False,
        "insert_column_paired": False,
        "anno_total": {
            "target_position": "373",
            "annot_position": "246"
        },
        "anno_id": "DISULFID | Test paired"
    }
    mock_map_filter_match_return = (True, mock_paired_result_dict)
    processed_annotations = set()

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch("transfer_annotations.map_and_filter_annot_pos", return_value=mock_map_filter_match_return) as mock_map_filter, \
         patch("transfer_annotations.add_to_transfer_dict") as mock_add_transfer:

        process_annotation(
            res_hit=True,
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID", "match_column")
        )

        expected_paired_key = (
            entry_mnemo_name_mock_Q9NU22,
            target_name_mock_Q9NU22,
            "246",  # paired_annot_pos_str
            "373",  # actual target position for match
            "DISULFID",
            "match_column"
        )

        assert expected_paired_key in processed_annotations
        assert mock_add_transfer.call_args[1]["paired_anno_total"] == mock_paired_result_dict["anno_total"]
        assert mock_add_transfer.call_args[1]["paired_anno_id"] == mock_paired_result_dict["anno_id"]
        expected_log_message = 'TRANSFER_ANNOTS --- PROCESS_ANNOT --- SUCCESS - ROUTE: Paired - Added to transfer_dict for target sp|Q9NU22|MDN1_HUMAN and annotated MCRB_ECOLI at target_caller 333 with paired 373 --- Caller-paired positions: 205-246 --- ORIGIN: Type DISULFID and Evidence ECO:0000269|PubMed:12345678'
        assert any(call_args[0][0] == expected_log_message for call_args in logger.debug.call_args_list)

###T validate_paired_annotations

def test_validate_paired_annotations_valid_case(
    logger,
    good_eco_codes_all,
    target_sequence_Q9NU22,
    annot_sequence_Q9NU22,
    entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
):
    """
    Test the case where paired annotations are valid and successfully validated.
    Here we assume this is running for the late pair from a successful run of the early pair:
    (validate_annotations() -> process_annotation() -> validate_paired_annotations()) == processed_annotations = set()
    """
    # Mock the output of map_and_filter_annot_pos for paired validation
    mock_paired_result_dict = {
        "annotation": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "evidence": "ECO:0000269|PubMed:12345678",
            "entry": "P15005",
            "aminoacid": "C",
            "paired_position": "205",
        },
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-205); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "373",
            "annot_position": "246",
            "index_position": "72",
            "annot_amino_acid": "C",
            "target_amino_acid": "C"
        },
        "paired_annot_pos_str": "205",
        "gapped_paired": False,
        "insert_column_paired": False
    }

    mock_result_tuple = (True, mock_paired_result_dict)

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_paired_result_dict) as mock_make_total:
        result_tuple = validate_paired_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            paired_annotation_dict=annotation_dict_246_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            counter_target_pos=counter_target_pos_mock, # Always None
            counter_annot_pos=counter_annot_pos_mock, # Always None
            paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22, # "246"
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22, # "333"
            processed_annotations=set()
        )

        # Verify the make_anno_total_dict call
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name="MCRB_ECOLI",
            annotation_dict=annotation_dict_246_Q9NU22,
            counter_target_pos_str="373", # 246 = 373 in target sequence numbering/col
            counter_annot_pos_str=paired_annot_pos_str_mock_Q9NU22, # "246"
            index=72,
            target_amino="C",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair,
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
        )

        # Assert that the result matches the mock data
        assert result_tuple == mock_result_tuple


def test_validate_paired_annotations_invalid_pair(
    logger,
    good_eco_codes_all,
    target_sequence_Q9NU22,
    annot_sequence_Q9NU22,
    entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
):
    """
    Test case where the paired annotation validation fails.
    Here we assume this is running for the late pair from a successful run of the early pair:
    (validate_annotations() -> process_annotation() -> validate_paired_annotations()) == processed_annotations = set()
    """

    # Modify target sequence to ensure failure at paired position
    target_sequence_mod_to_fail_73 = list(target_sequence_Q9NU22)
    target_sequence_mod_to_fail_73[72] = "R"
    target_sequence_mod_to_fail_73 = "".join(target_sequence_mod_to_fail_73)

    mock_paired_result_dict = {
        "annotation": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "evidence": "ECO:0000269|PubMed:12345678",
            "entry": "P15005",
            "aminoacid": "C",
            "paired_position": "205",
        },
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-205); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "373",
            "annot_position": "246",
            "index_position": "72",
            "annot_amino_acid": "C",
            "target_amino_acid": "R"
        },
        "paired_annot_pos_str": "205",
        "gapped_paired": False,
        "insert_column_paired": False
    }

    mock_result_tuple = (False, mock_paired_result_dict)

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_paired_result_dict) as mock_make_total:
        result_tuple = validate_paired_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_mod_to_fail_73,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            paired_annotation_dict=annotation_dict_246_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            counter_target_pos=counter_target_pos_mock, # Always None
            counter_annot_pos=counter_annot_pos_mock, # Always None
            paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22, # "246"
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22, # "333"
            processed_annotations=set()
        )

        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict_246_Q9NU22,
            counter_target_pos_str="373", # 246 = 373 in target sequence numbering/col
            counter_annot_pos_str=paired_annot_pos_str_mock_Q9NU22, # "246"
            index=72,
            target_amino="R",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair,
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
        )

        assert result_tuple == mock_result_tuple

def test_validate_paired_annotations_missing_paired_position(
    logger,
    good_eco_codes_all,
    target_sequence_Q9NU22,
    annot_sequence_Q9NU22,
    entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
):
    """
    Test case where the paired position is missing in the annotations.
    Here we assume this is running for the late pair from a successful run of the early pair:
    (validate_annotations() -> process_annotation() -> validate_paired_annotations()) == processed_annotations = set()
    """
    # Set paired position to a non-existent value
    annotation_dict_246_Q9NU22["paired_position"] = "999"

    result = validate_paired_annotations(
        logger=logger,
        good_eco_codes=good_eco_codes_all,
        target_sequence=target_sequence_Q9NU22,
        target_name=target_name_mock_Q9NU22,
        target_hit_start=target_hit_start_mock_Q9NU22,
        target_hit_end=target_hit_end_mock_Q9NU22,
        offset_start=offset_start_mock_Q9NU22,
        offset_end=offset_end_mock_Q9NU22,
        annot_sequence=annot_sequence_Q9NU22,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        paired_annotation_dict=annotation_dict_246_Q9NU22,
        entry_annotations=entry_annotations_disulfid_pair,
        counter_target_pos=counter_target_pos_mock,
        counter_annot_pos=counter_annot_pos_mock,
        paired_annot_pos_str="999",
        caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
        processed_annotations=set()
    )

    assert result == (False, {"gapped_paired": False, "insert_column_paired": False})

def test_validate_paired_annotations_earlier_gap(
    logger,
    good_eco_codes_all,
    target_sequence_Q9NU22,
    annot_sequence_Q9NU22,
    entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
):
    """Test case where first position of pair was a gap and is in processed_annotations"""

    processed_annotations = {
        (
            entry_mnemo_name_mock_Q9NU22,
            target_name_mock_Q9NU22,
            "246",  # paired_annot_pos_str
            "0",    # failure_target_pos
            "DISULFID",
            "gapped_target_position"
        )
    }

    result = validate_paired_annotations(
        logger=logger,
        good_eco_codes=good_eco_codes_all,
        target_sequence=target_sequence_Q9NU22,
        target_name=target_name_mock_Q9NU22,
        target_hit_start=target_hit_start_mock_Q9NU22,
        target_hit_end=target_hit_end_mock_Q9NU22,
        offset_start=offset_start_mock_Q9NU22,
        offset_end=offset_end_mock_Q9NU22,
        annot_sequence=annot_sequence_Q9NU22,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        paired_annotation_dict=annotation_dict_246_Q9NU22,
        entry_annotations=entry_annotations_disulfid_pair,
        counter_target_pos=counter_target_pos_mock,
        counter_annot_pos=counter_annot_pos_mock,
        paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22,
        caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
        processed_annotations=processed_annotations
    )

    assert result == (False, {"gapped_paired": True, "insert_column_paired": False})

def test_validate_paired_annotations_earlier_insert(
    logger,
    good_eco_codes_all,
    target_sequence_Q9NU22,
    annot_sequence_Q9NU22,
    entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
):
    """Test case where first position of pair was in insert column and is in processed_annotations"""

    processed_annotations = {
        (
            entry_mnemo_name_mock_Q9NU22,
            target_name_mock_Q9NU22,
            "246",  # paired_annot_pos_str
            "0",    # failure_target_pos
            "DISULFID",
            "insert_column"
        )
    }

    result = validate_paired_annotations(
        logger=logger,
        good_eco_codes=good_eco_codes_all,
        target_sequence=target_sequence_Q9NU22,
        target_name=target_name_mock_Q9NU22,
        target_hit_start=target_hit_start_mock_Q9NU22,
        target_hit_end=target_hit_end_mock_Q9NU22,
        offset_start=offset_start_mock_Q9NU22,
        offset_end=offset_end_mock_Q9NU22,
        annot_sequence=annot_sequence_Q9NU22,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        paired_annotation_dict=annotation_dict_246_Q9NU22,
        entry_annotations=entry_annotations_disulfid_pair,
        counter_target_pos=counter_target_pos_mock,
        counter_annot_pos=counter_annot_pos_mock,
        paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22,
        caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
        processed_annotations=processed_annotations
    )

    assert result == (False, {"gapped_paired": False, "insert_column_paired": True})

def test_validate_paired_annotations_insert_column_position(
    logger,
    good_eco_codes_all,
    annot_sequence_Q9NU22,
    target_sequence_Q9NU22,
    entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
):
    """Test case where paired position is in an insert column"""

    # Modify sequences to create insert column at paired position
    annot_sequence_modified = annot_sequence_Q9NU22.replace("GYCPN", "gycpn")
    target_sequence_modified = target_sequence_Q9NU22.replace("MYCCT", ".....")

    result = validate_paired_annotations(
        logger=logger,
        good_eco_codes=good_eco_codes_all,
        target_sequence=target_sequence_modified,
        target_name=target_name_mock_Q9NU22,
        target_hit_start=target_hit_start_mock_Q9NU22,
        target_hit_end=target_hit_end_mock_Q9NU22,
        offset_start=offset_start_mock_Q9NU22,
        offset_end=offset_end_mock_Q9NU22,
        annot_sequence=annot_sequence_modified,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        paired_annotation_dict=annotation_dict_246_Q9NU22,
        entry_annotations=entry_annotations_disulfid_pair,
        counter_target_pos=counter_target_pos_mock,
        counter_annot_pos=counter_annot_pos_mock,
        paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22,
        caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
        processed_annotations=set()
    )

    assert result == (False, {"gapped_paired": False, "insert_column_paired": True})

def test_validate_paired_annotations_gapped_position(
    logger,
    good_eco_codes_all,
    annot_sequence_Q9NU22,
    target_sequence_Q9NU22,
    entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
):
    """Test case where paired position aligns to a gap in target sequence"""

    # Modify target sequence to create gap at paired position
    target_sequence_modified = target_sequence_Q9NU22.replace("MYCCT", "M---T")

    result = validate_paired_annotations(
        logger=logger,
        good_eco_codes=good_eco_codes_all,
        target_sequence=target_sequence_modified,
        target_name=target_name_mock_Q9NU22,
        target_hit_start=target_hit_start_mock_Q9NU22,
        target_hit_end=target_hit_end_mock_Q9NU22,
        offset_start=offset_start_mock_Q9NU22,
        offset_end=offset_end_mock_Q9NU22,
        annot_sequence=annot_sequence_Q9NU22,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        paired_annotation_dict=annotation_dict_246_Q9NU22,
        entry_annotations=entry_annotations_disulfid_pair,
        counter_target_pos=counter_target_pos_mock,
        counter_annot_pos=counter_annot_pos_mock,
        paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22,
        caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
        processed_annotations=set()
    )

    assert result == (False, {"gapped_paired": True, "insert_column_paired": False})

def test_validate_paired_annotations_gap_in_first(
    logger,
    good_eco_codes_all,
    target_sequence_Q9NU22,
    annot_sequence_Q9NU22,
    entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
):
    """
    Test case where the first member of the pair had a gap
    (earlier_annotation_key present in processed_annotations).
    The paired validation itself succeeds (residues match), but caller_had_gap is True.
    """
    processed_annotations = {
        (
            entry_mnemo_name_mock_Q9NU22,  # "MCRB_ECOLI"
            target_name_mock_Q9NU22,       # "sp|Q9NU22|MDN1_HUMAN"
            "205",                         # paired_annot_pos_str
            None,                          # No residue in target position (gap)
            "DISULFID"                     # anno_type
        )
    }

    mock_paired_result_dict = {
        "annotation": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "evidence": "ECO:0000269|PubMed:12345678",
            "entry": "P15005",
            "aminoacid": "C",
            "paired_position": "205",
        },
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-205); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-205); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "373",
            "annot_position": "246",
            "index_position": "72",
            "annot_amino_acid": "C",
            "target_amino_acid": "C"
        },
        "paired_annot_pos_str": "205",
        "gapped_paired": True,  # What we're testing for
        "insert_column_paired": False
    }
    print(json.dumps(mock_paired_result_dict, indent = 4))
    mock_result_tuple = (True, mock_paired_result_dict)

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_paired_result_dict) as mock_make_total:
        result_tuple = validate_paired_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            paired_annotation_dict=annotation_dict_246_Q9NU22,
            entry_annotations=entry_annotations_disulfid_pair,
            counter_target_pos=counter_target_pos_mock,
            counter_annot_pos=counter_annot_pos_mock,
            paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22,
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
            processed_annotations=processed_annotations
        )

        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name="MCRB_ECOLI",
            annotation_dict=annotation_dict_246_Q9NU22,
            counter_target_pos_str="373",
            counter_annot_pos_str=paired_annot_pos_str_mock_Q9NU22,
            index=72,
            target_amino="C",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair,
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
        )

        assert result_tuple == mock_result_tuple
        # Additional specific assertion for what we're testing
        assert result_tuple[1]["gapped_paired"] is True

###T validate_annotations

def test_validate_annotations_process_annotation_parameters(
    logger, multi_logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, transfer_dict
):
    """Test that process_annotation is called with correct parameters"""

    with patch("transfer_annotations.process_annotation") as mock_process:
        validate_annotations(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            counter_target_pos=counter_target_pos_mock,  # Known position where we expect a match
            counter_annot_pos=counter_annot_pos_mock  # Known position with annotation
        )

        expected_calls = [
            call(
                res_hit=True,
                logger=logger,
                multi_logger=multi_logger,
                good_eco_codes=good_eco_codes_all,
                entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
                target_name=target_name_mock_Q9NU22,
                target_hit_start=target_hit_start_mock_Q9NU22,
                target_hit_end=target_hit_end_mock_Q9NU22,
                annotation_dict=entry_annotations_binding_only["201"][0],
                entry_annotations=entry_annotations_binding_only,
                transfer_dict=transfer_dict,
                target_sequence=target_sequence_Q9NU22,
                offset_start=offset_start_mock_Q9NU22,
                offset_end=offset_end_mock_Q9NU22,
                annot_sequence=annot_sequence_Q9NU22,
                counter_target_pos_str="329",
                counter_annot_pos_str="201",
                index=14,
                target_amino="G",
                processed_annotations=set(),
                annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "201", "329", "BINDING", "match_column")
            ),
            call(
                res_hit=True,
                logger=logger,
                multi_logger=multi_logger,
                good_eco_codes=good_eco_codes_all,
                entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
                target_name=target_name_mock_Q9NU22,
                target_hit_start=target_hit_start_mock_Q9NU22,
                target_hit_end=target_hit_end_mock_Q9NU22,
                annotation_dict=entry_annotations_binding_only["202"][0],
                entry_annotations=entry_annotations_binding_only,
                transfer_dict=transfer_dict,
                target_sequence=target_sequence_Q9NU22,
                offset_start=offset_start_mock_Q9NU22,
                offset_end=offset_end_mock_Q9NU22,
                annot_sequence=annot_sequence_Q9NU22,
                counter_target_pos_str="330",
                counter_annot_pos_str="202",
                index=15,
                target_amino="P",
                processed_annotations=set(),
                annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "202", "330", "BINDING", "match_column")
            )
        ]

        mock_process.assert_has_calls(expected_calls, any_order=False)

def test_validate_annotations_multiple_annotations(
    logger, multi_logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    transfer_dict
):
    """Test case where a single position has multiple annotations."""
    entry_annotations = {
        "204": [
            {"type": "BINDING", "description": "Interacts with ATP", "evidence": "ECO:0000255"},
            {"type": "ACT_SITE", "description": "Catalytic activity", "evidence": "ECO:0000269"}
        ]
    }

    with patch("transfer_annotations.process_annotation") as mock_process:
        validate_annotations(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            counter_target_pos=counter_target_pos_mock,
            counter_annot_pos=counter_annot_pos_mock
        )

        # Ensure all annotations at the position are processed
        assert mock_process.call_count == 2

def test_validate_annotations_skip_processed(
    logger, multi_logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, transfer_dict
):
    """Test case 1.2 - Skip when annotation was already processed"""

    processed_annotations = {
        (
            "MCRB_ECOLI",  # entry_mnemo_name
            "sp|Q9NU22|MDN1_HUMAN",  # target_name
            "201",  # counter_annot_pos_str
            "329",  # counter_target_pos_str
            "BINDING",  # anno_type
            "match_column"  # status
        ),
        (
            "MCRB_ECOLI",
            "sp|Q9NU22|MDN1_HUMAN",
            "202",
            "330",
            "BINDING",
            "match_column"
        )
    }

    with patch("transfer_annotations.process_annotation") as mock_process:
        validate_annotations(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            processed_annotations=processed_annotations,
            counter_target_pos=counter_target_pos_mock,
            counter_annot_pos=counter_annot_pos_mock
        )

        # Verify process_annotation was not called
        mock_process.assert_not_called()
        assert transfer_dict == {}

def test_validate_annotations_no_positions_in_range(
    logger, multi_logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, transfer_dict
):
    """Test case 2 - No annotations in offset range (201 and 202 don't fit in 500-600)"""
    with patch("transfer_annotations.process_annotation") as mock_process:
        validate_annotations(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_Q9NU22,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=500,
            offset_end=600,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            counter_target_pos=counter_target_pos_mock,
            counter_annot_pos=counter_annot_pos_mock
        )

        mock_process.assert_not_called()
        assert transfer_dict == {}

def test_validate_annotations_sequence_end_reached(
    logger, multi_logger, good_eco_codes_all,
    entry_annotations_binding_only, transfer_dict
):
    """Test case where sequence reaches end condition (line 633)"""
    # Use sequences guaranteed to hit end condition
    test_short_sequence = ".....ABC"  # Short sequence to ensure we hit end

    validate_annotations(
        logger=logger,
        multi_logger=multi_logger,
        good_eco_codes=good_eco_codes_all,
        target_sequence=test_short_sequence,
        target_name=target_name_mock_Q9NU22,
        target_hit_start=1,
        target_hit_end=3,  # Will trigger end condition when counter_target_pos == target_hit_end
        offset_start=1,
        offset_end=3,
        annot_sequence=test_short_sequence,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        entry_annotations=entry_annotations_binding_only,
        transfer_dict=transfer_dict,
        processed_annotations=set(),
        counter_target_pos=None,
        counter_annot_pos=None
    )

def test_validate_annotations_logging(
    logger, multi_logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, transfer_dict
):
    """Test that error logs are written when process_annotation fails."""
    mock_error = Exception("Test error")

    with patch("transfer_annotations.process_annotation", side_effect=mock_error):
        try:
            validate_annotations(
                logger=logger,
                multi_logger=multi_logger,
                good_eco_codes=good_eco_codes_all,
                target_sequence=target_sequence_Q9NU22,
                target_name=target_name_mock_Q9NU22,
                target_hit_start=target_hit_start_mock_Q9NU22,
                target_hit_end=target_hit_end_mock_Q9NU22,
                offset_start=201,
                offset_end=202,
                annot_sequence=annot_sequence_Q9NU22,
                entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
                entry_annotations=entry_annotations_binding_only,
                transfer_dict=transfer_dict,
                processed_annotations=set(),
                counter_target_pos=counter_target_pos_mock,
                counter_annot_pos=counter_annot_pos_mock
            )
        except Exception:
            pass

        # Verify the call was made
        multi_logger.assert_called_once()

        # Get the actual call arguments
        call_args = multi_logger.call_args

        # Verify the arguments individually
        assert call_args[0][0] == "error"  # First arg should be level
        assert call_args[0][1] == "TRANSFER_ANNOTS --- VAL_ANNOTS --- Downstream error in process_annotation: %s"  # Second arg should be message format
        assert str(call_args[0][2]) == str(mock_error)  # Compare string representations of exceptions

def test_validate_annotations_insert_column_position(
    logger, multi_logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22, transfer_dict, entry_annotations_binding_only
):
    """Test case for handling insert columns in annotated positions in source sequence"""
    # Modified sequence with lowercase letters to trigger insert column condition
    annot_sequence_mod_insert_14_to_16 = annot_sequence_Q9NU22.replace("GPP", "gpp")
    target_sequence_mod_insert_14_to_16 = target_sequence_Q9NU22.replace("GPI", "...")

    # Simulate a paired BINDING pair
    entry_annotations_binding_only_paired = copy.deepcopy(entry_annotations_binding_only)
    entry_annotations_binding_only_paired["201"][0]["paired_position"] = "202"
    entry_annotations_binding_only_paired["202"][0]["paired_position"] = "201"

    expected_insert_key = (
        entry_mnemo_name_mock_Q9NU22,
        target_name_mock_Q9NU22,
        "201",  # position matching the insert
        "0",    # failure_target_pos
        "BINDING",
        "insert_column"
    )

    processed_annotations_container = set()

    with patch("transfer_annotations.process_annotation") as mock_process:
        validate_annotations(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_mod_insert_14_to_16,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_mod_insert_14_to_16,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_binding_only_paired,
            transfer_dict=transfer_dict,
            processed_annotations=processed_annotations_container,
            counter_target_pos=counter_target_pos_mock,
            counter_annot_pos=counter_annot_pos_mock
        )

        # Verify insert column was handled (process_annotation shouldn't be called)
        mock_process.assert_not_called()
        assert expected_insert_key in processed_annotations_container

def test_validate_annotations_gapped_target_position(
    logger, multi_logger, target_sequence_Q9NU22, annot_sequence_Q9NU22, good_eco_codes_all, transfer_dict, entry_annotations_binding_only
):
    """Test case for handling annotations at positions that are gaps in target sequence"""

    # Placed gaps in target_seq at annotated position - no change needed in annot_seq
    target_sequence_mod_gapped_14_to_16 = target_sequence_Q9NU22.replace("GPI", "---")

    entry_annotations_binding_only_paired = copy.deepcopy(entry_annotations_binding_only)
    entry_annotations_binding_only_paired["201"][0]["paired_position"] = "202"
    entry_annotations_binding_only_paired["202"][0]["paired_position"] = "201"

    expected_insert_key = (
        entry_mnemo_name_mock_Q9NU22,
        target_name_mock_Q9NU22,
        "201",  # position matching the insert
        "0",    # failure_target_pos
        "BINDING",
        "gapped_target_position"
    )

    processed_annotations_container = set()

    with patch("transfer_annotations.process_annotation") as mock_process:
        validate_annotations(
            logger=logger,
            multi_logger=multi_logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_mod_gapped_14_to_16,
            target_name=target_name_mock_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            offset_start=offset_start_mock_Q9NU22,
            offset_end=offset_end_mock_Q9NU22,
            annot_sequence=annot_sequence_Q9NU22,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_annotations=entry_annotations_binding_only_paired,
            transfer_dict=transfer_dict,
            processed_annotations=processed_annotations_container,
            counter_target_pos=counter_target_pos_mock,
            counter_annot_pos=counter_annot_pos_mock
        )

        # Verify gap was handled (process_annotation shouldn't be called)
        mock_process.assert_not_called()
        assert expected_insert_key in processed_annotations_container

###T main

def test_main_success(logger, multi_logger, minimal_hmmalign_lines_fixture_Q9NU22, transfer_dict_populated_disulfid_list_Q9NU22, transfer_dict_populated_disulfid_post_gos_list_Q9NU22):
    """Test main function success path"""
    mock_args = Namespace(
        dom_align=hmmalign_result_mock,
        resource_dir=resource_dir_mock,
        domain_accession=domain_accession_mock,
        output_dir=output_dir_mock,
        eco_codes=good_eco_codes_mock,
        log=log_filepath_mock
    )

    with patch("transfer_annotations.parse_arguments", return_value=mock_args) as mock_parse, \
         patch("transfer_annotations.get_pfam_id_from_hmmalign_result", return_value="PF07728") as mock_get_pfam, \
         patch("transfer_annotations.read_files", return_value=(minimal_hmmalign_lines_fixture_Q9NU22, {})) as mock_read, \
         patch("transfer_annotations.find_and_map_annots", return_value=transfer_dict_populated_disulfid_list_Q9NU22) as mock_find, \
         patch("transfer_annotations.cleanup_improve_transfer_dict", return_value=transfer_dict_populated_disulfid_post_gos_list_Q9NU22) as mock_cleanup, \
         patch("transfer_annotations.write_reports") as mock_write, \
         patch("transfer_annotations.get_logger", return_value=(logger, None)) as mock_get_logger, \
         patch("transfer_annotations.get_multi_logger") as mock_get_multi_logger:

        mock_get_multi_logger.return_value = multi_logger

        main()

        mock_parse.assert_called_once()
        mock_get_pfam.assert_called_once_with(hmmalign_result_mock)
        mock_read.assert_called_once()
        mock_find.assert_called_once()
        mock_cleanup.assert_called_once_with(
            logger,
            multi_logger,
            transfer_dict_populated_disulfid_list_Q9NU22,
            "PF07728",
            minimal_hmmalign_lines_fixture_Q9NU22,
            os.path.join(resource_dir_mock, "PF07728", "conservations.json"),
            os.path.join(resource_dir_mock, "PF07728", "annotations.json"),
            output_dir_mock,
            os.path.join(resource_dir_mock, "mappings/interpro_pfam_accession_mapping.tsv")
        )
        mock_write.assert_called_once_with(
            logger,
            multi_logger,
            transfer_dict_populated_disulfid_post_gos_list_Q9NU22,
            output_dir_mock
        )
        logger.info.assert_any_call("TRANSFER_ANNOTS --- MAIN --- Transfer Dict FILLED")


# Integration tests
def test_main_integration_binding_Q9NU22_PF07728(
    minimal_hmmalign_lines_fixture_Q9NU22, transfer_dict_success_binding_Q9NU22, annotations_content_binding_fixture_Q9NU22_PF07728,
    conservations_content_Q9NU22_PF07728, mapping_content_Q9NU22_and_H0YB80_domains, iprscan_content_Q9NU22
    ):
    """Integration test for main function"""

    with TemporaryDirectory() as tmp_dir:
        # Define paths
        hmmalign_path = os.path.join(tmp_dir, "PF07728_hmmalign.sth")
        output_dir = os.path.join(tmp_dir, "output")
        resource_dir = os.path.join(tmp_dir, "resources")
        pfam_dir = os.path.join(resource_dir, "PF07728")
        mappings_dir = os.path.join(resource_dir, "mappings")
        target_dir = os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN")
        pfam_output_dir = os.path.join(output_dir, "PF07728")

        # Create needed directories
        for directory in [output_dir, resource_dir, pfam_dir, mappings_dir, target_dir, pfam_output_dir]:
            os.makedirs(directory, exist_ok=True)

        # Write the minimal HMM align file
        with open(hmmalign_path, "w", encoding="latin-1") as f:
            f.write(minimal_hmmalign_lines_fixture_Q9NU22)

        # Write annotation and conservation JSON
        with open(os.path.join(pfam_dir, "annotations.json"), "w", encoding="utf-8") as f:
            json.dump(annotations_content_binding_fixture_Q9NU22_PF07728, f)

        with open(os.path.join(pfam_dir, "conservations.json"), "w", encoding="utf-8") as f:
            json.dump(conservations_content_Q9NU22_PF07728, f)

        # Write mappings TSV
        with open(os.path.join(mappings_dir, "interpro_pfam_accession_mapping.tsv"), "w", encoding="utf-8") as f:
            f.write(mapping_content_Q9NU22_and_H0YB80_domains)

        # Write iprscan.tsv
        with open(os.path.join(target_dir, "iprscan.tsv"), "w", encoding="utf-8") as f:
            f.write(iprscan_content_Q9NU22)

        # Run main
        args = Namespace(
            dom_align=hmmalign_path,
            resource_dir=resource_dir,
            domain_accession=domain_accession_mock,
            output_dir=output_dir,
            eco_codes=good_eco_codes_mock,
            log=os.path.join(tmp_dir, "test.log")
        )

        with patch("transfer_annotations.parse_arguments", return_value=args):
            logger, _ = get_logger(args.log)
            main()

        # with open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json")) as f:
        #     transfer_dict = json.load(f)
        #     print(json.dumps(transfer_dict, indent=4))

        assert os.path.exists(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json"))
        assert transfer_dict_success_binding_Q9NU22 == json.load(open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json"), "r", encoding="utf-8"))

def test_main_integration_disulfid_Q9NU22_PF07728(
    minimal_hmmalign_lines_fixture_Q9NU22,
    transfer_dict_populated_disulfid_post_gos_list_Q9NU22,
    annotations_content_disulfid_fixture_Q9NU22_PF07728,
    conservations_content_Q9NU22_PF07728,
    mapping_content_Q9NU22_and_H0YB80_domains,
    iprscan_content_Q9NU22):
    """Integration test for main function"""

    with TemporaryDirectory() as tmp_dir:
        # Define paths
        hmmalign_path = os.path.join(tmp_dir, "PF07728_hmmalign.sth")
        output_dir = os.path.join(tmp_dir, "output")
        resource_dir = os.path.join(tmp_dir, "resources")
        pfam_dir = os.path.join(resource_dir, "PF07728")
        mappings_dir = os.path.join(resource_dir, "mappings")
        target_dir = os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN")
        pfam_output_dir = os.path.join(output_dir, "PF07728")

        # Create needed directories
        for directory in [output_dir, resource_dir, pfam_dir, mappings_dir, target_dir, pfam_output_dir]:
            os.makedirs(directory, exist_ok=True)

        # Write the minimal HMM align file
        with open(hmmalign_path, "w", encoding="latin-1") as f:
            f.write(minimal_hmmalign_lines_fixture_Q9NU22)

        # Write annotation and conservation JSON
        with open(os.path.join(pfam_dir, "annotations.json"), "w", encoding="utf-8") as f:
            json.dump(annotations_content_disulfid_fixture_Q9NU22_PF07728, f)

        with open(os.path.join(pfam_dir, "conservations.json"), "w", encoding="utf-8") as f:
            json.dump(conservations_content_Q9NU22_PF07728, f)

        # Write mappings TSV
        with open(os.path.join(mappings_dir, "interpro_pfam_accession_mapping.tsv"), "w", encoding="utf-8") as f:
            f.write(mapping_content_Q9NU22_and_H0YB80_domains)

        # Write iprscan.tsv
        with open(os.path.join(target_dir, "iprscan.tsv"), "w", encoding="utf-8") as f:
            f.write(iprscan_content_Q9NU22)

        # Run main
        args = Namespace(
            dom_align=hmmalign_path,
            resource_dir=resource_dir,
            domain_accession=domain_accession_mock,
            output_dir=output_dir,
            eco_codes=good_eco_codes_mock,
            log=os.path.join(tmp_dir, "test.log")
        )

        with patch("transfer_annotations.parse_arguments", return_value=args):
            logger, _ = get_logger(args.log)
            main()

        # with open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json")) as f:
        #     transfer_dict = json.load(f)
        #     print(json.dumps(transfer_dict, indent=4))

        assert os.path.exists(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json"))
        assert transfer_dict_populated_disulfid_post_gos_list_Q9NU22 == json.load(open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json"), "r", encoding="utf-8"))

def test_main_integration_all_types_H0YB80(
    minimal_hmmalign_lines_fixture_H0YB80,
    transfer_dict_success_all_types_H0YB80,
    annotations_content_all_types_fixture_H0YB80_PF00244,
    conservations_content_H0YB80_PF00244,
    mapping_content_Q9NU22_and_H0YB80_domains,
    iprscan_content_H0YB80):
    """Integration test for main function"""

    with TemporaryDirectory() as tmp_dir:
        # Define paths
        hmmalign_path = os.path.join(tmp_dir, "PF00244_hmmalign.sth")
        output_dir = os.path.join(tmp_dir, "output")
        resource_dir = os.path.join(tmp_dir, "resources")
        pfam_dir = os.path.join(resource_dir, "PF00244")
        mappings_dir = os.path.join(resource_dir, "mappings")
        target_dir = os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN")
        pfam_output_dir = os.path.join(output_dir, "PF00244")

        # Create needed directories
        for directory in [output_dir, resource_dir, pfam_dir, mappings_dir, target_dir, pfam_output_dir]:
            os.makedirs(directory, exist_ok=True)

        # Write the minimal HMM align file
        with open(hmmalign_path, "w", encoding="latin-1") as f:
            f.write(minimal_hmmalign_lines_fixture_H0YB80)

        # Write annotation and conservation JSON
        with open(os.path.join(pfam_dir, "annotations.json"), "w", encoding="utf-8") as f:
            json.dump(annotations_content_all_types_fixture_H0YB80_PF00244, f)

        with open(os.path.join(pfam_dir, "conservations.json"), "w", encoding="utf-8") as f:
            json.dump(conservations_content_H0YB80_PF00244, f)

        # Write mappings TSV
        with open(os.path.join(mappings_dir, "interpro_pfam_accession_mapping.tsv"), "w", encoding="utf-8") as f:
            f.write(mapping_content_Q9NU22_and_H0YB80_domains)

        # Write iprscan.tsv
        with open(os.path.join(target_dir, "iprscan.tsv"), "w", encoding="utf-8") as f:
            f.write(iprscan_content_H0YB80)

        # Run main
        args = Namespace(
            dom_align=hmmalign_path,
            resource_dir=resource_dir,
            domain_accession=domain_accession_mock,
            output_dir=output_dir,
            eco_codes=good_eco_codes_mock,
            log=os.path.join(tmp_dir, "test.log")
        )

        with patch("transfer_annotations.parse_arguments", return_value=args):
            logger, _ = get_logger(args.log)
            main()

        # with open(os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN", "PF00244_report.json")) as f:
        #     transfer_dict = json.load(f)
        #     print(json.dumps(transfer_dict, indent=4))

        assert os.path.exists(os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN", "PF00244_report.json"))
        assert transfer_dict_success_all_types_H0YB80 == json.load(open(os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN", "PF00244_report.json"), "r", encoding="utf-8"))
