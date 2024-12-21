"""
Unit tests for transfer_annotations.py
"""

import logging
import json
import sys
import os
import copy

from importlib import reload
from argparse import Namespace
from unittest.mock import patch, mock_open, MagicMock, ANY, call
from tempfile import TemporaryDirectory
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from transfer_annotations import (
    parse_arguments,
    configure_logging,
    get_pfam_id_from_hmmalign_result,
    get_annotation_filepath,
    read_files,
    find_and_map_annots,
    map_and_filter_annot_pos,
    validate_annotations,
    process_annotation,
    make_anno_total_dict,
    validate_paired_annotations,
    add_to_transfer_dict,
    write_report,
    main
)
import pytest


hmmalign_result_mock = "/home/user/results/human/PF07728_hmmalign.sth"
hmmalign_result_content_mock = """# STOCKHOLM 1.0

#=GS VWA8_HUMAN/105-261                              AC A3KMH1.2
#=GS Q8ZSL8_PYRAE/20-171                             AC Q8ZSL8.1

VWA8_HUMAN/105-261                                      .........DVFLIGPPGPLRRSIAM.QYLELT..............KREVEYIALSR.DT..TETDLKQRREIR............AGTAFYIDQCAVRAAT..................EGRTLILEGLEKAE.R........N....VLP....V......LNN.LLENR.E...MQLEDGRFLMSAERYD.kLLRDhtkkelds............wkivrvsenFRVIALGLPVP........rYSGNPLDPPLRSRF............................
Q19346_CAEEL/78-234                                     .........DVFLLGVPGKIRLELVL.RYLEAT..............NREFEYLPITR.DT..TEADIKQRREIR............DGTAYYTDLCAVRAAL..................KGRVLVIDGVERAE.R........N....VLP....I......LNN.LLENR.E...MQLDDGRFLMKHDKYD.eLKVKydeatlkk............mgmervsenFHVIALGLPVP........rFPGNSLDPPFRSRF............................
MCRB_ECOLI/196-350                                      .........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYRCN............GVGFRRKDGIFYNFCQqak...........eqpeKKYIFIIDEINRANlS........K....VFG....E......VMM.LMEHD.K...RGENWSVPLTYSENDE.eRFYVp..........................enVYIIGLMNTAD.........RSLAVVDYALRRRF............................
sp|Q9NU22|MDN1_HUMANtarget//325-451                     .........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYRCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................FQFFAT-----.........--------------rrllscgg....................
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

output_dir_mock = "/home/user/results/human/PF07728/"

good_eco_codes_mock = ["ECO:0000269", "ECO:0000255", "ECO:0000313", "ECO:0007669"]

log_filepath_mock = "/home/user/logs/transfer_annotations.log"


### Fixtures
@pytest.fixture
def annotations_content_binding_fixture():
    """Base annotations with only BINDING annotations"""
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
            ]
        }
    }

@pytest.fixture
# NOTE: Intrachain doesn't really appear in actual descriptions, added for clarity in my tests.
# Instead, you might see Interchain and, in these cases, there'd be no paired_position key.
def annotations_content_disulfid_fixture():
    """Separate fixture for DISULFID pair annotations"""
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
            ]
        }
    }

@pytest.fixture
def annotations_content_all_types_PF00244_fixture():
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
def logger():
    """Create a properly configured mock logger with all necessary methods"""
    mock_logger = MagicMock()
    mock_logger.error = MagicMock()
    mock_logger.info = MagicMock()
    mock_logger.warning = MagicMock()
    mock_logger.debug = MagicMock()
    return mock_logger

@pytest.fixture
def good_eco_codes_all():
    return ["ECO:0000255", "ECO:0000269", "ECO:0000313", "ECO:0007669"]

@pytest.fixture
def good_eco_codes_no_eco255():
    return ["ECO:0000269", "ECO:0000313", "ECO:0007669"]

@pytest.fixture
def target_sequence():
    return ".........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................FQFFAT-----.........--------------rrllscgg...................."

@pytest.fixture
def annot_sequence():
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
def entry_annotations_disulfid_pair_post_map_and_filter():
    """Just DISULFID pair annotations after map_and_filter_annot_pos in process_annotation"""
    return {
        "205": [
            {
                "type": "DISULFID",
                "description": "Intrachain (with C-246); in linked form",
                "evidence": "ECO:0000269|PubMed:12345678",
                "entry": "P15005",
                "aminoacid": "C", # NIILQGPPGC
                "paired_position": "246",
                "target_position": "333"
            }],
        "246": [
            {
                "type": "DISULFID",
                "description": "Intrachain (with C-205); in linked form",
                "evidence": "ECO:0000269|PubMed:12345678",
                "entry": "P15005",
                "aminoacid": "C", # SYEDFIQGYRCN
                "paired_position": "205",
                "target_position": "373"
            }
        ]
    }

@pytest.fixture
def entry_annotations_full(base_annotations, entry_annotations_disulfid_pair):
    """Complete annotations set"""
    return {**base_annotations, **entry_annotations_disulfid_pair}

@pytest.fixture
def transfer_dict():
    return {}

@pytest.fixture
def transfer_dict_populated():
    return {
        "sp|Q9NU22|MDN1_HUMAN": {
            "match": {
                "333": {
                    "DISULFID | Intrachain (with C-246); in linked form": {
                        "essentials": {
                            "type": "DISULFID",
                            "description": "Intrachain (with C-246); in linked form",
                            "count": 1
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
                            "count": 1
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
            "miss": {}
        }
    }

@pytest.fixture
def transfer_dict_success_disulfid_MDN1_HUMAN():
    return {
    "match": {
        "333": {
            "DISULFID | Intrachain (with C-246); in linked form": {
                "essentials": {
                    "type": "DISULFID",
                    "description": "Intrachain (with C-246); in linked form",
                    "count": 1
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
                    "count": 1
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
    "miss": {
        "681": {
            "DISULFID | Intrachain (with C-246); in linked form": {
                "essentials": {
                    "type": "DISULFID",
                    "description": "Intrachain (with C-246); in linked form",
                    "count": 1
                },
                "evidence": {
                    "ECO:0000269|PubMed:12345678": {
                        "rep_primary_accession": "P15005",
                        "rep_mnemo_name": "MCRB_ECOLI",
                        "count": 1
                    }
                },
                "paired_position": {
                    "717": {
                        "rep_primary_accession": "P15005",
                        "rep_mnemo_name": "MCRB_ECOLI",
                        "count": 1
                    }
                },
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
        "717": {
            "DISULFID | Intrachain (with C-205); in linked form": {
                "essentials": {
                    "type": "DISULFID",
                    "description": "Intrachain (with C-205); in linked form",
                    "count": 1
                },
                "evidence": {
                    "ECO:0000269|PubMed:12345678": {
                        "rep_primary_accession": "P15005",
                        "rep_mnemo_name": "MCRB_ECOLI",
                        "count": 1
                    }
                },
                "paired_position": {
                    "681": {
                        "rep_primary_accession": "P15005",
                        "rep_mnemo_name": "MCRB_ECOLI",
                        "count": 1
                    }
                },
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
    }
}

@pytest.fixture
def transfer_dict_success_binding_MDN1_HUMAN():
    return {
    "match": {
        "329": {
            "BINDING | Interacts with GTP": {
                "essentials": {
                    "type": "BINDING",
                    "description": "Interacts with GTP",
                    "count": 1
                },
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
                }
            }
        },
        "330": {
            "BINDING | Interacts with GTP": {
                "essentials": {
                    "type": "BINDING",
                    "description": "Interacts with GTP",
                    "count": 1
                },
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
                }
            }
        },
        "677": {
            "BINDING | Interacts with GTP": {
                "essentials": {
                    "type": "BINDING",
                    "description": "Interacts with GTP",
                    "count": 1
                },
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
                }
            }
        }
    },
    "miss": {
        "678": {
            "BINDING | Interacts with GTP": {
                "essentials": {
                    "type": "BINDING",
                    "description": "Interacts with GTP",
                    "count": 1
                },
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
                }
            }
        }
    }
}

@pytest.fixture
def transfer_dict_success_all_types_H0YB80_HUMAN():
    return {
    }

## Added in Process_annotation
@pytest.fixture
def get_annotation_dict():
    """Helper fixture that provides individual annotation dictionaries"""
    def get_annotation(fixture, pos, index=0):
        """Extract specific annotation dict from fixture
        Args:
            fixture: The source fixture (for now only
            annotations_content_binding_fixture or annotations_content_disulfid_fixture)
            pos: Position key (e.g. "201")
            index: Index in annotation list (default 0)
        """
        return fixture["MCRB_ECOLI"][pos][index]
    return get_annotation

@pytest.fixture
def annotation_dict_205(annotations_content_disulfid_fixture, get_annotation_dict):
    return get_annotation_dict(annotations_content_disulfid_fixture, "205")

@pytest.fixture
def annotation_dict_246(annotations_content_disulfid_fixture, get_annotation_dict):
    return get_annotation_dict(annotations_content_disulfid_fixture, "246")

@pytest.fixture
def annotation_dict_201(annotations_content_binding_fixture, get_annotation_dict):
    return get_annotation_dict(annotations_content_binding_fixture, "201")

@pytest.fixture
def annotation_dict_202(annotations_content_binding_fixture, get_annotation_dict):
    return get_annotation_dict(annotations_content_binding_fixture, "202")

@pytest.fixture
def annotation_dict_post_make_anno_total(annotation_dict_205):
    annotation_dict = copy.deepcopy(annotation_dict_205)
    return annotation_dict

@pytest.fixture
def paired_annotation_dict_post_map_filter(annotation_dict_246):
    annotation_dict = copy.deepcopy(annotation_dict_246)
    return annotation_dict

@pytest.fixture
def mock_make_anno_total_disulfid_return_205_P15005(annotation_dict_post_make_anno_total):
    return {
        'annotation': annotation_dict_post_make_anno_total,
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-246); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-246); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678',
            'target_position': '333',
            'annot_position': '205'
        },
        'paired_annot_pos_str': '246'
    }

@pytest.fixture
def mock_make_anno_total_disulfid_return_246_P15005(annotation_dict_post_make_anno_total):
   return {
        'annotation': annotation_dict_post_make_anno_total,
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-205); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-205); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678',
            'target_position': '373',
            'annot_position': '246'
        },
        'paired_annot_position': '205'
    }

@pytest.fixture
def mock_make_anno_total_disulfid_return_205_P00750():
    return {
        'annotation': {
                    "type": "DISULFID",
                    "description": "Intrachain (with C-246); in linked form",
                    "evidence": "ECO:0000250", # Changed from 269 to 250, testing different ECOs
                    "entry": "P00750", # Random Mock
                    "aminoacid": "C",
                    "paired_position": "246"
                },
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-246); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-246); in linked form',
            'count': 1,
            "evidence": "ECO:0000250", # Changed from 269 to 250, testing different ECOs
            'target_position': '333',
            'annot_position': '205',
            'paired_target_position': '373'
        },
        'paired_annot_position': '246'
    }

@pytest.fixture
def mock_make_anno_total_disulfid_return_246_P00750():
    return {
        'annotation': {
                    "type": "DISULFID",
                    "description": "Intrachain (with C-205); in linked form",
                    "evidence": "ECO:0000250", # Changed from 269 to 250, testing different ECOs
                    "entry": "P00750", # Random Mock
                    "aminoacid": "C",
                    "paired_position": "205"
                },
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-205); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-205); in linked form',
            'count': 1,
            "evidence": "ECO:0000250", # Changed from 269 to 250, testing different ECOs
            'target_position': '373',
            'annot_position': '246',
            'paired_target_position': '333'
        },
        'paired_annot_position': '205'
    }


@pytest.fixture
def mock_map_filter_disulfid_return(paired_annotation_dict_post_map_filter):
    return (True, {
        'annotation': paired_annotation_dict_post_map_filter,
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-205); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-205); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678',
            'target_position': '373',
            'annot_position': '246'
        },
        'paired_annot_position': '205'
    })

# add_transfer_dict
@pytest.fixture
def mock_transfer_dict():
    """Fixture for initializing a transfer_dict"""
    return {}

@pytest.fixture
def mock_anno_total_binding():
    return {
            'type': 'BINDING',
            'description': 'Interacts with GTP',
            'count': 1,
            'evidence': 'ECO:0000255',
            'ligand_id': 'ChEBI:CHEBI:37565',
            'target_position': '329',
            'annot_position': '201'
    }

@pytest.fixture
def mock_anno_total_disulfid():
    """Fixture for mock annotation data."""
    return {
        "type": "DISULFID",
        "description": "Intrachain (with C-246); in linked form",
        "count": 1,
        "evidence": "ECO:0000269|PubMed:12345678",
        "target_position": "333",
        "annot_position": "205",
        "paired_target_position": "373"
    }

@pytest.fixture
def mock_paired_anno_total_disulfid():
    """Fixture for mock paired annotation data."""
    return {
        "type": "DISULFID",
        "description": "Intrachain (with C-205); in linked form",
        "count": 1,
        "evidence": "ECO:0000269|PubMed:12345678",
        "target_position": "373",
        "annot_position": "246",
        "paired_target_position": "333"
    }

###### Configuration Constants
### Added in Map_and_filter_annot_pos
target_name_mock = "sp|Q9NU22|MDN1_HUMAN"
target_hit_start_mock = 325
target_hit_end_mock = 451
offset_start_mock = 196
offset_end_mock = 350
entry_mnemo_name_mock = "MCRB_ECOLI"

## Extra variables for paired position testing
paired_annot_pos_str_mock = "246"
caller_target_pos_str_mock = "333"
caller_annot_pos_str_mock = "205"
# annotation_dict = use the get_annotation_dict fixture
## Result for running with variables for paired positions
paired_position_res_hit_plus_paired_result_dict_tuple_mock = (True, {'annotation': {'type': 'DISULFID', 'description': 'Intrachain (with C-205); in linked form', 'evidence': 'ECO:0000269|PubMed:12345678', 'entry': 'P15005', 'aminoacid': 'C', 'paired_position': '333', 'target_position': '373'}, 'anno_type': 'DISULFID', 'anno_id': 'DISULFID | Intrachain (with C-205); in linked form', 'anno_total': {'type': 'DISULFID', 'description': 'Intrachain (with C-205); in linked form', 'count': 1, 'evidence': 'ECO:0000269|PubMed:12345678', 'paired_position': '333'}, 'paired_position': '333'})


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
entry_primary_accession_mock = "P15005"

# Repeated Anno ID test - Using _return_348 and _return_725 fixtures
entry_primary_accession_mock_repeated = "P00720"
entry_mnemo_name_mock_repeated = "TPA_HUMAN"

# anno_id = Define in the test, new
# anno_total = Define in the test, new

## If called from a paireable type in process_annotation, also note that:
# paired_annot_pos = Define in the test
# caller_target_pos = Define in the test
# annotation_dict = use the get_annotation_dict fixture

### Added/Modified in _add_single_annotation
# paired_additional_keys/additional_keys = Define in the test
# paired_position_res_hit/hit = Define in the test ???
# paired_anno_id = Define in the test ???
# paired_anno_total = Define in the test ???

### Added/Modified in main
minimal_hmmalign_content_mock = """# STOCKHOLM 1.0
sp|Q9NU22|MDN1_HUMANtarget//325-451                     .........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYRCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................-FQFFAT-----.........--------------rrllscgg....................
sp|Q9NU22|MDN1_HUMANtarget//672-755                     .........PVLLVGETGTGKTSTIQ.YLAHIT..............GHRLRVVNMNQ.QS..DTADLLGGYKP-............-VDHKLIWLPLREAFE..................--------------.-........-....---....-......---.-----.-...----------------..----.............................-----------.........--------------elfaqtfskkqnftflghiqtc......
MCRB_ECOLI/196-350                                      .........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYRCN............GVGFRRKDGIFYNFCQqak...........eqpeKKYIFIIDEINRANlS........K....VFG....E......VMM.LMEHD.K...RGENWSVPLTYSENDE.eRFYVp..........................enVYIIGLMNTAD.........RSLAVVDYALRRRF............................
//
"""



def test_parse_arguments_required():
    test_args = [
        "-iA", hmmalign_result_mock,
        "-r", resource_dir_mock,
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
        output_dir="/home/user/results/human/PF07728/",
        eco_codes=[],
        log="logs/transfer_annotations.log"
    )
    assert vars(args) == vars(expected)

def test_parse_arguments_optional():
    test_args = [
        "-iA", hmmalign_result_mock,
        "-r", resource_dir_mock,
        "-o", output_dir_mock,
        "-e", *good_eco_codes_mock, # Unpacking list of ECO codes to simupaired space-separated input
        "-l", log_filepath_mock
    ]

    with pytest.MonkeyPatch().context() as mp:
        mp.setattr("sys.argv", ["script_name.py"] + test_args)
        args = parse_arguments()

    assert args.log == log_filepath_mock
    assert args.eco_codes == good_eco_codes_mock

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

def test_get_pfam_id_from_hmmalign_result():
    assert get_pfam_id_from_hmmalign_result(hmmalign_result_mock) == "PF07728"
    hmmalign_result_mock_2 = "/home/user/results/human/PF00001_hmmalign.sth"
    assert get_pfam_id_from_hmmalign_result(hmmalign_result_mock_2) == "PF00001"

def test_get_annotation_filepath():
    assert get_annotation_filepath(resource_dir_mock, 'PF07728') == "/home/user/resources/PF07728/annotations.json"

def test_read_files(annotations_content_binding_fixture):
    annotations_content = json.dumps(annotations_content_binding_fixture)

    mock_files = {
        "dummy_hmmalign": hmmalign_result_content_mock,
        "dummy_annotations": annotations_content
    }

    def mock_open_files(file, *args, **kwargs):
        if file in mock_files:
            return mock_open(read_data=mock_files[file])()
        raise FileNotFoundError(f"File not found: {file}")

    with patch("builtins.open", new=mock_open_files):
        hmmalign_lines, annotations = read_files("dummy_hmmalign", "dummy_annotations")
        assert hmmalign_lines == hmmalign_result_content_mock.splitlines(keepends=True)
        assert annotations == annotations_content_binding_fixture

###T find_and_map_annots
def test_find_and_map_annots_basic_case(annotations_content_binding_fixture, logger, target_sequence, annot_sequence):
    """Test basic case with single binding annotation"""
    mock_lines = [
        "# STOCKHOLM 1.0\n",
        "MCRB_ECOLI/196-350      .........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYCPN............GVGFRRKDGIFYNFCQqak...........eqpeKKYIFIIDEINRANlS........K....VFG....E......VMM.LMEHD.K...RGENWSVPLTYSENDE.eRFYVp..........................enVYIIGLMNTAD.........RSLAVVDYALRRRF............................\n",
        "sp|Q9NU22|MDN1_HUMANtarget//325-451                     .........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................FQFFAT-----.........--------------rrllscgg....................\n"
    ]

    with patch('transfer_annotations.map_and_filter_annot_pos') as mock_map_and_filter:
        find_and_map_annots(
            logger=logger,
            hmmalign_lines=mock_lines,
            annotations=annotations_content_binding_fixture,
            good_eco_codes=[]
        )

        mock_map_and_filter.assert_called_with(
            logger=logger,
            good_eco_codes=[],
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_annotations=annotations_content_binding_fixture["MCRB_ECOLI"],
            transfer_dict={},
            processed_annotations=set(),
        )

def test_find_and_map_annots_no_target_sequence(annotations_content_binding_fixture, logger):
    """Test case with no matching annotations"""
    mock_lines = [
        "# STOCKHOLM 1.0\n",
        "VWA8_HUMAN/105-261    .........DVFLIGPPGPLRRSIAM.QYLELT..............KREVEYIALSR.DT..TETDLKQRREIR\n"
    ]

    transfer_dict = find_and_map_annots(
        logger=logger,
        hmmalign_lines=mock_lines,
        annotations=annotations_content_binding_fixture,
        good_eco_codes=[]
    )

    assert not transfer_dict

    # Verify the logger captured the correct error message
    logger.error.assert_called_once_with(
        "---> ERROR --- FIND_AND_MAP --- No target sequences found in hmmalign lines!"
    )


###T map_and_filter_annot_pos
def test_map_and_filter_with_paired_positions(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_disulfid_pair, transfer_dict,
    get_annotation_dict, annotations_content_disulfid_fixture
):
    """Unit test for map_and_filter_annot_pos with paired positions."""
    annotation_dict = get_annotation_dict(annotations_content_disulfid_fixture, "246")

    with patch('transfer_annotations.validate_paired_annotations') as mock_validate_paired_annotations:
        # Configure the mock to return expected values
        mock_validate_paired_annotations.return_value = (True, {'dummy': 'result'})

        # Call the function under test
        result = map_and_filter_annot_pos(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            paired_annot_pos_str=paired_annot_pos_str_mock,
            caller_target_pos_str=caller_target_pos_str_mock,
            paired_annotation_dict=annotation_dict
        )

        # Assert that validate_paired_annotations was called with expected arguments
        mock_validate_paired_annotations.assert_called_once_with(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            paired_annotation_dict=annotation_dict,
            entry_annotations=entry_annotations_disulfid_pair,
            counter_target_pos=ANY,  # Since this is set inside the function
            counter_annot_pos=ANY,  # Since this is set inside the function
            paired_annot_pos_str=paired_annot_pos_str_mock,
            caller_target_pos_str=caller_target_pos_str_mock,
        )

        # Verify the result
        assert result == (True, {'dummy': 'result'})

def test_map_and_filter_without_paired_positions(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_binding_only, transfer_dict
):
    """Unit test for map_and_filter_annot_pos without paired positions."""

    with patch('transfer_annotations.validate_annotations') as mock_validate_annotations:
        # Call the function under test
        result = map_and_filter_annot_pos(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            # No paired position parameters
        )

        # Assert that validate_annotations was called with expected arguments
        mock_validate_annotations.assert_called_once_with(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            counter_target_pos=ANY,  # These are set within the function
            counter_annot_pos=ANY,
        )

        # Since the function doesn't return anything in this case
        assert result is None

###T validate_annotations
def test_validate_annotations_process_annotation_parameters(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_binding_only, transfer_dict
):
    """Test that process_annotation is called with correct parameters"""

    with patch('transfer_annotations.process_annotation') as mock_process:
        validate_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
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
                good_eco_codes=good_eco_codes_all,
                entry_mnemo_name=entry_mnemo_name_mock,
                target_name=target_name_mock,
                target_hit_start=target_hit_start_mock,
                target_hit_end=target_hit_end_mock,
                annotation_dict=entry_annotations_binding_only["201"][0],
                entry_annotations=entry_annotations_binding_only,
                transfer_dict=transfer_dict,
                target_sequence=target_sequence,
                offset_start=offset_start_mock,
                offset_end=offset_end_mock,
                annot_sequence=annot_sequence,
                counter_target_pos_str='329',
                counter_annot_pos_str="201",
                processed_annotations=set(),
                annotation_key=('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '201', '329', 'BINDING')
            ),
            call(
                res_hit=True,
                logger=logger,
                good_eco_codes=good_eco_codes_all,
                entry_mnemo_name=entry_mnemo_name_mock,
                target_name=target_name_mock,
                target_hit_start=target_hit_start_mock,
                target_hit_end=target_hit_end_mock,
                annotation_dict=entry_annotations_binding_only["202"][0],
                entry_annotations=entry_annotations_binding_only,
                transfer_dict=transfer_dict,
                target_sequence=target_sequence,
                offset_start=offset_start_mock,
                offset_end=offset_end_mock,
                annot_sequence=annot_sequence,
                counter_target_pos_str='330',
                counter_annot_pos_str="202",
                processed_annotations=set(),
                annotation_key=('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '202', '330', 'BINDING')
            )
        ]

        mock_process.assert_has_calls(expected_calls, any_order=False)

def test_validate_annotations_multiple_annotations(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    transfer_dict
):
    """Test case where a single position has multiple annotations."""
    entry_annotations = {
        "204": [
            {"type": "BINDING", "description": "Interacts with ATP", "evidence": "ECO:0000255"},
            {"type": "ACT_SITE", "description": "Catalytic activity", "evidence": "ECO:0000269"}
        ]
    }

    with patch('transfer_annotations.process_annotation') as mock_process:
        validate_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=204,
            offset_end=204,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_annotations=entry_annotations,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            counter_target_pos=counter_target_pos_mock,
            counter_annot_pos=counter_annot_pos_mock
        )

        # Ensure all annotations at the position are processed
        assert mock_process.call_count == 2

def test_validate_annotations_skip_processed(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_binding_only, transfer_dict
):
    """Test case 1.2 - Skip when annotation was already processed"""

    processed_annotations = {
        (
            "MCRB_ECOLI",  # entry_mnemo_name
            "sp|Q9NU22|MDN1_HUMAN",  # target_name
            "201",  # counter_annot_pos_str
            "329",  # counter_target_pos_str
            "BINDING"  # anno_type
        ),
        (
            "MCRB_ECOLI",
            "sp|Q9NU22|MDN1_HUMAN",
            "202",
            "330",
            "BINDING"
        )
    }

    with patch('transfer_annotations.process_annotation') as mock_process:
        validate_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
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
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_binding_only, transfer_dict
):
    """Test case 2 - No annotations in offset range (201 and 202 don't fit in 500-600)"""
    with patch('transfer_annotations.process_annotation') as mock_process:
        validate_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=500,
            offset_end=600,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            counter_target_pos=counter_target_pos_mock,
            counter_annot_pos=counter_annot_pos_mock
        )

        mock_process.assert_not_called()
        assert transfer_dict == {}

def test_validate_annotations_sequence_end_reached(
    logger, good_eco_codes_all,
    entry_annotations_binding_only, transfer_dict
):
    """Test case where sequence reaches end condition (line 633)"""
    # Use sequences guaranteed to hit end condition
    test_short_sequence = ".....ABC"  # Short sequence to ensure we hit end

    validate_annotations(
        logger=logger,
        good_eco_codes=good_eco_codes_all,
        target_sequence=test_short_sequence,
        target_name=target_name_mock,
        target_hit_start=1,
        target_hit_end=3,  # Will trigger end condition when counter_target_pos == target_hit_end
        offset_start=1,
        offset_end=3,
        annot_sequence=test_short_sequence,
        entry_mnemo_name=entry_mnemo_name_mock,
        entry_annotations=entry_annotations_binding_only,
        transfer_dict=transfer_dict,
        processed_annotations=set(),
        counter_target_pos=None,
        counter_annot_pos=None
    )

def test_validate_annotations_logging(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_binding_only, transfer_dict
):
    """Test that error logs are written when process_annotation fails."""
    mock_error = Exception("Test error")

    with patch('transfer_annotations.process_annotation', side_effect=mock_error):
        try:
            validate_annotations(
                logger=logger,
                good_eco_codes=good_eco_codes_all,
                target_sequence=target_sequence,
                target_name=target_name_mock,
                target_hit_start=target_hit_start_mock,
                target_hit_end=target_hit_end_mock,
                offset_start=201,
                offset_end=202,
                annot_sequence=annot_sequence,
                entry_mnemo_name=entry_mnemo_name_mock,
                entry_annotations=entry_annotations_binding_only,
                transfer_dict=transfer_dict,
                processed_annotations=set(),
                counter_target_pos=counter_target_pos_mock,
                counter_annot_pos=counter_annot_pos_mock
            )
        except Exception:
            pass

        # Verify error log was written
        logger.error.assert_any_call("---> ERROR --- VAL_ANNOTS --- Error in validate_annotations: Test error")


###T process_annotation
def test_process_annotation_basic_single(
    logger, good_eco_codes_all, transfer_dict,
    target_sequence, annot_sequence,
    entry_annotations_binding_only, get_annotation_dict,
    annotations_content_binding_fixture
):
    """Test processing of a single BINDING annotation"""
    annotation_dict = get_annotation_dict(annotations_content_binding_fixture, "201")
    mock_make_anno_total = {
        'annotation': annotation_dict,
        'anno_type': 'BINDING',
        'anno_id': 'BINDING | Interacts with GTP',
        'anno_total': {
            'type': 'BINDING',
            'description': 'Interacts with GTP',
            'count': 1,
            'evidence': 'ECO:0000255',
            'ligand_id': 'ChEBI:CHEBI:37565',
            'target_position': '329',
            'annot_position': '201'
        },
        'paired_annot_pos_str': None
    }
    processed_annotations = set()

    with patch('transfer_annotations.make_anno_total_dict', return_value=mock_make_anno_total) as mock_make_total, \
         patch('transfer_annotations.add_to_transfer_dict') as mock_add_dict:
        process_annotation(
            res_hit=True,
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            annotation_dict=annotation_dict,
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            counter_target_pos_str="329",
            counter_annot_pos_str="201",
            processed_annotations=processed_annotations,
            annotation_key=('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '201', '329', 'BINDING')
        )

        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock,
            annotation_dict=annotation_dict,
            counter_target_pos_str='329',
            counter_annot_pos_str="201",
            logger=logger,
            entry_annotations=entry_annotations_binding_only
        )

        mock_add_dict.assert_called_once_with(
            hit=True,
            logger=logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock,
            anno_id="BINDING | Interacts with GTP",
            anno_total=mock_make_anno_total["anno_total"],
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_primary_accession=entry_primary_accession_mock,
        )

        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '201', '329', 'BINDING') in processed_annotations

def test_process_annotation_paired_failure(
    logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair,
    target_sequence, annot_sequence,
    mock_make_anno_total_disulfid_return_205_P15005,
    mock_map_filter_disulfid_return,
    annotation_dict_205, annotation_dict_246
):
    """
    Test processing of paired DISULFID annotations
    where paired pair fails to match, different residues
    """

    # Modify target sequence to ensure failure at paired position
    target_sequence_mod_to_fail_72 = list(target_sequence)
    target_sequence_mod_to_fail_72[72] = "R"
    target_sequence_mod_to_fail_72 = "".join(target_sequence_mod_to_fail_72)
    _, mock_paired_result_dict  = mock_map_filter_disulfid_return
    mock_map_filter_disulfid_fail_return = (False, mock_paired_result_dict)

    processed_annotations = set()

    with patch('transfer_annotations.make_anno_total_dict', return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch('transfer_annotations.map_and_filter_annot_pos', return_value=mock_map_filter_disulfid_fail_return) as mock_map_filter, \
         patch('transfer_annotations.add_to_transfer_dict') as mock_add_transfer:
        process_annotation(
            res_hit=True,
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            annotation_dict=annotation_dict_205,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence_mod_to_fail_72,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID")
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock,
            annotation_dict=annotation_dict_205,
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair
        )

        # Verify map_and_filter_annot_pos was called correctly
        mock_map_filter.assert_called_once_with(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_mod_to_fail_72,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            paired_annot_pos_str="246",
            caller_target_pos_str="333",
            paired_annotation_dict=annotation_dict_246
        )

        # Verify remove_failed_annotations was called correctly
        mock_add_transfer.assert_called_once_with(
            hit=True,
            logger=logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock,
            anno_id='DISULFID | Intrachain (with C-246); in linked form',
            anno_total=mock_make_anno_total_disulfid_return_205_P15005['anno_total'],
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_primary_accession=entry_primary_accession_mock,
            paired_position_res_hit=False,
            paired_anno_id='DISULFID | Intrachain (with C-205); in linked form',
            paired_anno_total=mock_map_filter_disulfid_return[1]['anno_total']
        )

        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '205', '333', 'DISULFID') in processed_annotations
        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '246', '373', 'DISULFID') in processed_annotations


def test_process_annotation_paired_success(
    logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair, annotation_dict_205,
    annotation_dict_246, mock_make_anno_total_disulfid_return_205_P15005,
    mock_map_filter_disulfid_return,
    target_sequence, annot_sequence
):
    """Test processing of paired DISULFID annotations where paired pair succeeds"""

    processed_annotations = set()

    with patch('transfer_annotations.make_anno_total_dict', return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch('transfer_annotations.map_and_filter_annot_pos', return_value=mock_map_filter_disulfid_return) as mock_map_filter, \
         patch('transfer_annotations.add_to_transfer_dict') as mock_add_transfer:
        process_annotation(
            res_hit=True,
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            annotation_dict=annotation_dict_205,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            counter_target_pos_str="333",
            counter_annot_pos_str="205",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID")
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock,
            annotation_dict=annotation_dict_205,
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair
        )

        # Verify map_and_filter_annot_pos was called correctly
        mock_map_filter.assert_called_once_with(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            paired_annot_pos_str="246",
            caller_target_pos_str="333",
            paired_annotation_dict=annotation_dict_246
        )

        # Verify remove_failed_annotations was called correctly
        mock_add_transfer.assert_called_once_with(
            hit=True,
            logger=logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock,
            anno_id='DISULFID | Intrachain (with C-246); in linked form',
            anno_total=mock_make_anno_total_disulfid_return_205_P15005['anno_total'],
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_primary_accession=entry_primary_accession_mock,
            paired_position_res_hit=True,
            paired_anno_id='DISULFID | Intrachain (with C-205); in linked form',
            paired_anno_total=mock_map_filter_disulfid_return[1]['anno_total']
        )

        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '205', '333', 'DISULFID') in processed_annotations
        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '246', '373', 'DISULFID') in processed_annotations


def test_process_annotation_paired_no_paired_dict(
    logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair, annotation_dict_205,
    mock_make_anno_total_disulfid_return_205_P15005,
    target_sequence, annot_sequence
):
    """Test processing of paired DISULFID annotations where no paired annotation dict is found"""

    # Create a modified entry_annotations without the paired annotation
    modified_entry_annotations = {
        "205": entry_annotations_disulfid_pair["205"],
        "246": []  # Empty list to trigger the else branch
    }

    processed_annotations = set()

    with patch('transfer_annotations.make_anno_total_dict', return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch('transfer_annotations.add_to_transfer_dict') as mock_add_transfer:
        process_annotation(
            res_hit=True,
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            annotation_dict=annotation_dict_205,
            entry_annotations=modified_entry_annotations,
            transfer_dict=transfer_dict,
            target_sequence=target_sequence,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID")
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock,
            annotation_dict=annotation_dict_205,
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            logger=logger,
            entry_annotations=modified_entry_annotations
        )

        # Verify remove_failed_annotations was called correctly
        mock_add_transfer.assert_called_once_with(
            hit=True,
            logger=logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock,
            anno_id='DISULFID | Intrachain (with C-246); in linked form',
            anno_total=mock_make_anno_total_disulfid_return_205_P15005['anno_total'],
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_primary_accession=entry_primary_accession_mock,
            paired_position_res_hit=False,
            paired_anno_id=None,
            paired_anno_total=None
        )

        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '205', '333', 'DISULFID') in processed_annotations



###T make_anno_total_dict
def test_make_anno_total_dict_basic(
    logger, good_eco_codes_all, entry_annotations_disulfid_pair,
    annotation_dict_205
):
    """Test basic case with valid DISULFID annotation"""
    result = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock,
        annotation_dict=annotation_dict_205,
        counter_target_pos_str="333",
        counter_annot_pos_str="205",
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair
    )

    expected = {
        'annotation': {
            **annotation_dict_205,
        },
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-246); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-246); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678',
            'target_position': '333',
            'annot_position': '205',
        },
        'paired_annot_pos_str': '246'
    }

    assert result == expected

def test_make_anno_total_dict_multiple_evidence(
    logger, good_eco_codes_all, entry_annotations_disulfid_pair,
    annotation_dict_205
):
    """Test case with multiple evidences in the string"""
    annotation_dict_205["evidence"] = ["ECO:0000269|PubMed:12345678, ECO:0007744|PDB:1ABC"]

    result_single_pass = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock,
        annotation_dict=annotation_dict_205,
        counter_target_pos_str='333',
        counter_annot_pos_str="205",
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair
    )

    expected_single_pass = {
        'annotation': {
            **annotation_dict_205,
        },
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-246); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-246); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678',
            'target_position': '333',
            'annot_position': '205',
        },
        'paired_annot_pos_str': '246'
    }

    expected_multiple_pass = copy.deepcopy(expected_single_pass)
    expected_multiple_pass["anno_total"]["evidence"] = "ECO:0000269|PubMed:12345678, ECO:0007744|PDB:1ABC"

    assert result_single_pass == expected_single_pass

    good_eco_codes_all_modded = good_eco_codes_all + ["ECO:0007744"]

    result_multiple_pass = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all_modded,
        entry_mnemo_name=entry_mnemo_name_mock,
        annotation_dict=annotation_dict_205,
        counter_target_pos_str='333',
        counter_annot_pos_str="205",
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair
    )

    assert result_multiple_pass == expected_multiple_pass

def test_make_anno_total_dict_eco_filtering(
    logger, entry_annotations_binding_only,
    get_annotation_dict, annotations_content_binding_fixture
):
    """Test ECO code filtering with BINDING annotation"""
    annotation_dict = get_annotation_dict(annotations_content_binding_fixture, "201")

    result = make_anno_total_dict(
        good_eco_codes=["ECO:0000269"],  # Exclude all other codes, including current: ECO:0000255
        entry_mnemo_name=entry_mnemo_name_mock,
        annotation_dict=annotation_dict,
        counter_target_pos_str='329',
        counter_annot_pos_str="201",
        logger=logger,
        entry_annotations=entry_annotations_binding_only
    )

    assert result['anno_total'] is None

def test_make_anno_total_dict_binding_type(
    logger, good_eco_codes_all, entry_annotations_binding_only,
    get_annotation_dict, annotations_content_binding_fixture
):
    """Test BINDING type annotation with additional keys"""
    annotation_dict = get_annotation_dict(annotations_content_binding_fixture, "201")

    result = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock,
        annotation_dict=annotation_dict,
        counter_target_pos_str='329',
        counter_annot_pos_str="201",
        logger=logger,
        entry_annotations=entry_annotations_binding_only
    )

    expected = {
        'annotation': {
            **annotation_dict,
        },
        'anno_type': 'BINDING',
        'anno_id': 'BINDING | Interacts with GTP',
        'anno_total': {
            'type': 'BINDING',
            'description': 'Interacts with GTP',
            'count': 1,
            'evidence': 'ECO:0000255',
            'ligand_id': 'ChEBI:CHEBI:37565',
            'target_position': '329',
            'annot_position': '201'
        },
        'paired_annot_pos_str': None
    }

    assert result['annotation'] is annotation_dict
    assert result['anno_type'] == 'BINDING'
    assert result['anno_id'] == 'BINDING | Interacts with GTP'
    assert result['paired_annot_pos_str'] is None
    assert result == expected

def test_make_anno_total_dict_with_caller_target_pos(
    logger, good_eco_codes_all, entry_annotations_disulfid_pair,
    annotation_dict_246
):
    """Test with caller_target_pos provided (paired annotation case)"""
    result = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock,
        annotation_dict=annotation_dict_246,
        counter_target_pos_str='373',
        counter_annot_pos_str="246",
        caller_target_pos_str='333',
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair
    )

    expected = {
        'annotation': {
            **annotation_dict_246,
        },
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-205); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-205); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678',
            'target_position': '373',
            'annot_position': '246',
            'paired_target_position': '333'
        },
        'paired_annot_pos_str': '205'
    }

    assert result['annotation'] is annotation_dict_246
    assert result['anno_type'] == 'DISULFID'
    assert result['anno_id'] == 'DISULFID | Intrachain (with C-205); in linked form'
    assert result['paired_annot_pos_str'] == '205'
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
        entry_mnemo_name=entry_mnemo_name_mock,
        annotation_dict=annotation_without_evidence,
        counter_target_pos_str='100',
        counter_annot_pos_str="201",
        logger=logger,
        entry_annotations={"201": [annotation_without_evidence]}
    )

    assert result['annotation'] is annotation_without_evidence
    assert result['anno_type'] is 'BINDING'
    assert result['anno_id'] == 'BINDING | Test binding'
    assert result['paired_annot_pos_str'] is None
    assert result['anno_total'] is None


###T validate_paired_annotations
def test_validate_paired_annotations_valid_case(
    logger,
    good_eco_codes_all,
    target_sequence,
    annot_sequence,
    entry_annotations_disulfid_pair,
    annotation_dict_246
):
    """
    Test the case where paired annotations are valid and successfully validated.
    """
    # Mock the output of map_and_filter_annot_pos for paired validation
    mock_paired_result_dict = {
        "annotation": {
            "type": "DISULFID",
            "description": "Intrachain (with C-246); in linked form",
            "evidence": "ECO:0000269|PubMed:12345678",
            "entry": "P15005",
            "aminoacid": "C",
            "paired_position": "205",
        },
        "anno_type": "DISULFID",
        "anno_id": "DISULFID | Intrachain (with C-246); in linked form",
        "anno_total": {
            "type": "DISULFID",
            "description": "Intrachain (with C-246); in linked form",
            "count": 1,
            "evidence": "ECO:0000269|PubMed:12345678",
            "target_position": "373",
            "annot_position": "246",
        },
        "paired_annot_pos_str": "205",
    }

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_paired_result_dict) as mock_make_total:
        result_tuple = validate_paired_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            paired_annotation_dict=annotation_dict_246,
            entry_annotations=entry_annotations_disulfid_pair,
            counter_target_pos=counter_target_pos_mock, # Always None
            counter_annot_pos=counter_annot_pos_mock, # Always None
            paired_annot_pos_str=paired_annot_pos_str_mock, # '246'
            caller_target_pos_str=caller_target_pos_str_mock, # '333'
        )

        # Verify the make_anno_total_dict call
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name="MCRB_ECOLI",
            annotation_dict=annotation_dict_246,
            counter_target_pos_str='373', # 246 = 373 in target sequence numbering/col
            counter_annot_pos_str=paired_annot_pos_str_mock, # '246'
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair,
            caller_target_pos_str=caller_target_pos_str_mock,
        )

        # Assert that the result matches the mock data
        assert result_tuple == (True, mock_paired_result_dict)

#RESOLVER
def test_validate_paired_annotations_invalid_pair(
    logger,
    good_eco_codes_all,
    target_sequence,
    annot_sequence,
    entry_annotations_disulfid_pair,
    annotation_dict_246
):
    """
    Test case where the paired annotation validation fails.
    """

    # Modify target sequence to ensure failure at paired position
    target_sequence_mod_to_fail_72 = list(target_sequence)
    target_sequence_mod_to_fail_72[72] = "R"
    target_sequence_mod_to_fail_72 = "".join(target_sequence_mod_to_fail_72)

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
        },
        "paired_annot_pos_str": "205",
    }

    mock_result_tuple = (False, mock_paired_result_dict)

    with patch("transfer_annotations.make_anno_total_dict", return_value=mock_paired_result_dict) as mock_make_total:
        result_tuple = validate_paired_annotations(
            logger=logger,
            good_eco_codes=good_eco_codes_all,
            target_sequence=target_sequence_mod_to_fail_72,
            target_name=target_name_mock,
            target_hit_start=target_hit_start_mock,
            target_hit_end=target_hit_end_mock,
            offset_start=offset_start_mock,
            offset_end=offset_end_mock,
            annot_sequence=annot_sequence,
            entry_mnemo_name=entry_mnemo_name_mock,
            paired_annotation_dict=annotation_dict_246,
            entry_annotations=entry_annotations_disulfid_pair,
            counter_target_pos=counter_target_pos_mock, # Always None
            counter_annot_pos=counter_annot_pos_mock, # Always None
            paired_annot_pos_str=paired_annot_pos_str_mock, # '246'
            caller_target_pos_str=caller_target_pos_str_mock, # '333'
        )

        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock,
            annotation_dict=annotation_dict_246,
            counter_target_pos_str='373', # 246 = 373 in target sequence numbering/col
            counter_annot_pos_str=paired_annot_pos_str_mock, # '246'
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair,
            caller_target_pos_str=caller_target_pos_str_mock,
        )

        assert result_tuple == mock_result_tuple

def test_validate_paired_annotations_missing_paired_position(
    logger,
    good_eco_codes_all,
    target_sequence,
    annot_sequence,
    entry_annotations_disulfid_pair,
    annotation_dict_246
):
    """
    Test case where the paired position is missing in the annotations.
    """
    # Set paired position to a non-existent value
    annotation_dict_246["paired_position"] = "999"

    result = validate_paired_annotations(
        logger=logger,
        good_eco_codes=good_eco_codes_all,
        target_sequence=target_sequence,
        target_name=target_name_mock,
        target_hit_start=target_hit_start_mock,
        target_hit_end=target_hit_end_mock,
        offset_start=offset_start_mock,
        offset_end=offset_end_mock,
        annot_sequence=annot_sequence,
        entry_mnemo_name=entry_mnemo_name_mock,
        paired_annotation_dict=annotation_dict_246,
        entry_annotations=entry_annotations_disulfid_pair,
        counter_target_pos=counter_target_pos_mock,
        counter_annot_pos=counter_annot_pos_mock,
        paired_annot_pos_str="999",
        caller_target_pos_str=caller_target_pos_str_mock,
    )

    assert result == (False, {})




###T add_to_transfer_dict
def test_add_to_transfer_dict_disulfid_single_first_addition(
    logger, transfer_dict,
    mock_make_anno_total_disulfid_return_205_P15005,
    mock_anno_total_disulfid
):
    """Test add_to_transfer_dict with successful annotation addition"""

    anno_id = mock_make_anno_total_disulfid_return_205_P15005['anno_id']
    anno_total = mock_anno_total_disulfid

    add_to_transfer_dict(
        hit=True,
        logger=logger,
        transfer_dict=transfer_dict,
        target_name=target_name_mock,
        anno_id=anno_id,
        anno_total=anno_total,
        entry_mnemo_name=entry_mnemo_name_mock,
        entry_primary_accession=entry_primary_accession_mock,
    )

    assert target_name_mock in transfer_dict
    assert 'match' in transfer_dict[target_name_mock]
    assert anno_total['target_position'] in transfer_dict[target_name_mock]['match']
    assert anno_id in transfer_dict[target_name_mock]['match'][anno_total['target_position']]
    assert transfer_dict[target_name_mock]['match'][anno_total['target_position']][anno_id]['essentials']['type'] == 'DISULFID'

def test_add_to_transfer_dict_paired_disulfide(
    logger, transfer_dict,
    mock_make_anno_total_disulfid_return_205_P15005,
    mock_make_anno_total_disulfid_return_246_P15005,
    mock_anno_total_disulfid,
    mock_paired_anno_total_disulfid
):
    """Test paired disulfide annotation processing"""
    anno_id = mock_make_anno_total_disulfid_return_205_P15005['anno_id']
    anno_total = mock_anno_total_disulfid.copy()
    paired_anno_total = mock_paired_anno_total_disulfid.copy()
    paired_anno_id = mock_make_anno_total_disulfid_return_246_P15005['anno_id']

    add_to_transfer_dict(
        hit=True,
        logger=logger,
        transfer_dict=transfer_dict,
        target_name=target_name_mock,
        anno_id=anno_id,
        anno_total=anno_total,
        entry_mnemo_name=entry_mnemo_name_mock,
        entry_primary_accession=entry_primary_accession_mock,
        paired_position_res_hit=True,
        paired_anno_id=paired_anno_id,
        paired_anno_total=paired_anno_total
    )

    # Assert main position
    assert '333' in transfer_dict[target_name_mock]['match']
    assert anno_id in transfer_dict[target_name_mock]['match']['333']

    # Assert paired position
    assert '373' in transfer_dict[target_name_mock]['match']
    assert paired_anno_id in transfer_dict[target_name_mock]['match']['373']

    # Assert cross-references
    assert '373' in transfer_dict[target_name_mock]['match']['333'][anno_id]['paired_position']
    assert '333' in transfer_dict[target_name_mock]['match']['373'][paired_anno_id]['paired_position']





def test_add_to_transfer_dict_paired_repeated_anno_id(
    logger, transfer_dict_populated,
    mock_make_anno_total_disulfid_return_205_P00750,
    mock_make_anno_total_disulfid_return_246_P00750,
):
    """Test paired disulfide annotation processing"""
    anno_id = mock_make_anno_total_disulfid_return_205_P00750['anno_id']
    anno_total = mock_make_anno_total_disulfid_return_205_P00750['anno_total']
    paired_anno_id = mock_make_anno_total_disulfid_return_246_P00750['anno_id']
    paired_anno_total = mock_make_anno_total_disulfid_return_246_P00750['anno_total']

    add_to_transfer_dict(
        hit=True,
        logger=logger,
        transfer_dict=transfer_dict_populated,
        target_name=target_name_mock,
        anno_id=anno_id,
        anno_total=anno_total,
        entry_mnemo_name=entry_mnemo_name_mock_repeated,
        entry_primary_accession=entry_primary_accession_mock_repeated,
        paired_position_res_hit=True,
        paired_anno_id=paired_anno_id,
        paired_anno_total=paired_anno_total
    )

    print(json.dumps(transfer_dict_populated, indent=4))

    # Helper function to access nested dictionary paths
    def get_annotation(pos):
        return transfer_dict_populated[target_name_mock]['match'][pos]

    # Test data
    positions = {
        '333': {'partner': '373', 'annot_pos': '205'},
        '373': {'partner': '333', 'annot_pos': '246'}
    }

    for pos, data in positions.items():
        anno = get_annotation(pos)['DISULFID | Intrachain (with C-246); in linked form' if pos == '333' else 'DISULFID | Intrachain (with C-205); in linked form']

        # Count assertions
        assert anno['essentials']['count'] == 2
        assert anno['paired_position'][data['partner']]['count'] == 2
        assert anno['additional_keys']['annot_position'][data['annot_pos']]['count'] == 2

        # Evidence assertions
        for eco_code in ['ECO:0000250', 'ECO:0000269|PubMed:12345678']:
            assert eco_code in anno['evidence']
            if eco_code == 'ECO:0000250':
                evidence = anno['evidence'][eco_code]
                assert evidence['rep_mnemo_name'] == 'TPA_HUMAN'
                assert evidence['rep_primary_accession'] == 'P00720'
                assert evidence['count'] == 1

    # Verify no misses
    assert transfer_dict_populated[target_name_mock]['miss'] == {}

def test_add_to_transfer_dict_hit_miss_pair(
    logger, transfer_dict,
    mock_make_anno_total_disulfid_return_205_P15005,
    mock_make_anno_total_disulfid_return_246_P15005,
    mock_anno_total_disulfid,
    mock_paired_anno_total_disulfid
):
    """Test disulfide pair where one residue hits and other misses"""
    anno_id = mock_make_anno_total_disulfid_return_205_P15005['anno_id']
    anno_total = mock_anno_total_disulfid.copy()
    paired_anno_total = mock_paired_anno_total_disulfid.copy()
    paired_anno_id = mock_make_anno_total_disulfid_return_246_P15005['anno_id']

    add_to_transfer_dict(
        hit=True,
        logger=logger,
        transfer_dict=transfer_dict,
        target_name=target_name_mock,
        anno_id=anno_id,
        anno_total=anno_total,
        entry_mnemo_name=entry_mnemo_name_mock,
        entry_primary_accession=entry_primary_accession_mock,
        paired_position_res_hit=False,
        paired_anno_id=paired_anno_id,
        paired_anno_total=paired_anno_total
    )

    # Assert hit position
    assert '333' in transfer_dict[target_name_mock]['match']

    # Assert missed position
    assert '373' in transfer_dict[target_name_mock]['miss']



def test_add_to_transfer_dict_disulfid_single_mocked_helper(
    logger, transfer_dict,
    mock_make_anno_total_disulfid_return_205_P15005,
    mock_anno_total_disulfid
):
    """Test add_to_transfer_dict with successful annotation addition"""

    anno_id = mock_make_anno_total_disulfid_return_205_P15005['anno_id']
    anno_total = mock_anno_total_disulfid

    # Call the function
    with patch('transfer_annotations._add_single_annotation') as mock_add_single:
        add_to_transfer_dict(
            hit=True,
            logger=logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock,
            anno_id=anno_id,
            anno_total=anno_total,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_primary_accession=entry_primary_accession_mock,
        )

        # Assertions to verify correct addition to transfer_dict
        mock_add_single.assert_called_once_with(
            hit=True,
            transfer_dict=transfer_dict,
            target_name=target_name_mock,
            anno_id=anno_id,
            anno_total=anno_total,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_primary_accession=entry_primary_accession_mock,
            additional_keys={}
        )

def test_add_to_transfer_dict_disulfid_paired_mocked_helper(
    logger, transfer_dict,
    mock_make_anno_total_disulfid_return_205_P15005,
    mock_make_anno_total_disulfid_return_246_P15005,
    mock_anno_total_disulfid,
    mock_paired_anno_total_disulfid
):
    """
    Test add_to_transfer_dict with successful annotation pair
    addition. Intent on seeing that expected calls went through.
    """

    anno_id = mock_make_anno_total_disulfid_return_205_P15005['anno_id']
    anno_total = mock_anno_total_disulfid.copy()
    paired_anno_total = mock_paired_anno_total_disulfid.copy()
    paired_anno_id = mock_make_anno_total_disulfid_return_246_P15005['anno_id']

    # Call the function
    with patch('transfer_annotations._add_single_annotation') as mock_add_single:
        add_to_transfer_dict(
            hit=True,
            logger=logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock,
            anno_id=anno_id,
            anno_total=anno_total,
            entry_mnemo_name=entry_mnemo_name_mock,
            entry_primary_accession=entry_primary_accession_mock,
            paired_position_res_hit=True,
            paired_anno_id=paired_anno_id,
            paired_anno_total=paired_anno_total
        )

        expected_calls = [
            call(
                hit=True,
                transfer_dict=transfer_dict,
                target_name=target_name_mock,
                anno_id=anno_id,
                anno_total=anno_total,
                entry_mnemo_name=entry_mnemo_name_mock,
                entry_primary_accession=entry_primary_accession_mock,
                additional_keys={}
            ),
            call(
                hit=True,
                transfer_dict=transfer_dict,
                target_name=target_name_mock,
                anno_id=paired_anno_id,
                anno_total=paired_anno_total,
                entry_mnemo_name=entry_mnemo_name_mock,
                entry_primary_accession=entry_primary_accession_mock,
                additional_keys={}
            )
        ]

        mock_add_single.assert_has_calls(expected_calls, any_order=False)
        assert mock_add_single.call_count == 2

###T write_report

def test_write_report_with_data(
    tmp_path, logger,
    transfer_dict_populated
    ):
    """Test write_report with actual annotation data"""

    with patch('builtins.open', mock_open()) as mock_file:
        write_report(logger, transfer_dict_populated, "PF00244", str(tmp_path))

        mock_file.assert_called_once_with(ANY, 'w', encoding="utf-8")
        written_content = mock_file().write.call_args[0][0]
        assert json.loads(written_content) == transfer_dict_populated['sp|Q9NU22|MDN1_HUMAN']
        assert logger.debug.call_count == 2

def test_write_report_empty_data(tmp_path, logger):
    """Test write_report with empty transfer dict"""
    transfer_dict = {"sp|Q9NU22|MDN1_HUMAN": {}}

    with patch('builtins.open', mock_open()) as mock_file:
        write_report(logger, transfer_dict, "PF00244", str(tmp_path))

        mock_file.assert_called_once_with(ANY, 'w', encoding="utf-8")
        written_content = mock_file().write.call_args[0][0]
        assert json.loads(written_content) == {"Annotations": "None"}
        assert logger.debug.call_count == 2

def test_write_report_multiple_targets(tmp_path, logger):
    """Test write_report with multiple target sequences"""
    transfer_dict = {
        "TEST1": {"match": {}},
        "TEST2": {"match": {}}
    }

    with patch('builtins.open', mock_open()) as mock_file:
        write_report(logger, transfer_dict, "PF00244", str(tmp_path))

        assert mock_file.call_count == 2
        assert logger.debug.call_count == 4
        calls = mock_file.call_args_list
        assert all(call_args[1] == {'encoding': 'utf-8'} for call_args in calls)

def test_write_report_pipe_in_target_name(tmp_path, logger):
    """Test write_report handles target names with pipes"""
    transfer_dict = {
        "TEST|A": {"match": {}}
    }

    with patch('builtins.open', mock_open()) as mock_file:
        write_report(logger, transfer_dict, "PF00244", str(tmp_path))

        mock_file.assert_called_once_with(ANY, 'w', encoding="utf-8")
        # Verify target name was transformed
        called_path = mock_file.call_args[0][0]
        assert "TEST-A" in called_path


###T main
def test_main_success(logger):
    """Test main function success path"""
    mock_args = Namespace(
        dom_align=hmmalign_result_mock,
        resource_dir=resource_dir_mock,
        output_dir=output_dir_mock,
        eco_codes=good_eco_codes_mock,
        log=log_filepath_mock
    )

    mock_transfer_dict = {"TEST1": {"match": {}}}

    with patch('transfer_annotations.parse_arguments', return_value=mock_args) as mock_parse, \
         patch('transfer_annotations.get_pfam_id_from_hmmalign_result', return_value="PF00244") as mock_get_pfam, \
         patch('transfer_annotations.read_files', return_value=([], {})) as mock_read, \
         patch('transfer_annotations.find_and_map_annots', return_value=mock_transfer_dict) as mock_find, \
         patch('transfer_annotations.write_report') as mock_write:

        main(logger)

        mock_parse.assert_called_once()
        mock_get_pfam.assert_called_once_with(hmmalign_result_mock)
        mock_read.assert_called_once()
        mock_find.assert_called_once()
        mock_write.assert_called_once_with(logger, mock_transfer_dict, "PF00244", output_dir_mock)
        logger.info.assert_any_call("---> MAIN --- Transfer Dict FILLED")

# Integration tests
def test_main_integration_binding_Q9NU22_MDN1_HUMAN(transfer_dict_success_binding_MDN1_HUMAN, annotations_content_binding_fixture):
    """Integration test for main function"""

    with TemporaryDirectory() as tmp_dir:
        # Setup test files and directories
        hmmalign_path = os.path.join(tmp_dir, "PF00244_hmmalign.sth")
        with open(hmmalign_path, 'w', encoding='latin-1') as f:
            f.write(minimal_hmmalign_content_mock)

        # Create full directory structure
        resource_dir = os.path.join(tmp_dir, "resources")
        pfam_dir = os.path.join(resource_dir, "PF00244")
        os.makedirs(pfam_dir)

        output_dir = os.path.join(tmp_dir, "output")
        # Create output directories for expected targets
        target_dirs = [
            os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN"),
            os.path.join(output_dir, "PF00244")
        ]
        for dir_path in [output_dir] + target_dirs:
            os.makedirs(dir_path, exist_ok=True)

        # Write annotations to proper path
        with open(os.path.join(pfam_dir, "annotations.json"), 'w', encoding='utf-8') as f:
            json.dump(annotations_content_binding_fixture, f)

        # Run main
        args = Namespace(
            dom_align=hmmalign_path,
            resource_dir=resource_dir,
            output_dir=output_dir,
            eco_codes=good_eco_codes_mock,
            log=os.path.join(tmp_dir, "test.log")
        )

        with patch('transfer_annotations.parse_arguments', return_value=args):
            logger = configure_logging(args.log)
            main(logger)

        assert os.path.exists(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF00244_report.json"))
        assert transfer_dict_success_binding_MDN1_HUMAN == json.load(open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF00244_report.json"), 'r', encoding='utf-8'))

def test_main_integration_disulfid_Q9NU22_MDN1_HUMAN(transfer_dict_success_disulfid_MDN1_HUMAN, annotations_content_disulfid_fixture):
    """Integration test for main function"""

    with TemporaryDirectory() as tmp_dir:
        # Setup test files and directories
        hmmalign_path = os.path.join(tmp_dir, "PF00244_hmmalign.sth")
        with open(hmmalign_path, 'w', encoding='latin-1') as f:
            f.write(minimal_hmmalign_content_mock)

        # Create full directory structure
        resource_dir = os.path.join(tmp_dir, "resources")
        pfam_dir = os.path.join(resource_dir, "PF00244")
        os.makedirs(pfam_dir)

        output_dir = os.path.join(tmp_dir, "output")
        # Create output directories for expected targets
        target_dirs = [
            os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN"),
            os.path.join(output_dir, "PF00244")
        ]
        for dir_path in [output_dir] + target_dirs:
            os.makedirs(dir_path, exist_ok=True)

        # Write annotations to proper path
        with open(os.path.join(pfam_dir, "annotations.json"), 'w', encoding='utf-8') as f:
            json.dump(annotations_content_disulfid_fixture, f)

        # Run main
        args = Namespace(
            dom_align=hmmalign_path,
            resource_dir=resource_dir,
            output_dir=output_dir,
            eco_codes=good_eco_codes_mock,
            log=os.path.join(tmp_dir, "test.log")
        )

        with patch('transfer_annotations.parse_arguments', return_value=args):
            logger = configure_logging(args.log)
            main(logger)

        assert os.path.exists(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF00244_report.json"))
        assert transfer_dict_success_disulfid_MDN1_HUMAN == json.load(open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF00244_report.json"), 'r', encoding='utf-8'))

def test_main_integration_all_types_H0YB80_HUMAN(transfer_dict_success_all_types_H0YB80_HUMAN, annotations_content_all_types_PF00244_fixture):
    """Integration test for main function"""

    with TemporaryDirectory() as tmp_dir:
        # Setup test files and directories
        hmmalign_path = os.path.join(tmp_dir, "PF00244_hmmalign.sth")
        with open(hmmalign_path, 'w', encoding='latin-1') as f:
            f.write(minimal_hmmalign_content_mock)

        # Create full directory structure
        resource_dir = os.path.join(tmp_dir, "resources")
        pfam_dir = os.path.join(resource_dir, "PF00244")
        os.makedirs(pfam_dir)

        output_dir = os.path.join(tmp_dir, "output")

        # Create output directories for expected targets
        target_dirs = [
            os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN"),
            os.path.join(output_dir, "PF00244")
        ]
        for dir_path in [output_dir] + target_dirs:
            os.makedirs(dir_path, exist_ok=True)

        # Write annotations to proper path
        with open(os.path.join(pfam_dir, "annotations.json"), 'w', encoding='utf-8') as f:
            json.dump(annotations_content_all_types_PF00244_fixture, f)

        # Run main
        args = Namespace(
            dom_align=hmmalign_path,
            resource_dir=resource_dir,
            output_dir=output_dir,
            eco_codes=good_eco_codes_mock,
            log=os.path.join(tmp_dir, "test.log")
        )

        with patch('transfer_annotations.parse_arguments', return_value=args):
            logger = configure_logging(args.log)
            main(logger)

        assert os.path.exists(os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN", "PF00244_report.json"))
        assert transfer_dict_success_all_types_H0YB80_HUMAN == json.load(open(os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN", "PF00244_report.json"), 'r', encoding='utf-8'))



        # with open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF00244_report.json")) as f:
        #     transfer_dict = json.load(f)
        #     print(json.dumps(transfer_dict, indent=4))