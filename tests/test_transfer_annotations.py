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
from unittest.mock import patch, ANY, call, mock_open, MagicMock
from tempfile import TemporaryDirectory
# Add the parent directory to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from transfer_annotations import (
    parse_arguments,
    configure_logging,
    get_pfam_id_from_hmmalign_result,
    get_annotation_filepath,
    read_files,
    iterate_aligned_sequences,
    find_and_map_annots,
    read_conservations_and_annotations,
    parse_go_annotations,
    gather_go_terms_for_target,
    get_alignment_sequences,
    populate_conservation,
    populate_go_data_for_annotations,
    cleanup_improve_transfer_dict,
    convert_sets_to_lists,
    write_reports,
    map_and_filter_annot_pos,
    add_to_transfer_dict,
    make_anno_total_dict,
    validate_annotations,
    process_annotation,
    validate_paired_annotations,
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

go_terms_mock = "/home/user/results/human/"

good_eco_codes_mock = ["ECO:0000269", "ECO:0000255", "ECO:0000313", "ECO:0007669"]

log_filepath_mock = "/home/user/logs/transfer_annotations.log"


### Fixtures
@pytest.fixture
def annotations_content_binding_fixture_Q9NU22_PF07728():
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
def annotations_content_disulfid_fixture_Q9NU22_PF07728():
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
        "143": 0.9853259562184268,
        "146": 0.8806005719733079,
        "149": 0.956244087038789,
        "267": 0.8021297672114909,
        "283": 0.959562479284057,
        "284": 0.8452762923351159
    }
}

@pytest.fixture
def conservations_content_H0YB80_PF00244():
    return {
    "C1E2K1_MICCC/19-248": {
        "20": 0.9802880970432145,
        "21": 0.8013100436681223,
        "22": 0.83632447954056,
        "23": 0.8230659025787965,
        "24": 0.8863473909935669,
        "25": 0.8404558404558404,
        "26": 0.8196605374823197,
        "27": 0.8921775898520085,
        "28": 0.9500351864883885,
        "29": 0.8625525946704067,
        "32": 0.9007936507936508,
        "36": 0.8561151079136691,
        "51": 0.9382239382239382,
        "54": 0.8360025624599615,
        "55": 0.9296225207933462,
        "56": 0.9481765834932822,
        "57": 0.8817891373801917,
        "58": 0.9156549520766774,
        "59": 0.8792332268370607,
        "60": 0.9060702875399361,
        "61": 0.8446291560102301,
        "62": 0.8378033205619413,
        "63": 0.8959795788130185,
        "64": 0.9527760051052967,
        "65": 0.8864795918367347,
        "68": 0.8354591836734694,
        "70": 0.8068833652007649,
        "71": 0.9700255102040817,
        "74": 0.8662420382165605,
        "75": 0.9292993630573249,
        "81": 0.8132646490663232,
        "83": 0.8134615384615385,
        "101": 0.8401768793430195,
        "106": 0.8635502210991788,
        "108": 0.9348513598987982,
        "109": 0.8949367088607595,
        "113": 0.9358322744599746,
        "124": 0.9197452229299363,
        "132": 0.8320561941251596,
        "135": 0.8855498721227621,
        "136": 0.9213051823416507,
        "137": 0.8509277031349968,
        "139": 0.9430217669654289,
        "140": 0.8856416772554002,
        "141": 0.8162093171665603,
        "142": 0.8750796685787126,
        "143": 0.9859783301465902,
        "144": 0.9315856777493606,
        "146": 0.9775928297055058,
        "147": 0.9814696485623003,
        "149": 0.8355640535372849,
        "150": 0.9681528662420382,
        "168": 0.9736180904522613,
        "171": 0.9656893325015595,
        "175": 0.8119873817034701,
        "182": 0.8391739674593242,
        "184": 0.9311639549436797,
        "186": 0.9043151969981238,
        "187": 0.9743107769423559,
        "188": 0.9191222570532915,
        "189": 0.9228356336260979,
        "190": 0.8745294855708908,
        "191": 0.9780288763339611,
        "192": 0.964824120603015,
        "194": 0.9050911376492772,
        "195": 0.9376966645689113,
        "196": 0.9346323067253299,
        "198": 0.8826498422712934,
        "199": 0.9215189873417722,
        "200": 0.833859759949463,
        "207": 0.9721695129664769,
        "208": 0.8228969006957622,
        "211": 0.9375394321766561,
        "214": 0.8852564102564102,
        "215": 0.87001287001287,
        "216": 0.8009020618556701,
        "218": 0.9090322580645162,
        "222": 0.8463047743623283,
        "223": 0.8356344510190664,
        "227": 0.8168187744458931,
        "230": 0.8641732283464567,
        "231": 0.8045901639344263,
        "232": 0.8844386080105056,
        "233": 0.8355263157894737,
        "234": 0.8320158102766798,
        "235": 0.8148880105401844,
        "236": 0.9163372859025033,
        "237": 0.8081740276862228,
        "238": 0.8982826948480845,
        "239": 0.9404367968232958,
        "240": 0.9178263750828363,
        "241": 0.8887408394403731,
        "242": 0.9231283422459893,
        "243": 0.9685408299866132,
        "244": 0.9009370816599732,
        "245": 0.8181208053691276,
        "246": 0.9040163376446563,
        "247": 0.9993174061433447,
        "248": 0.8651196519216824
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
def target_sequence_Q9NU22():
    return ".........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................FQFFAT-----.........--------------rrllscgg...................."

@pytest.fixture
def target_sequence_continuous_Q9NU22(target_sequence_Q9NU22):
    target_sequence_to_process = target_sequence_Q9NU22
    target_sequence_continuous_Q9NU22 = ''.join([char.upper() for char in target_sequence_to_process if char.isalpha()])
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
def transfer_dict_populated_disulfid_Q9NU22():
    return {
    "DOMAIN": {
        "sequence_id": {
            "sp|Q9NU22|MDN1_HUMAN": {
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
                        "matches": set(["333", "373"]),
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
                        "373": "72"
                        },
                    "aln_to_target": {
                        "18": "333",
                        "72": "373"
                        }
                    }
            }
        }
    }
}

@pytest.fixture
def transfer_dict_populated_disulfid_list_Q9NU22(transfer_dict_populated_disulfid_Q9NU22):
    """Same dictionary but with lists instead of sets"""
    return convert_sets_to_lists(copy.deepcopy(transfer_dict_populated_disulfid_Q9NU22))

@pytest.fixture
def transfer_dict_populated_disulfid_post_processed_Q9NU22(transfer_dict_populated_disulfid_Q9NU22):
    """Creates an extended version of transfer_dict_populated_disulfid_Q9NU22 with GO and conservation data."""
    extended_dict = copy.deepcopy(transfer_dict_populated_disulfid_Q9NU22)
    extended_dict["PF07728"] = extended_dict.pop("DOMAIN")
    extended_dict["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["conservations"]["positions"].update({
        "137": {
            "conservation": 0.956244087038789,
            "hit": True
        },
        "179": {
            "conservation": 0.8452762923351159,
            "hit": True
        }
    })
    extended_dict["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["conservations"]["indices"]["matches"].update(["137", "179"])
    go_data = {
        "GO": {
            "MCRB_ECOLI": {
                "terms": "GO:006305",
                "jaccard_index": "0.45"
            }
        }
    }
    for pos in ["333", "373"]:
        anno_key = next(iter(extended_dict["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["annotations"]["positions"][pos].keys()))
        extended_dict["PF07728"]["sequence_id"]["sp|Q9NU22|MDN1_HUMAN"]["annotations"]["positions"][pos][anno_key].update(go_data)
    return extended_dict

@pytest.fixture
def transfer_dict_populated_disulfid_post_processed_list_Q9NU22(transfer_dict_populated_disulfid_post_processed_Q9NU22):
    """Post-processed version with lists"""
    return convert_sets_to_lists(copy.deepcopy(transfer_dict_populated_disulfid_post_processed_Q9NU22))

# @pytest.fixture
# def transfer_dict_success_disulfid_Q9NU22():
#     return {
#     "match": {
#         "333": {
#             "DISULFID | Intrachain (with C-246); in linked form": {
#                 "essentials": {
#                     "type": "DISULFID",
#                     "description": "Intrachain (with C-246); in linked form",
#                     "count": 1,
#                     "annot_amino_acid": "C",
#                     "target_amino_acid": "C"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:12345678": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "paired_position": {
#                     "373": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "annot_position": {
#                         "205": {
#                             "rep_primary_accession": "P15005",
#                             "rep_mnemo_name": "MCRB_ECOLI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         },
#         "373": {
#             "DISULFID | Intrachain (with C-205); in linked form": {
#                 "essentials": {
#                     "type": "DISULFID",
#                     "description": "Intrachain (with C-205); in linked form",
#                     "count": 1,
#                     "annot_amino_acid": "C",
#                     "target_amino_acid": "R"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:12345678": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "paired_position": {
#                     "333": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "annot_position": {
#                         "246": {
#                             "rep_primary_accession": "P15005",
#                             "rep_mnemo_name": "MCRB_ECOLI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         }
#     },
#     "miss": {
#         "681": {
#             "DISULFID | Intrachain (with C-246); in linked form": {
#                 "essentials": {
#                     "type": "DISULFID",
#                     "description": "Intrachain (with C-246); in linked form",
#                     "count": 1,
#                     "annot_amino_acid": "C",
#                     "target_amino_acid": "T"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:12345678": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "paired_position": {
#                     "717": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "annot_position": {
#                         "205": {
#                             "rep_primary_accession": "P15005",
#                             "rep_mnemo_name": "MCRB_ECOLI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         },
#         "717": {
#             "DISULFID | Intrachain (with C-205); in linked form": {
#                 "essentials": {
#                     "type": "DISULFID",
#                     "description": "Intrachain (with C-205); in linked form",
#                     "count": 1,
#                     "annot_amino_acid": "C",
#                     "target_amino_acid": "K"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:12345678": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "paired_position": {
#                     "681": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "annot_position": {
#                         "246": {
#                             "rep_primary_accession": "P15005",
#                             "rep_mnemo_name": "MCRB_ECOLI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         }
#     }
# }

# @pytest.fixture
# def transfer_dict_success_binding_Q9NU22():
#     return {
#     "match": {
#         "329": {
#             "BINDING | Interacts with GTP": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with GTP",
#                     "count": 1,
#                     "annot_amino_acid": "G",
#                     "target_amino_acid": "G"
#                 },
#                 "evidence": {
#                     "ECO:0000255": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_id": {
#                         "ChEBI:CHEBI:37565": {
#                             "rep_primary_accession": "P15005",
#                             "rep_mnemo_name": "MCRB_ECOLI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         },
#         "330": {
#             "BINDING | Interacts with GTP": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with GTP",
#                     "count": 1,
#                     "annot_amino_acid": "P",
#                     "target_amino_acid": "P"
#                 },
#                 "evidence": {
#                     "ECO:0000255": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_id": {
#                         "ChEBI:CHEBI:37565": {
#                             "rep_primary_accession": "P15005",
#                             "rep_mnemo_name": "MCRB_ECOLI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         },
#         "677": {
#             "BINDING | Interacts with GTP": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with GTP",
#                     "count": 1,
#                     "annot_amino_acid": "G",
#                     "target_amino_acid": "G"
#                 },
#                 "evidence": {
#                     "ECO:0000255": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_id": {
#                         "ChEBI:CHEBI:37565": {
#                             "rep_primary_accession": "P15005",
#                             "rep_mnemo_name": "MCRB_ECOLI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         }
#     },
#     "miss": {
#         "678": {
#             "BINDING | Interacts with GTP": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with GTP",
#                     "count": 1,
#                     "annot_amino_acid": "P",
#                     "target_amino_acid": "E"
#                 },
#                 "evidence": {
#                     "ECO:0000255": {
#                         "rep_primary_accession": "P15005",
#                         "rep_mnemo_name": "MCRB_ECOLI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_id": {
#                         "ChEBI:CHEBI:37565": {
#                             "rep_primary_accession": "P15005",
#                             "rep_mnemo_name": "MCRB_ECOLI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         }
#     }
# }

# @pytest.fixture
# def transfer_dict_success_all_types_H0YB80():
#     return {
#     "match": {
#         "12": {
#             "BINDING | Interacts with peptide": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with peptide",
#                     "count": 1,
#                     "annot_amino_acid": "R",
#                     "target_amino_acid": "R"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:21853016": {
#                         "rep_primary_accession": "Q5CUW0",
#                         "rep_mnemo_name": "Q5CUW0_CRYPI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_ccd_id": {
#                         "peptide": {
#                             "rep_primary_accession": "Q5CUW0",
#                             "rep_mnemo_name": "Q5CUW0_CRYPI",
#                             "count": 1
#                         }
#                     }
#                 }
#             },
#             "BINDING | Interacts with O-phospho-L-serine": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with O-phospho-L-serine",
#                     "count": 1,
#                     "annot_amino_acid": "R",
#                     "target_amino_acid": "R"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:26551337": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 },
#                 "paired_position": {
#                     "13": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_id": {
#                         "ChEBI:CHEBI:57524": {
#                             "rep_primary_accession": "E2RU97",
#                             "rep_mnemo_name": "1433_GIAIC",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         },
#         "13": {
#             "BINDING | Interacts with peptide": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with peptide",
#                     "count": 1,
#                     "annot_amino_acid": "Y",
#                     "target_amino_acid": "Y"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:21853016": {
#                         "rep_primary_accession": "Q5CUW0",
#                         "rep_mnemo_name": "Q5CUW0_CRYPI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_ccd_id": {
#                         "peptide": {
#                             "rep_primary_accession": "Q5CUW0",
#                             "rep_mnemo_name": "Q5CUW0_CRYPI",
#                             "count": 1
#                         }
#                     }
#                 }
#             },
#             "BINDING | Interacts with O-phospho-L-serine": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with O-phospho-L-serine",
#                     "count": 1,
#                     "annot_amino_acid": "Y",
#                     "target_amino_acid": "Y"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:26551337": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 },
#                 "paired_position": {
#                     "12": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_id": {
#                         "ChEBI:CHEBI:57524": {
#                             "rep_primary_accession": "E2RU97",
#                             "rep_mnemo_name": "1433_GIAIC",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         },
#         "57": {
#             "BINDING | Interacts with peptide": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with peptide",
#                     "count": 2,
#                     "annot_amino_acid": "L",
#                     "target_amino_acid": "L"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:21853016": {
#                         "rep_primary_accession": "Q5CUW0",
#                         "rep_mnemo_name": "Q5CUW0_CRYPI",
#                         "count": 1
#                     },
#                     "ECO:0000269|PubMed:26551337": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_ccd_id": {
#                         "peptide": {
#                             "rep_primary_accession": "Q5CUW0",
#                             "rep_mnemo_name": "Q5CUW0_CRYPI",
#                             "count": 2
#                         }
#                     }
#                 }
#             }
#         },
#         "58": {
#             "BINDING | Interacts with peptide": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with peptide",
#                     "count": 2,
#                     "annot_amino_acid": "N",
#                     "target_amino_acid": "N"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:21853016": {
#                         "rep_primary_accession": "Q5CUW0",
#                         "rep_mnemo_name": "Q5CUW0_CRYPI",
#                         "count": 1
#                     },
#                     "ECO:0000269|PubMed:26551337": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_ccd_id": {
#                         "peptide": {
#                             "rep_primary_accession": "Q5CUW0",
#                             "rep_mnemo_name": "Q5CUW0_CRYPI",
#                             "count": 2
#                         }
#                     }
#                 }
#             }
#         },
#         "65": {
#             "BINDING | Interacts with peptide": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with peptide",
#                     "count": 1,
#                     "annot_amino_acid": "E",
#                     "target_amino_acid": "E"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:21853016": {
#                         "rep_primary_accession": "Q5CUW0",
#                         "rep_mnemo_name": "Q5CUW0_CRYPI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_ccd_id": {
#                         "peptide": {
#                             "rep_primary_accession": "Q5CUW0",
#                             "rep_mnemo_name": "Q5CUW0_CRYPI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         },
#         "109": {
#             "BINDING | Interacts with peptide": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with peptide",
#                     "count": 2,
#                     "annot_amino_acid": "N",
#                     "target_amino_acid": "N"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:21853016": {
#                         "rep_primary_accession": "Q5CUW0",
#                         "rep_mnemo_name": "Q5CUW0_CRYPI",
#                         "count": 1
#                     },
#                     "ECO:0000269|PubMed:26551337": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_ccd_id": {
#                         "peptide": {
#                             "rep_primary_accession": "Q5CUW0",
#                             "rep_mnemo_name": "Q5CUW0_CRYPI",
#                             "count": 2
#                         }
#                     }
#                 }
#             }
#         },
#         "113": {
#             "BINDING | Interacts with peptide": {
#                 "essentials": {
#                     "type": "BINDING",
#                     "description": "Interacts with peptide",
#                     "count": 1,
#                     "annot_amino_acid": "W",
#                     "target_amino_acid": "W"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:21853016": {
#                         "rep_primary_accession": "Q5CUW0",
#                         "rep_mnemo_name": "Q5CUW0_CRYPI",
#                         "count": 1
#                     }
#                 },
#                 "additional_keys": {
#                     "ligand_ccd_id": {
#                         "peptide": {
#                             "rep_primary_accession": "Q5CUW0",
#                             "rep_mnemo_name": "Q5CUW0_CRYPI",
#                             "count": 1
#                         }
#                     }
#                 }
#             }
#         },
#         "61": {
#             "MUTAGEN | V->D: Loss of binding to difopein.": {
#                 "essentials": {
#                     "type": "MUTAGEN",
#                     "description": "V->D: Loss of binding to difopein.",
#                     "count": 1,
#                     "annot_amino_acid": "V",
#                     "target_amino_acid": "V"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:19733174": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 }
#             }
#         }
#     },
#     "miss": {
#         "78": {
#             "MUTAGEN | R->K: Increased oligomerization.": {
#                 "essentials": {
#                     "type": "MUTAGEN",
#                     "description": "R->K: Increased oligomerization.",
#                     "count": 1,
#                     "annot_amino_acid": "R",
#                     "target_amino_acid": "K"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:24658679": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 }
#             }
#         },
#         "86": {
#             "MUTAGEN | T->A: Slightly decreased oligomerization.": {
#                 "essentials": {
#                     "type": "MUTAGEN",
#                     "description": "T->A: Slightly decreased oligomerization.",
#                     "count": 1,
#                     "annot_amino_acid": "T",
#                     "target_amino_acid": "A"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:24658679": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 }
#             }
#         },
#         "92": {
#             "MOD_RES | Phosphothreonine": {
#                 "essentials": {
#                     "type": "MOD_RES",
#                     "description": "Phosphothreonine",
#                     "count": 1,
#                     "annot_amino_acid": "T",
#                     "target_amino_acid": "S"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:16368691, ECO:0000269|PubMed:19733174, ECO:0000269|PubMed:24147113, ECO:0000269|PubMed:24658679": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 }
#             },
#             "MUTAGEN | T->A: Loss of phosphorylation by a protein kinase. No effect on subcellular localization. Dramatic decrease in the number of encysting parasites and cysts, but a large increase in the number of trophozoites. In encysting cells of 12 hours, significantly slower cyst conversion rate compared to the wild-type. No effect on binding to difopein. Decreased binding to a number of synthetic phosphopeptides.": {
#                 "essentials": {
#                     "type": "MUTAGEN",
#                     "description": "T->A: Loss of phosphorylation by a protein kinase. No effect on subcellular localization. Dramatic decrease in the number of encysting parasites and cysts, but a large increase in the number of trophozoites. In encysting cells of 12 hours, significantly slower cyst conversion rate compared to the wild-type. No effect on binding to difopein. Decreased binding to a number of synthetic phosphopeptides.",
#                     "count": 1,
#                     "annot_amino_acid": "T",
#                     "target_amino_acid": "S"
#                 },
#                 "evidence": {
#                     "ECO:0000269|PubMed:16368691, ECO:0000269|PubMed:19733174": {
#                         "rep_primary_accession": "E2RU97",
#                         "rep_mnemo_name": "1433_GIAIC",
#                         "count": 1
#                     }
#                 }
#             }
#         }
#     }
# }

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
def annotation_dict_201_Q9NU22(annotations_content_binding_fixture_Q9NU22_PF07728, get_annotation_dict):
    return get_annotation_dict(annotations_content_binding_fixture_Q9NU22_PF07728, "201")

@pytest.fixture
def annotation_dict_202_Q9NU22(annotations_content_binding_fixture_Q9NU22_PF07728, get_annotation_dict):
    return get_annotation_dict(annotations_content_binding_fixture_Q9NU22_PF07728, "202")

@pytest.fixture
def mock_make_anno_total_disulfid_return_205_P15005(annotation_dict_205_Q9NU22):
    return {
        'annotation': annotation_dict_205_Q9NU22,
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-246); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-246); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678',
            'target_position': '333',
            'annot_position': '205',
            'index_position': '18',
            'annot_amino_acid': 'C',
            'target_amino_acid': 'C',
            'paired_target_position': '373'
        },
        'paired_annot_pos_str': '246'
    }

@pytest.fixture
def mock_make_anno_total_disulfid_return_246_P15005(annotation_dict_205_Q9NU22):
   return {
        'annotation': annotation_dict_205_Q9NU22,
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-205); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-205); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678',
            'target_position': '373',
            'annot_position': '246',
            'index_position': '72',
            'annot_amino_acid': 'C',
            'target_amino_acid': 'C',
            'paired_target_position': '333'
        },
        'paired_annot_pos_str': '205'
    }

@pytest.fixture
def mock_make_anno_total_disulfid_return_205_Q9NU22():
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
            'index_position': '18',
            'annot_amino_acid': 'C',
            'target_amino_acid': 'C',
            'paired_target_position': '373'
        },
        'paired_annot_pos_str': '246'
    }

@pytest.fixture
def mock_make_anno_total_disulfid_return_246_Q9NU22():
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
            'index_position': '72',
            'annot_amino_acid': 'C',
            'target_amino_acid': 'C',
            'paired_target_position': '333'
        },
        'paired_annot_position': '205'
    }


@pytest.fixture
def mock_map_filter_disulfid_return(annotation_dict_246_Q9NU22):
    return (True, {
        'annotation': annotation_dict_246_Q9NU22,
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Intrachain (with C-205); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Intrachain (with C-205); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678',
            'target_position': '373',
            'annot_position': '246',
            'index_position': '72',
            'annot_amino_acid': 'C',
            'target_amino_acid': 'C',
            'paired_target_position': '333'
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

@pytest.fixture
def minimal_hmmalign_lines_fixture_Q9NU22():
    return """# STOCKHOLM 1.0
Q7US48_RHOBA/138-284                          .GLMMVGEPGTAKSMLGELLAVAIS....GTGSLAVQGTAGTTEDQIKYGWNYAmll......drgPVHEALVPSPVMTAMR..........DGKVVRFEEITRCL.PEVQDALISILSER.R...MMIPEMDQAGDDNANSvFASP..................gFSIIATANLRD.......RGVSEMSAALKRRF.......................
MCRB_ECOLI/196-350                                      .........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYRCN............GVGFRRKDGIFYNFCQqak...........eqpeKKYIFIIDEINRANlS........K....VFG....E......VMM.LMEHD.K...RGENWSVPLTYSENDE.eRFYVp..........................enVYIIGLMNTAD.........RSLAVVDYALRRRF............................
sp|Q9NU22|MDN1_HUMANtarget//325-451                     .........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYRCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................-FQFFAT-----.........--------------rrllscgg....................
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
@pytest.fixture
def iprscan_tsv_mock(tmp_path):
    """Creates a mock InterProScan TSV file without header."""
    iprscan_file = tmp_path / "target_name" / "iprscan.tsv"
    iprscan_file.parent.mkdir(parents=True, exist_ok=True)
    iprscan_file.write_text(
        "sp|Q9NU22|MDN1_HUMAN\t7b57368205d60b6e7b538d2181ad7f2b\t5596\tPfam\tPF07728\t"
        "AAA domain (dynein-related subfamily)\t325\t448\t2.7E-13\tT\t07-01-2025\tIPR011704\t"
        "ATPase, dynein-related, AAA domain\tGO:0005524(InterPro)|GO:0016887(InterPro)\t-\n"
    )
    return str(tmp_path)

@pytest.fixture
def iprscan_tsv_mock_no_go(tmp_path):
    """Creates a mock InterProScan TSV file with '-' for GO terms."""
    iprscan_file = tmp_path / "target_name" / "iprscan.tsv"
    iprscan_file.parent.mkdir(parents=True, exist_ok=True)
    iprscan_file.write_text(
        "sp|Q9NU22|MDN1_HUMAN\t7b57368205d60b6e7b538d2181ad7f2b\t5596\tSMART\tSM00327\t"
        "VWA_4\t5382\t5563\t6.1E-10\tT\t07-01-2025\tIPR002035\t"
        "von Willebrand factor, type A\t-\t-\n"
    )
    return str(tmp_path)

@pytest.fixture
def target_id_plus_seq_mock(minimal_hmmalign_lines_fixture_Q9NU22):
    """Returns target ID and sequence as tuple."""
    lines = minimal_hmmalign_lines_fixture_Q9NU22.splitlines()
    try:
        line = next(line for line in lines if line.startswith("sp|Q9NU22|MDN1_HUMANtarget//325-451"))
    except StopIteration:
        raise ValueError("Target ID 'sp|Q9NU22|MDN1_HUMANtarget//325-451' not found in hmmalign_lines.")
    id_seq = line.split()[:2]  # Returns [id, seq]
    return id_seq

@pytest.fixture
def conservation_id_plus_seq_mock(minimal_hmmalign_lines_fixture_Q9NU22):
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
# paired_position_res_hit/hit = Define in the test ???
# paired_anno_id = Define in the test ???
# paired_anno_total = Define in the test ???

### Added in clenup_improve_transfer_dict and helper functions
pfam_id_mock_Q9NU22 = "PF07728"

###T parse_arguments
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

###T get_pfam_id_from_hmmalign_result
def test_get_pfam_id_from_hmmalign_result():
    assert get_pfam_id_from_hmmalign_result(hmmalign_result_mock) == "PF07728"
    hmmalign_result_mock_2 = "/home/user/results/human/PF00001_hmmalign.sth"
    assert get_pfam_id_from_hmmalign_result(hmmalign_result_mock_2) == "PF00001"

###T get_annotation_filepath
def test_get_annotation_filepath():
    assert get_annotation_filepath(resource_dir_mock, 'PF07728') == ("/home/user/resources/PF07728/annotations.json", "/home/user/resources/PF07728/conservations.json")

###T read_files
def test_read_files(annotations_content_binding_fixture_Q9NU22_PF07728):
    annotations_content = json.dumps(annotations_content_binding_fixture_Q9NU22_PF07728)

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
        assert hmmalign_lines == hmmalign_result_content_mock.splitlines()
        assert annotations == annotations_content_binding_fixture_Q9NU22_PF07728

###T iterate_aligned_sequences

@pytest.fixture
def basic_seq_pair():
    return "ABC", "XYZ", 1, 10, 3, 12

@pytest.fixture
def gapped_seq_pair():
    return "A-BC", "X-YZ", 1, 10, 3, 12

def test_iterate_aligned_sequences_basic(basic_seq_pair):
    seq_a, seq_b, start_a, start_b, end_a, end_b = basic_seq_pair
    result = list(iterate_aligned_sequences(seq_a, seq_b, start_a, start_b, end_a, end_b))

    expected = [
        (0, 1, 10, 'A', 'X'),
        (1, 2, 11, 'B', 'Y'),
        (2, 3, 12, 'C', 'Z')
    ]

    assert result == expected

def test_iterate_aligned_sequences_with_gaps(gapped_seq_pair):
    seq_a, seq_b, start_a, start_b, end_a, end_b = gapped_seq_pair
    result = list(iterate_aligned_sequences(seq_a, seq_b, start_a, start_b, end_a, end_b))

    expected = [
        (0, 1, 10, 'A', 'X'),
        (1, 1, 10, '-', '-'),
        (2, 2, 11, 'B', 'Y'),
        (3, 3, 12, 'C', 'Z')
    ]

    assert result == expected

def test_iterate_aligned_sequences_early_stop():
    result = list(iterate_aligned_sequences(
        "ABCDEF", "XYZUVW",
        start_a=1, start_b=10,
        end_a=3, end_b=12  # Should stop after 'C'/'Z'
    ))

    expected = [
        (0, 1, 10, 'A', 'X'),
        (1, 2, 11, 'B', 'Y'),
        (2, 3, 12, 'C', 'Z')
    ]

    assert result == expected

def test_iterate_aligned_sequences_empty():
    result = list(iterate_aligned_sequences("", "", 0, 0, 1, 1))
    assert result == []


###T find_and_map_annots
def test_find_and_map_annots_basic_case(annotations_content_binding_fixture_Q9NU22_PF07728, logger, target_sequence_Q9NU22, annot_sequence_Q9NU22):
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
            annotations=annotations_content_binding_fixture_Q9NU22_PF07728,
            good_eco_codes=[]
        )

        mock_map_and_filter.assert_called_with(
            logger=logger,
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

def test_find_and_map_annots_no_target_sequence(annotations_content_binding_fixture_Q9NU22_PF07728, logger):
    """Test case with no matching annotations"""
    mock_lines = [
        "# STOCKHOLM 1.0\n",
        "VWA8_HUMAN/105-261    .........DVFLIGPPGPLRRSIAM.QYLELT..............KREVEYIALSR.DT..TETDLKQRREIR\n"
    ]

    transfer_dict = find_and_map_annots(
        logger=logger,
        hmmalign_lines=mock_lines,
        annotations=annotations_content_binding_fixture_Q9NU22_PF07728,
        good_eco_codes=[]
    )

    assert not transfer_dict

    # Verify the logger captured the correct error message
    logger.error.assert_called_once_with(
        "---> ERROR --- FIND_AND_MAP --- No target sequences found in hmmalign lines!"
    )

###T read_conservations_and_annotations

def test_read_conservations_and_annotations_success(
    tmp_path, mock_json_filepaths, conservations_content_Q9NU22_PF07728, annotations_content_disulfid_fixture_Q9NU22_PF07728):
    # Create temporary files
    cons_file, annot_file = mock_json_filepaths

    with open(cons_file, 'w') as f:
        json.dump(conservations_content_Q9NU22_PF07728, f)
    with open(annot_file, 'w') as f:
        json.dump(annotations_content_disulfid_fixture_Q9NU22_PF07728, f)

    # Test function
    conservations, annotations = read_conservations_and_annotations(str(cons_file), str(annot_file))

    assert conservations == conservations_content_Q9NU22_PF07728
    assert annotations == annotations_content_disulfid_fixture_Q9NU22_PF07728

def test_read_conservations_and_annotations_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_conservations_and_annotations(
            str(tmp_path / "nonexistent.json"),
            str(tmp_path / "annotations.json")
        )

def test_read_conservations_and_annotations_invalid_json(tmp_path, mock_json_filepaths):
    # Create file with invalid JSON
    cons_file, annot_file = mock_json_filepaths

    with open(cons_file, 'w') as f:
        f.write("invalid json")
    with open(annot_file, 'w') as f:
        json.dump({}, f)

    with pytest.raises(json.JSONDecodeError):
        read_conservations_and_annotations(str(cons_file), str(annot_file))

def test_read_conservations_and_annotations_empty_files(mock_json_filepaths):
    """Test handling of empty but valid JSON files"""
    cons_file, annot_file = mock_json_filepaths

    with open(cons_file, 'w', encoding='utf-8') as f:
        json.dump({}, f)
    with open(annot_file, 'w', encoding='utf-8') as f:
        json.dump({}, f)

    conservations, annotations = read_conservations_and_annotations(str(cons_file), str(annot_file))
    assert conservations == {}
    assert annotations == {}

###T parse_go_annotations



###T map_and_filter_annot_pos
def test_map_and_filter_with_paired_positions(
    logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_disulfid_pair, transfer_dict,
    get_annotation_dict, annotations_content_disulfid_fixture_Q9NU22_PF07728
):
    """Unit test for map_and_filter_annot_pos with paired positions."""
    annotation_dict = get_annotation_dict(annotations_content_disulfid_fixture_Q9NU22_PF07728, "246")

    with patch('transfer_annotations.validate_paired_annotations') as mock_validate_paired_annotations:
        # Configure the mock to return expected values
        mock_validate_paired_annotations.return_value = (True, {'dummy': 'result'})

        # Call the function under test
        result = map_and_filter_annot_pos(
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
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
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
        )

        # Verify the result
        assert result == (True, {'dummy': 'result'})

def test_map_and_filter_without_paired_positions(
    logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, transfer_dict
):
    """Unit test for map_and_filter_annot_pos without paired positions."""

    with patch('transfer_annotations.validate_annotations') as mock_validate_annotations:
        # Call the function under test
        result = map_and_filter_annot_pos(
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
            entry_annotations=entry_annotations_binding_only,
            transfer_dict=transfer_dict,
            processed_annotations=set(),
            # No paired position parameters
        )

        # Assert that validate_annotations was called with expected arguments
        mock_validate_annotations.assert_called_once_with(
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
    logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, transfer_dict
):
    """Test that process_annotation is called with correct parameters"""

    with patch('transfer_annotations.process_annotation') as mock_process:
        validate_annotations(
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
                annotation_key=('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '201', '329', 'BINDING')
            ),
            call(
                res_hit=True,
                logger=logger,
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
                counter_target_pos_str='330',
                counter_annot_pos_str="202",
                index=15,
                target_amino="P",
                processed_annotations=set(),
                annotation_key=('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '202', '330', 'BINDING')
            )
        ]

        mock_process.assert_has_calls(expected_calls, any_order=False)

def test_validate_annotations_multiple_annotations(
    logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
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
    logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
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
    logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, transfer_dict
):
    """Test case 2 - No annotations in offset range (201 and 202 don't fit in 500-600)"""
    with patch('transfer_annotations.process_annotation') as mock_process:
        validate_annotations(
            logger=logger,
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
    logger, good_eco_codes_all, target_sequence_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, transfer_dict
):
    """Test that error logs are written when process_annotation fails."""
    mock_error = Exception("Test error")

    with patch('transfer_annotations.process_annotation', side_effect=mock_error):
        try:
            validate_annotations(
                logger=logger,
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

        # Verify error log was written
        logger.error.assert_any_call("---> ERROR --- VAL_ANNOTS --- Error in validate_annotations: Test error")


###T process_annotation
def test_process_annotation_basic_single(
    logger, good_eco_codes_all, transfer_dict,
    target_sequence_Q9NU22, target_sequence_continuous_Q9NU22, annot_sequence_Q9NU22,
    entry_annotations_binding_only, get_annotation_dict,
    annotations_content_binding_fixture_Q9NU22_PF07728
):
    """Test processing of a single BINDING annotation"""
    annotation_dict = get_annotation_dict(annotations_content_binding_fixture_Q9NU22_PF07728, "201")
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
            'annot_position': '201',
            'index_position': '14',
            'annot_amino_acid': 'G',
            'target_amino_acid': 'G',
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
            annotation_key=('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '201', '329', 'BINDING')
        )

        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict,
            counter_target_pos_str='329',
            counter_annot_pos_str="201",
            index=14,
            target_amino="G",
            logger=logger,
            entry_annotations=entry_annotations_binding_only
        )

        mock_add_dict.assert_called_once_with(
            hit=True,
            logger=logger,
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

        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '201', '329', 'BINDING') in processed_annotations

def test_process_annotation_paired_failure(
    logger, good_eco_codes_all, transfer_dict,
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

    with patch('transfer_annotations.make_anno_total_dict', return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch('transfer_annotations.map_and_filter_annot_pos', return_value=mock_map_filter_disulfid_fail_return) as mock_map_filter, \
         patch('transfer_annotations.add_to_transfer_dict') as mock_add_transfer:
        process_annotation(
            res_hit=True,
            logger=logger,
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
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID")
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair
        )

        # Verify map_and_filter_annot_pos was called correctly
        mock_map_filter.assert_called_once_with(
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
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            paired_annot_pos_str="246",
            caller_target_pos_str="333",
            paired_annotation_dict=annotation_dict_246_Q9NU22
        )

        # Verify remove_failed_annotations was called correctly
        mock_add_transfer.assert_called_once_with(
            hit=True,
            logger=logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_mod_to_fail_73_continuous,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id='DISULFID | Intrachain (with C-246); in linked form',
            anno_total=mock_make_anno_total_disulfid_return_205_P15005['anno_total'],
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            paired_position_res_hit=False,
            paired_anno_id='DISULFID | Intrachain (with C-205); in linked form',
            paired_anno_total=mock_map_filter_disulfid_return[1]['anno_total']
        )

        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '205', '333', 'DISULFID') in processed_annotations
        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '246', '373', 'DISULFID') in processed_annotations


def test_process_annotation_paired_success(
    logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair, annotation_dict_205_Q9NU22,
    annotation_dict_246_Q9NU22, mock_make_anno_total_disulfid_return_205_P15005,
    mock_map_filter_disulfid_return,
    target_sequence_Q9NU22, target_sequence_continuous_Q9NU22, annot_sequence_Q9NU22
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
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID")
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair
        )

        # Verify map_and_filter_annot_pos was called correctly
        mock_map_filter.assert_called_once_with(
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
            entry_annotations=entry_annotations_disulfid_pair,
            transfer_dict=transfer_dict,
            paired_annot_pos_str="246",
            caller_target_pos_str="333",
            paired_annotation_dict=annotation_dict_246_Q9NU22
        )

        # Verify remove_failed_annotations was called correctly
        mock_add_transfer.assert_called_once_with(
            hit=True,
            logger=logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_continuous_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id='DISULFID | Intrachain (with C-246); in linked form',
            anno_total=mock_make_anno_total_disulfid_return_205_P15005['anno_total'],
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            paired_position_res_hit=True,
            paired_anno_id='DISULFID | Intrachain (with C-205); in linked form',
            paired_anno_total=mock_map_filter_disulfid_return[1]['anno_total']
        )

        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '205', '333', 'DISULFID') in processed_annotations
        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '246', '373', 'DISULFID') in processed_annotations


def test_process_annotation_paired_no_paired_dict(
    logger, good_eco_codes_all, transfer_dict,
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

    with patch('transfer_annotations.make_anno_total_dict', return_value=mock_make_anno_total_disulfid_return_205_P15005) as mock_make_total, \
         patch('transfer_annotations.add_to_transfer_dict') as mock_add_transfer:
        process_annotation(
            res_hit=True,
            logger=logger,
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
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            processed_annotations=processed_annotations,
            annotation_key=("MCRB_ECOLI", "sp|Q9NU22|MDN1_HUMAN", "205", "333", "DISULFID")
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict_205_Q9NU22,
            counter_target_pos_str='333',
            counter_annot_pos_str="205",
            index=18,
            target_amino="C",
            logger=logger,
            entry_annotations=modified_entry_annotations
        )

        # Verify remove_failed_annotations was called correctly
        mock_add_transfer.assert_called_once_with(
            hit=True,
            logger=logger,
            transfer_dict=transfer_dict,
            target_name=target_name_mock_Q9NU22,
            target_sequence_continuous=target_sequence_continuous_Q9NU22,
            target_hit_start=target_hit_start_mock_Q9NU22,
            target_hit_end=target_hit_end_mock_Q9NU22,
            anno_id='DISULFID | Intrachain (with C-246); in linked form',
            anno_total=mock_make_anno_total_disulfid_return_205_P15005['anno_total'],
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            paired_position_res_hit=False,
            paired_anno_id=None,
            paired_anno_total=None
        )

        assert ('MCRB_ECOLI', 'sp|Q9NU22|MDN1_HUMAN', '205', '333', 'DISULFID') in processed_annotations



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
        'annotation': {
            **annotation_dict_205_Q9NU22,
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
            'index_position': '18',
            'annot_amino_acid': 'C',
            'target_amino_acid': 'C',
        },
        'paired_annot_pos_str': '246'
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
        counter_target_pos_str='333',
        counter_annot_pos_str="205",
        index=18,
        target_amino="C",
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair
    )

    expected_single_pass = {
        'annotation': {
            **annotation_dict_205_Q9NU22,
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
            'index_position': '18',
            'annot_amino_acid': 'C',
            'target_amino_acid': 'C',
        },
        'paired_annot_pos_str': '246'
    }

    expected_multiple_pass = copy.deepcopy(expected_single_pass)
    expected_multiple_pass["anno_total"]["evidence"] = "ECO:0000269|PubMed:12345678, ECO:0007744|PDB:1ABC"

    assert result_single_pass == expected_single_pass

    good_eco_codes_all_modded = good_eco_codes_all + ["ECO:0007744"]

    result_multiple_pass = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all_modded,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_dict_205_Q9NU22,
        counter_target_pos_str='333',
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
        counter_target_pos_str='329',
        counter_annot_pos_str="201",
        index=14,
        target_amino="G",
        logger=logger,
        entry_annotations=entry_annotations_binding_only
    )

    assert result['anno_total'] is None

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
        counter_target_pos_str='329',
        counter_annot_pos_str="201",
        index=14,
        target_amino="G",
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
            'annot_position': '201',
            'index_position': '14',
            'annot_amino_acid': 'G',
            'target_amino_acid': 'G'
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
    annotation_dict_246_Q9NU22
):
    """Test with caller_target_pos provided (paired annotation case)"""
    result = make_anno_total_dict(
        good_eco_codes=good_eco_codes_all,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_dict_246_Q9NU22,
        counter_target_pos_str='373',
        counter_annot_pos_str="246",
        caller_target_pos_str='333',
        index=18,
        target_amino="C",
        logger=logger,
        entry_annotations=entry_annotations_disulfid_pair
    )

    expected = {
        'annotation': {
            **annotation_dict_246_Q9NU22,
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
            'index_position': '18',
            'annot_amino_acid': 'C',
            'target_amino_acid': 'C',
            'paired_target_position': '333',
        },
        'paired_annot_pos_str': '205'
    }

    assert result['annotation'] is annotation_dict_246_Q9NU22
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
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        annotation_dict=annotation_without_evidence,
        counter_target_pos_str='100',
        counter_annot_pos_str="201",
        index=14,
        target_amino='G',
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
    target_sequence_Q9NU22,
    annot_sequence_Q9NU22,
    entry_annotations_disulfid_pair,
    annotation_dict_246_Q9NU22
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
            "index_position": "72",
            "annot_amino_acid": "C",
            "target_amino_acid": "C"
        },
        "paired_annot_pos_str": "205",
    }

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
            paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22, # '246'
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22, # '333'
        )

        # Verify the make_anno_total_dict call
        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name="MCRB_ECOLI",
            annotation_dict=annotation_dict_246_Q9NU22,
            counter_target_pos_str='373', # 246 = 373 in target sequence numbering/col
            counter_annot_pos_str=paired_annot_pos_str_mock_Q9NU22, # '246'
            index=72,
            target_amino="C",
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair,
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22,
        )

        # Assert that the result matches the mock data
        assert result_tuple == (True, mock_paired_result_dict)

#RESOLVER
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
            paired_annot_pos_str=paired_annot_pos_str_mock_Q9NU22, # '246'
            caller_target_pos_str=caller_target_pos_str_mock_Q9NU22, # '333'
        )

        mock_make_total.assert_called_once_with(
            good_eco_codes=good_eco_codes_all,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            annotation_dict=annotation_dict_246_Q9NU22,
            counter_target_pos_str='373', # 246 = 373 in target sequence numbering/col
            counter_annot_pos_str=paired_annot_pos_str_mock_Q9NU22, # '246'
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
    )

    assert result == (False, {})

###T add_to_transfer_dict
def test_add_to_transfer_dict_disulfid_single_first_addition(
    logger, transfer_dict,
    target_sequence_continuous_Q9NU22,
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
        target_name=target_name_mock_Q9NU22,
        target_sequence_continuous=target_sequence_continuous_Q9NU22,
        target_hit_start=target_hit_start_mock_Q9NU22,
        target_hit_end=target_hit_end_mock_Q9NU22,
        anno_id=anno_id,
        anno_total=anno_total,
        entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
        entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
    )
    ### DELETE
    # transfo_dicto = convert_sets_to_lists(transfer_dict)
    # print(json.dumps(transfo_dicto, indent=4))
    ### DELETE
    assert target_name_mock_Q9NU22 in transfer_dict['DOMAIN']['sequence_id']
    assert 'positions' in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']
    assert '333' in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']
    assert anno_total['target_position'] in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']
    assert anno_id in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions'][anno_total['target_position']]
    assert transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions'][anno_total['target_position']][anno_id]['essentials']['type'] == 'DISULFID'

def test_add_to_transfer_dict_paired_disulfide(
    logger, transfer_dict,
    target_sequence_continuous_Q9NU22,
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

    # Assert main position
    assert '333' in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']
    assert anno_id in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']['333']

    # Assert paired position
    assert '373' in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']
    assert paired_anno_id in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']['373']

    # Assert cross-references
    assert '373' in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']['333'][anno_id]['paired_position']
    assert '333' in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']['373'][paired_anno_id]['paired_position']


def test_add_to_transfer_dict_paired_repeated_anno_id(
    logger, transfer_dict_populated_disulfid_Q9NU22,
    target_sequence_continuous_Q9NU22,
    mock_make_anno_total_disulfid_return_205_Q9NU22,
    mock_make_anno_total_disulfid_return_246_Q9NU22,
):
    """Test paired disulfide annotation processing"""
    anno_id = mock_make_anno_total_disulfid_return_205_Q9NU22['anno_id']
    anno_total = mock_make_anno_total_disulfid_return_205_Q9NU22['anno_total']
    paired_anno_id = mock_make_anno_total_disulfid_return_246_Q9NU22['anno_id']
    paired_anno_total = mock_make_anno_total_disulfid_return_246_Q9NU22['anno_total']

    add_to_transfer_dict(
        hit=True,
        logger=logger,
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
    def get_annotation(target_position_str):
        return transfer_dict_populated_disulfid_Q9NU22['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions'][target_position_str]

    # Test data
    positions = {
        '333': {'partner': '373', 'annot_pos': '205'},
        '373': {'partner': '333', 'annot_pos': '246'}
    }

    for target_position_str, data in positions.items():
        anno = get_annotation(target_position_str)['DISULFID | Intrachain (with C-246); in linked form' if target_position_str == '333' else 'DISULFID | Intrachain (with C-205); in linked form']

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
    # print(json.dumps(transfer_dict_populated_disulfid_Q9NU22, indent=4))
    assert transfer_dict_populated_disulfid_Q9NU22['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']['333'][anno_id]['hit'] is True
    assert transfer_dict_populated_disulfid_Q9NU22['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']['373'][paired_anno_id]['hit'] is True

def test_add_to_transfer_dict_hit_miss_pair(
    logger, transfer_dict,
    target_sequence_continuous_Q9NU22,
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

    # Assert hit position
    assert '333' in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']
    assert transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']['333'][anno_id]['hit'] is True

    # Assert missed position
    assert '373' in transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']
    assert transfer_dict['DOMAIN']['sequence_id'][target_name_mock_Q9NU22]['annotations']['positions']['373'][paired_anno_id]['hit'] is False



def test_add_to_transfer_dict_disulfid_single_mocked_helper(
    logger, transfer_dict,
    target_sequence_continuous_Q9NU22,
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
            anno_id=anno_id,
            anno_total=anno_total,
            entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
            entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
            additional_keys={}
        )

def test_add_to_transfer_dict_disulfid_paired_mocked_helper(
    logger, transfer_dict,
    target_sequence_continuous_Q9NU22,
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
                anno_id=paired_anno_id,
                anno_total=paired_anno_total,
                entry_mnemo_name=entry_mnemo_name_mock_Q9NU22,
                entry_primary_accession=entry_primary_accesion_mock_Q9NU22,
                additional_keys={}
            )
        ]

        mock_add_single.assert_has_calls(expected_calls, any_order=False)
        assert mock_add_single.call_count == 2

###T gather_go_terms_for_target
# YAY
def test_gather_go_terms_for_target(iprscan_tsv_mock):
    go_terms = gather_go_terms_for_target("target_name", iprscan_tsv_mock)
    assert go_terms == {"GO:0005524", "GO:0016887"}

def test_gather_go_terms_for_target_no_go(iprscan_tsv_mock_no_go):
    go_terms = gather_go_terms_for_target("target_name", iprscan_tsv_mock_no_go)
    assert go_terms == set()

###T get_alignment_sequences

def test_get_alignment_sequences(minimal_hmmalign_lines_fixture_Q9NU22, target_id_plus_seq_mock, conservation_id_plus_seq_mock):
    hmmalign_lines_list = minimal_hmmalign_lines_fixture_Q9NU22.splitlines()
    target_seq, conservation_seq = get_alignment_sequences(
        hmmalign_lines_list,
        "sp|Q9NU22|MDN1_HUMANtarget//325-451",
        "Q7US48_RHOBA/138-284"
    )

    _, target_seq_expected = target_id_plus_seq_mock
    _, conservation_seq_expected = conservation_id_plus_seq_mock

    assert target_seq == target_seq_expected
    assert conservation_seq == conservation_seq_expected

def test_get_alignment_sequence_no_target(minimal_hmmalign_lines_fixture_Q9NU22):
    hmmalign_lines_list = minimal_hmmalign_lines_fixture_Q9NU22.splitlines()

    # Remove the line with the target ID
    hmmalign_lines_list = [line for line in hmmalign_lines_list if not line.startswith("sp|Q9NU22|MDN1_HUMANtarget//325-451")]

    with pytest.raises(ValueError, match="Target ID 'sp|Q9NU22|MDN1_HUMANtarget//325-451' not found in hmmalign_lines."):
        get_alignment_sequences(hmmalign_lines_list, "sp|Q9NU22|MDN1_HUMANtarget//325-451", "Q7US48_RHOBA/138-284")

def test_get_alignment_sequence_no_conservation(minimal_hmmalign_lines_fixture_Q9NU22):
    hmmalign_lines_list = minimal_hmmalign_lines_fixture_Q9NU22.splitlines()

    # Remove the line with the conservation ID
    hmmalign_lines_list = [line for line in hmmalign_lines_list if not line.startswith("Q7US48_RHOBA/138-284")]

    with pytest.raises(ValueError, match="Conservation ID 'Q7US48_RHOBA/138-284' not found in hmmalign_lines."):
        get_alignment_sequences(hmmalign_lines_list, "sp|Q9NU22|MDN1_HUMANtarget//325-451", "Q7US48_RHOBA/138-284")



###T populate_conservation

# def test_populate_conservation(
#     transfer_dict_populated_disulfid_Q9NU22, pfam_id_mock, target_name_mock_Q9NU22, target_sequence_Q9NU22, conservation_seq_mock, conservation_key_mock, conservations_mock, conservations_positions_mock, target_hit_start_mock_Q9NU22, target_hit_end_mock_Q9NU22, conservation_start_mock, conservation_end_mock
#     logger):
#     target_id, target_seq = 
#     conservation_seq = "GPPGCGKTFVAR"
#     conservations = {"conservation_key": {"1": {"score": 5}, "3": {"score": 10}}}
#     iterate_func = lambda seq_a, seq_b, start_a, start_b, end_a, end_b: [
#         (0, 1, 1, "G", "G"),
#         (1, 2, 2, "P", "P"),
#         (2, 3, 3, "P", "Q"),
#     ]
#     populate_conservation(
#         transfer_dict=transfer_dict,
#         pfam_id="PF00001",
#         target_name="target",
#         target_seq=target_sequence_Q9NU22,
#         conservation_seq=conservation_seq,
#         conservation_key="conservation_key",
#         conservations=conservations,
#         conservations_positions=["1", "3"],
#         target_hit_start=1,
#         target_hit_end=3,
#         conservation_start=1,
#         conservation_end=3,
#     )
#     assert transfer_dict["PF00001"]["sequence_id"]["target"]["conservations"]["positions"]["1"] == {
#         "conservation": {"score": 5},
#         "hit": {True},
#     }


###T populate_go_data_for_annotations



###T cleanup_improve_transfer_dict


# @pytest.mark.usefixtures("mock_read_conservations_and_annotations")
# def test_cleanup_improve_transfer_dict(
#     minimal_transfer_dict,
#     mock_iterate_aligned_sequences,
#     monkeypatch
# ):
#     """
#     Tests cleanup_improve_transfer_dict moves 'DOMAIN' to [pfam_id], reads data, etc.
#     """
#     pfam_id = "PFTEST"
#     # Mocks
#     monkeypatch.setattr("transfer_annotations.get_alignment_sequences", lambda *args, **kw: ("AAA", "AAA"))
#     monkeypatch.setattr("transfer_annotations.gather_go_terms_for_target", lambda *args, **kw: {"GO:000xyz"})

#     result = cleanup_improve_transfer_dict(
#         transfer_dict=minimal_transfer_dict,
#         pfam_id=pfam_id,
#         hmmalign_lines=[],
#         conservations_filepath="fake.json",
#         annotations_filepath="fake.json",
#         output_dir="fake_dir",
#     )

#     assert pfam_id in result
#     assert "DOMAIN" not in result
#     target_data = result[pfam_id]["sequence_id"]["TestTarget"]
#     assert "conservations" in target_data
#     assert "annotations" in target_data


###T convert_sets_to_lists


###T write_reports
def test_write_reports_with_data(
    tmp_path, logger,
    transfer_dict_populated_disulfid_post_processed_Q9NU22
):
    """Test write_reports with actual annotation data"""
    mock_dict = transfer_dict_populated_disulfid_post_processed_Q9NU22

    with patch('builtins.open', mock_open()) as mock_file:
        write_reports(logger, mock_dict, str(tmp_path))

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

def test_write_reports_empty_data(tmp_path, logger):
    """Test write_reports with empty transfer dict"""
    transfer_dict = {
        "PF07728": {
            "sequence_id": {
                "sp|Q9NU22|MDN1_HUMAN": {}
            }
        }
    }

    with patch('builtins.open', mock_open()) as mock_file:
        write_reports(logger, transfer_dict, str(tmp_path))

        calls = mock_file.call_args_list
        assert len(calls) == 2  # Main report + target report
        assert logger.debug.call_count == 3

def test_write_reports_multiple_targets(tmp_path, logger):
    """Test write_reports with multiple target sequences"""
    transfer_dict = {
        "PF07728": {
            "sequence_id": {
                "TEST1": {"some_data": "value1"},
                "TEST2": {"some_data": "value2"}
            }
        }
    }

    with patch('builtins.open', mock_open()) as mock_file:
        write_reports(logger, transfer_dict, str(tmp_path))

        calls = mock_file.call_args_list
        assert len(calls) == 3  # Main report + 2 target reports
        assert logger.debug.call_count == 5  # 1 main + 2 data writes + 2 confirmations

def test_write_reports_pipe_in_target_name(tmp_path, logger):
    """Test write_reports handles target names with pipes"""
    transfer_dict = {
        "PF07728": {
            "sequence_id": {
                "TEST|A": {"some_data": "value"}
            }
        }
    }

    with patch('builtins.open', mock_open()) as mock_file:
        write_reports(logger, transfer_dict, str(tmp_path))

        calls = mock_file.call_args_list
        assert len(calls) == 2  # Main report + target report

        # Check paths
        paths = [call[0][0] for call in calls]
        assert any("TEST-A" in path for path in paths)
        assert logger.debug.call_count == 3

def test_write_reports_completely_empty_dict(tmp_path, logger):
    """Test write_reports with completely empty transfer dict"""
    transfer_dict = {}

    with pytest.raises(StopIteration):
        write_reports(logger, transfer_dict, str(tmp_path))

def test_write_reports_empty_sequence_id(tmp_path, logger):
    """Test write_reports with empty sequence_id dict"""
    transfer_dict = {
        "PF07728": {
            "sequence_id": {}
        }
    }

    with patch('builtins.open', mock_open()) as mock_file:
        write_reports(logger, transfer_dict, str(tmp_path))

        calls = mock_file.call_args_list
        assert len(calls) == 1  # Only main report
        assert logger.debug.call_count == 1

def test_write_reports_none_target_data(tmp_path, logger):
    """Test write_reports with None as target data"""
    transfer_dict = {
        "PF07728": {
            "sequence_id": {
                "TEST": None
            }
        }
    }

    with patch('json.dump') as mock_json_dump:
        write_reports(logger, transfer_dict, str(tmp_path))

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


###T main
def test_main_success(logger, minimal_hmmalign_lines_fixture_Q9NU22, transfer_dict_populated_disulfid_list_Q9NU22, transfer_dict_populated_disulfid_post_processed_list_Q9NU22):
    """Test main function success path"""
    mock_args = Namespace(
        dom_align=hmmalign_result_mock,
        resource_dir=resource_dir_mock,
        output_dir=output_dir_mock,
        eco_codes=good_eco_codes_mock,
        log=log_filepath_mock
    )

    with patch('transfer_annotations.parse_arguments', return_value=mock_args) as mock_parse, \
         patch('transfer_annotations.get_pfam_id_from_hmmalign_result', return_value="PF07728") as mock_get_pfam, \
         patch('transfer_annotations.read_files', return_value=(minimal_hmmalign_lines_fixture_Q9NU22, {})) as mock_read, \
         patch('transfer_annotations.find_and_map_annots', return_value=transfer_dict_populated_disulfid_list_Q9NU22) as mock_find, \
         patch('transfer_annotations.cleanup_improve_transfer_dict', return_value=transfer_dict_populated_disulfid_post_processed_list_Q9NU22) as mock_cleanup, \
         patch('transfer_annotations.write_reports') as mock_write:

        main(logger)

        mock_parse.assert_called_once()
        mock_get_pfam.assert_called_once_with(hmmalign_result_mock)
        mock_read.assert_called_once()
        mock_find.assert_called_once()
        mock_cleanup.assert_called_once_with(
            transfer_dict_populated_disulfid_list_Q9NU22,
            "PF07728",
            minimal_hmmalign_lines_fixture_Q9NU22,
            os.path.join(resource_dir_mock, "PF07728", "conservations.json"),
            os.path.join(resource_dir_mock, "PF07728", "annotations.json"),
            output_dir_mock
        )
        mock_write.assert_called_once_with(
            logger,
            transfer_dict_populated_disulfid_post_processed_list_Q9NU22,
            output_dir_mock
        )
        logger.info.assert_any_call("---> MAIN --- Transfer Dict FILLED")




# Integration tests
def test_main_integration_binding_Q9NU22_PF07728(minimal_hmmalign_lines_fixture_Q9NU22, transfer_dict_populated_binding_Q9NU22, annotations_content_binding_fixture_Q9NU22_PF07728, conservations_content_Q9NU22_PF07728):
    """Integration test for main function"""

    with TemporaryDirectory() as tmp_dir:
        # Setup test files and directories
        hmmalign_path = os.path.join(tmp_dir, "PF07728_hmmalign.sth")
        with open(hmmalign_path, 'w', encoding='latin-1') as f:
            f.write(minimal_hmmalign_lines_fixture_Q9NU22)

        # Create full directory structure
        resource_dir = os.path.join(tmp_dir, "resources")
        pfam_dir = os.path.join(resource_dir, "PF07728")
        os.makedirs(pfam_dir)

        output_dir = os.path.join(tmp_dir, "output")
        # Create output directories for expected targets
        target_dirs = [
            os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN"),
            os.path.join(output_dir, "PF07728")
        ]
        for dir_path in [output_dir] + target_dirs:
            os.makedirs(dir_path, exist_ok=True)

        # Write annotations to proper path
        with open(os.path.join(pfam_dir, "annotations.json"), 'w', encoding='utf-8') as f:
            json.dump(annotations_content_binding_fixture_Q9NU22_PF07728, f)

        with open(os.path.join(pfam_dir, "conservations.json"), 'w', encoding='utf-8') as f:
            json.dump(conservations_content_Q9NU22_PF07728, f)

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

        with open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json")) as f:
            transfer_dict = json.load(f)
            print(json.dumps(transfer_dict, indent=4))

        # assert os.path.exists(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json"))
        # assert transfer_dict_success_binding_Q9NU22 == json.load(open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json"), 'r', encoding='utf-8'))

def test_main_integration_disulfid_Q9NU22_PF07728(minimal_hmmalign_lines_fixture_Q9NU22, transfer_dict_populated_disulfid_post_processed_list_Q9NU22, annotations_content_disulfid_fixture_Q9NU22_PF07728, conservations_content_Q9NU22_PF07728):
    """Integration test for main function"""

    with TemporaryDirectory() as tmp_dir:
        # Setup test files and directories
        hmmalign_path = os.path.join(tmp_dir, "PF07728_hmmalign.sth")
        with open(hmmalign_path, 'w', encoding='latin-1') as f:
            f.write(minimal_hmmalign_lines_fixture_Q9NU22)

        # Create full directory structure
        resource_dir = os.path.join(tmp_dir, "resources")
        pfam_dir = os.path.join(resource_dir, "PF07728")
        os.makedirs(pfam_dir)

        output_dir = os.path.join(tmp_dir, "output")
        # Create output directories for expected targets
        target_dirs = [
            os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN"),
            os.path.join(output_dir, "PF07728")
        ]
        for dir_path in [output_dir] + target_dirs:
            os.makedirs(dir_path, exist_ok=True)

        # Write annotations to proper path
        with open(os.path.join(pfam_dir, "annotations.json"), 'w', encoding='utf-8') as f:
            json.dump(annotations_content_disulfid_fixture_Q9NU22_PF07728, f)

        with open(os.path.join(pfam_dir, "conservations.json"), 'w', encoding='utf-8') as f:
            json.dump(conservations_content_Q9NU22_PF07728, f)

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

        with open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json")) as f:
            transfer_dict = json.load(f)
            print("AYOOOOOOOO &&&&& ")
            print(json.dumps(transfer_dict, indent=4))

        assert os.path.exists(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json"))
        assert transfer_dict_populated_disulfid_post_processed_list_Q9NU22 == json.load(open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF07728_report.json"), 'r', encoding='utf-8'))

def test_main_integration_all_types_H0YB80(minimal_hmmalign_lines_fixture_H0YB80, transfer_dict_success_all_types_H0YB80, annotations_content_all_types_fixture_H0YB80_PF00244, conservations_content_H0YB80_PF00244):
    """Integration test for main function"""

    with TemporaryDirectory() as tmp_dir:
        # Setup test files and directories
        hmmalign_path = os.path.join(tmp_dir, "PF00244_hmmalign.sth")
        with open(hmmalign_path, 'w', encoding='latin-1') as f:
            f.write(minimal_hmmalign_lines_fixture_H0YB80)

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
            json.dump(annotations_content_all_types_fixture_H0YB80_PF00244, f)

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

        # with open(os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN", "PF00244_report.json")) as f:
        #     transfer_dict = json.load(f)
        #     print(json.dumps(transfer_dict, indent=4))

        assert os.path.exists(os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN", "PF00244_report.json"))
        assert transfer_dict_success_all_types_H0YB80 == json.load(open(os.path.join(output_dir, "tr-H0YB80-H0YB80_HUMAN", "PF00244_report.json"), 'r', encoding='utf-8'))



        # with open(os.path.join(output_dir, "sp-Q9NU22-MDN1_HUMAN", "PF00244_report.json")) as f:
        #     transfer_dict = json.load(f)
        #     print(json.dumps(transfer_dict, indent=4))
