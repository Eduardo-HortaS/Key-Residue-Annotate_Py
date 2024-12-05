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
from unittest.mock import patch, mock_open, MagicMock
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
    _add_single_annotation,
    remove_failed_annotations,
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
def annotations_content_disulfid_fixture():
    """Separate fixture for DISULFID pair annotations"""
    return {
        "MCRB_ECOLI": {
            "205": [
                {
                    "type": "DISULFID",
                    "description": "Interchain (with C-246); in linked form",
                    "evidence": "ECO:0000269|PubMed:12345678",
                    "entry": "P15005",
                    "aminoacid": "C",
                    "paired_position": "246"
                }
            ],
            "246": [
                {
                    "type": "DISULFID",
                    "description": "Interchain (with C-205); in linked form",
                    "evidence": "ECO:0000269|PubMed:12345678",
                    "entry": "P15005",
                    "aminoacid": "C",
                    "paired_position": "205"
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
                "description": "Interchain (with C-246); in linked form",
                "evidence": "ECO:0000269|PubMed:12345678",
                "entry": "P15005",
                "aminoacid": "C", # NIILQGPPGC
                "paired_position": "246"
            }],
        "246": [
            {
                "type": "DISULFID",
                "description": "Interchain (with C-205); in linked form",
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
                "description": "Interchain (with C-246); in linked form",
                "evidence": "ECO:0000269|PubMed:12345678",
                "entry": "P15005",
                "aminoacid": "C", # NIILQGPPGC
                "paired_position": "246",
                "target_position": "333"
            }],
        "246": [
            {
                "type": "DISULFID",
                "description": "Interchain (with C-205); in linked form",
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
def annot_pos_paireable_type_hit_bools():
    return {"201": {"BINDING": False}, "202": {"BINDING": False}, '205': {'DISULFID': False}, '246': {'DISULFID': False}}

@pytest.fixture
def transfer_dict():
    return {}

@pytest.fixture
def transfer_fict_success_binding():
    return {
    "sp|Q9NU22|MDN1_HUMAN": {
        'match': {
            329: {
                'BINDING | Interacts with GTP': {
                    'essentials': {
                        'type': 'BINDING',
                        'description': 'Interacts with GTP',
                        'count': 1
                    },
                    'evidence': {
                        'ECO:0000255': {
                            'rep_mnemo_name': 'MCRB_ECOLI',
                            'rep_primary_accession': 'P15005',
                            'count': 1
                        }
                    },
                    'additional_keys': {
                        'ligand_id': {
                            'count': 1,
                            'rep_mnemo_name': 'MCRB_ECOLI',
                            'rep_primary_accession': 'P15005',
                        }
                    }
                }
            },
            330: {
                'BINDING | Interacts with GTP': {
                    'essentials': {
                        'type': 'BINDING',
                        'description': 'Interacts with GTP',
                        'count': 1
                    },
                    'evidence': {
                        'ECO:0000255': {
                            'rep_mnemo_name': 'MCRB_ECOLI',
                            'rep_primary_accession': 'P15005',
                            'count': 1
                        }
                    },
                    'additional_keys': {
                        'ligand_id': {
                            'count': 1,
                            'rep_mnemo_name': 'MCRB_ECOLI',
                            'rep_primary_accession': 'P15005',
                        }
                    }
                }
            }
        },
        'miss': {}
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
def annotation_dict_post_make_anno_total(annotation_dict_205):
    annotation_dict = copy.deepcopy(annotation_dict_205)
    annotation_dict['target_position'] = 333
    annotation_dict['paired_position'] = 246
    return annotation_dict

@pytest.fixture
def paired_annotation_dict_post_map_filter(annotation_dict_246):
    annotation_dict = copy.deepcopy(annotation_dict_246)
    annotation_dict['target_position'] = 205
    annotation_dict['paired_position'] = 373
    return annotation_dict

@pytest.fixture
def mock_make_anno_total_disulfid_return(annotation_dict_post_make_anno_total):
    return {
        'annotation': annotation_dict_post_make_anno_total,
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Interchain (with C-246); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Interchain (with C-246); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678'
        },
        'paired_position': '246'
    }

@pytest.fixture
def mock_map_filter_disulfid_return(paired_annotation_dict_post_map_filter):
    return (True, {
        'annotation': paired_annotation_dict_post_map_filter,
        'anno_type': 'DISULFID',
        'anno_id': 'DISULFID | Interchain (with C-205); in linked form',
        'anno_total': {
            'type': 'DISULFID',
            'description': 'Interchain (with C-205); in linked form',
            'count': 1,
            'evidence': 'ECO:0000269|PubMed:12345678'
        },
        'paired_position': '205'
    })


###### Configuration Constants
### Added in Map_and_filter_annot_pos
target_name_mock = "sp|Q9NU22|MDN1_HUMAN"
target_hit_start_mock = 325
target_hit_end_mock = 451
offset_start_mock = 196
offset_end_mock = 350
entry_mnemo_name_mock = "MCRB_ECOLI"

## Extra variables for paired position testing
target_paired_uniprot_pos_mock = "246"
caller_target_pos_mock = 333
# annotation_dict = use the get_annotation_dict fixture
## Result for running with variables for paired positions
paired_position_res_hit_plus_late_pair_result_dict_tuple_mock = (True, {'annotation': {'type': 'DISULFID', 'description': 'Interchain (with C-205); in linked form', 'evidence': 'ECO:0000269|PubMed:12345678', 'entry': 'P15005', 'aminoacid': 'C', 'paired_position': '333', 'target_position': 373}, 'anno_type': 'DISULFID', 'anno_id': 'DISULFID | Interchain (with C-205); in linked form', 'anno_total': {'type': 'DISULFID', 'description': 'Interchain (with C-205); in linked form', 'count': 1, 'evidence': 'ECO:0000269|PubMed:12345678', 'paired_position': '333'}, 'paired_position': '333'})


### Added Defaults for use in validate_annotations and validate_paired_annotations
### (values change in downstream functions from these cited)
counter_target_pos_mock = None
counter_uniprot_pos_mock = None

### Added/Modified in Process_annotation
# res_hit = Define in the test
# annotation_dict = use the get_annotation_dict fixture
# counter_target_pos = Define in the test
# counter_uniprot_pos_str = Define in the test

### Added/Modified in Make_anno_total_dict
# annotation_dict = use the get_annotation_dict fixture
# counter_target_pos = Define in the test
# counter_uniprot_pos_str = Define in the test
# annot_pos_paireable_type_hit_bools = use fixture with the same name

## If called from validate_paired_annotations:
# caller_target_pos = Define in the test
# logger = use the logger fixture
# entry_annotations = use the entry_annotations_* fixture, according to the test

### Added/Modified in Add_to_transfer_dict
entry_primary_accession_mock = "P15005"
# counter_target_pos = Define in the test
# anno_id = Define in the test, new
# anno_total = Define in the test, new

## If called from a paireable type in process_annotation, also note that:
# target_paired_uniprot_pos = Define in the test
# caller_target_pos = Define in the test
# annotation_dict = use the get_annotation_dict fixture

### Added/Modified in _add_single_annotation
# paired_additional_keys/additional_keys = Define in the test
# paired_position_res_hit/hit = Define in the test ???
# late_pair_anno_id = Define in the test ???
# late_pair_anno_total = Define in the test ???

### Added/Modified in Remove_failed_annotations
# entry_annotations = use the entry_annotations_* fixture, according to the test
# counter_unipro_pos_str = Define in the test
# paired_position = Define in the test
# annotation_type = Define in the test


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
        "-e", *good_eco_codes_mock, # Unpacking list of ECO codes to simulate space-separated input
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



def test_find_and_map_annots_basic_case(annotations_content_binding_fixture, logger, target_sequence, annot_sequence):
    """Test basic case with single binding annotation"""
    mock_lines = [
        "# STOCKHOLM 1.0\n",
        "MCRB_ECOLI/196-350      .........NIILQGPPGCGKTFVAR.RLAYLLtg.........ekaPQRVNMVQFHQ.SY..SYEDFIQGYCPN............GVGFRRKDGIFYNFCQqak...........eqpeKKYIFIIDEINRANlS........K....VFG....E......VMM.LMEHD.K...RGENWSVPLTYSENDE.eRFYVp..........................enVYIIGLMNTAD.........RSLAVVDYALRRRF............................\n",
        "sp|Q9NU22|MDN1_HUMANtarget//325-451                     .........-VLLEGPIGCGKTSLVE.YLAAVTgr..........tkPPQLLKVQLGD.QT..DSKMLLGMYCCTd..........vPGEFVWQPGTLTQAAT..................MGHWILLEDIDYAP.L........D....VVS....V......LIP.LLENG.E...LLIPGRGDCLKVAPG-..----.............................FQFFAT-----.........--------------rrllscgg....................\n"
    ]

    with patch('transfer_annotations.map_and_filter_annot_pos') as mock_map_and_filter:
        find_and_map_annots(
            logger,
            mock_lines,
            annotations_content_binding_fixture,
            []  # No ECO code filtering
        )

        # assert transfer_dict == transfer_fict_success_binding, "Transfer dictionary does not match the expected structure."
        mock_map_and_filter.assert_called_with(
            logger,
            [],  # Empty good_eco_codes
            target_sequence,
            target_name_mock,
            target_hit_start_mock,  # target_hit_start
            target_hit_end_mock,  # target_hit_end
            offset_start_mock,  # offset_start
            offset_end_mock,  # offset_end
            annot_sequence,
            entry_mnemo_name_mock,
            annotations_content_binding_fixture["MCRB_ECOLI"],
            {},  # Initial transfer_dict
            {'201': {'BINDING': False}, '202': {'BINDING': False}},  # Initial annot_pos_paireable_type_hit_bools
        )

def test_find_and_map_annots_no_target_sequence(annotations_content_binding_fixture, logger):
    """Test case with no matching annotations"""
    mock_lines = [
        "# STOCKHOLM 1.0\n",
        "VWA8_HUMAN/105-261    .........DVFLIGPPGPLRRSIAM.QYLELT..............KREVEYIALSR.DT..TETDLKQRREIR\n"
    ]

    transfer_dict = find_and_map_annots(
        logger,
        mock_lines,
        annotations_content_binding_fixture,
        []
    )

    assert transfer_dict == {}

    # Verify the logger captured the correct error message
    logger.error.assert_called_once_with(
        "---> DEBUG --- FIND_AND_MAP --- No target sequences found in hmmalign lines!"
    )

def test_map_and_filter_with_paired_positions(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_disulfid_pair, transfer_dict,
    annot_pos_paireable_type_hit_bools, get_annotation_dict,
    annotations_content_disulfid_fixture
):
    annotation_dict = get_annotation_dict(annotations_content_disulfid_fixture, "246")
    # Scenario 1: Paired positions provided
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
        annot_pos_paireable_type_hit_bools=annot_pos_paireable_type_hit_bools,
        target_paired_uniprot_pos=target_paired_uniprot_pos_mock,
        caller_target_pos=caller_target_pos_mock,
        annotation_dict=annotation_dict,
    )
    assert result == paired_position_res_hit_plus_late_pair_result_dict_tuple_mock
    logger.info.assert_called()


def test_map_and_filter_without_paired_positions(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_disulfid_pair, transfer_dict,
    annot_pos_paireable_type_hit_bools
):
    # Scenario 2: No paired positions
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
        annot_pos_paireable_type_hit_bools=annot_pos_paireable_type_hit_bools,
    )
    assert result is None  # Modify this based on expected output



def test_validate_annotations_process_annotation_called(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_binding_only, transfer_dict,
    annot_pos_paireable_type_hit_bools
):
    """Test case 1.1 - process_annotation is called when conditions are met"""

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
            annot_pos_paireable_type_hit_bools=annot_pos_paireable_type_hit_bools,
            counter_target_pos=counter_target_pos_mock,
            counter_uniprot_pos=counter_uniprot_pos_mock
        )

        # Assert process_annotation was called
        mock_process.assert_called()

def test_validate_annotations_skip_processed(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_binding_only, transfer_dict
):
    """Test case 1.2 - Skip when annotation was already processed"""

    test_hit_bools = {
        "201": {"BINDING": True},
        "202": {"BINDING": True}
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
            annot_pos_paireable_type_hit_bools=test_hit_bools,
            counter_target_pos=counter_target_pos_mock,
            counter_uniprot_pos=counter_uniprot_pos_mock
        )

        # Verify process_annotation was not called
        mock_process.assert_not_called()
        assert transfer_dict == {}

def test_validate_annotations_no_positions_in_range(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_binding_only, transfer_dict,
    annot_pos_paireable_type_hit_bools
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
            annot_pos_paireable_type_hit_bools=annot_pos_paireable_type_hit_bools,
            counter_target_pos=counter_target_pos_mock,
            counter_uniprot_pos=counter_uniprot_pos_mock
        )

        mock_process.assert_not_called()
        assert transfer_dict == {}

def test_validate_annotations_sequence_end_reached(
    logger, good_eco_codes_all, target_sequence, annot_sequence,
    entry_annotations_binding_only, transfer_dict,
    annot_pos_paireable_type_hit_bools
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
        annot_pos_paireable_type_hit_bools=annot_pos_paireable_type_hit_bools,
        counter_target_pos=None,
        counter_uniprot_pos=None
    )


def test_process_annotation_basic_single(
    logger, good_eco_codes_all, transfer_dict,
    target_sequence, annot_sequence,
    entry_annotations_binding_only, get_annotation_dict,
    annotations_content_binding_fixture,
    annot_pos_paireable_type_hit_bools
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
            'evidence': 'ECO:0000255'
        },
        'paired_position': None
    }

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
            counter_target_pos=329,
            counter_uniprot_pos_str="201",
            annot_pos_paireable_type_hit_bools=annot_pos_paireable_type_hit_bools
        )

        mock_make_total.assert_called_once_with(
            good_eco_codes_all,
            entry_mnemo_name_mock,
            annotation_dict,
            329,
            "201",
            annot_pos_paireable_type_hit_bools,
            logger=logger,
            entry_annotations=entry_annotations_binding_only
        )

        mock_add_dict.assert_called_once_with(
            True,
            logger,
            transfer_dict,
            target_name_mock,
            329,
            "BINDING | Interacts with GTP",
            mock_make_anno_total["anno_total"],
            entry_mnemo_name_mock,
            entry_primary_accession_mock,
        )

        assert annot_pos_paireable_type_hit_bools["201"]["BINDING"] is True


def test_process_annotation_paired_failure(
    logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair,
    target_sequence, annot_sequence,
    mock_make_anno_total_disulfid_return,
    annotation_dict_205, annotation_dict_246,
    annot_pos_paireable_type_hit_bools
):
    """Test processing of paired DISULFID annotations where late pair fails"""

    # Modify target sequence to ensure failure at paired position
    target_sequence_mod_to_fail_72 = list(target_sequence)
    target_sequence_mod_to_fail_72[72] = "R"
    target_sequence_mod_to_fail_72 = "".join(target_sequence_mod_to_fail_72)

    mock_map_filter_fail_return = (True, {})

    with patch('transfer_annotations.make_anno_total_dict', return_value=mock_make_anno_total_disulfid_return) as mock_make_total, \
         patch('transfer_annotations.map_and_filter_annot_pos', return_value=mock_map_filter_fail_return) as mock_map_filter, \
         patch('transfer_annotations.remove_failed_annotations') as mock_remove:
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
            counter_target_pos=333,
            counter_uniprot_pos_str="205",
            annot_pos_paireable_type_hit_bools=annot_pos_paireable_type_hit_bools
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes_all,
            entry_mnemo_name_mock,
            annotation_dict_205,
            333,
            "205",
            annot_pos_paireable_type_hit_bools,
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair
        )

        # Verify map_and_filter_annot_pos was called correctly
        mock_map_filter.assert_called_once_with(
            logger,
            good_eco_codes_all,
            target_sequence_mod_to_fail_72,
            target_name_mock,
            target_hit_start_mock,
            target_hit_end_mock,
            offset_start_mock,
            offset_end_mock,
            annot_sequence,
            entry_mnemo_name_mock,
            entry_annotations_disulfid_pair,
            transfer_dict,
            annot_pos_paireable_type_hit_bools,
            target_paired_uniprot_pos="246",
            caller_target_pos=333,
            annotation_dict=annotation_dict_246
        )

        # Verify remove_failed_annotations was called correctly
        mock_remove.assert_called_once_with(
            entry_annotations_disulfid_pair,
            "205",
            "246",
            "DISULFID"
        )

        # Verify end state
        assert annot_pos_paireable_type_hit_bools["205"]["DISULFID"] is True
        assert annot_pos_paireable_type_hit_bools["246"]["DISULFID"] is True


### NEWsss
def test_process_annotation_paired_success(
    logger, good_eco_codes_all, transfer_dict,
    entry_annotations_disulfid_pair, annotation_dict_205,
    annotation_dict_246, annotation_dict_post_make_anno_total,
    paired_annotation_dict_post_map_filter, mock_make_anno_total_disulfid_return,
    mock_map_filter_disulfid_return, annot_pos_paireable_type_hit_bools, target_sequence, annot_sequence,
):
    """Test processing of paired DISULFID annotations where late pair succeeds"""

    with patch('transfer_annotations.make_anno_total_dict', return_value=mock_make_anno_total_disulfid_return) as mock_make_total, \
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
            counter_target_pos=333,
            counter_uniprot_pos_str="205",
            annot_pos_paireable_type_hit_bools=annot_pos_paireable_type_hit_bools
        )

        # Verify make_anno_total_dict was called correctly
        mock_make_total.assert_called_once_with(
            good_eco_codes_all,
            entry_mnemo_name_mock,
            annotation_dict_205,
            333,
            "205",
            annot_pos_paireable_type_hit_bools,
            logger=logger,
            entry_annotations=entry_annotations_disulfid_pair
        )

        # Verify map_and_filter_annot_pos was called correctly
        mock_map_filter.assert_called_once_with(
            logger,
            good_eco_codes_all,
            target_sequence,
            target_name_mock,
            target_hit_start_mock,
            target_hit_end_mock,
            offset_start_mock,
            offset_end_mock,
            annot_sequence,
            entry_mnemo_name_mock,
            entry_annotations_disulfid_pair,
            transfer_dict,
            annot_pos_paireable_type_hit_bools,
            target_paired_uniprot_pos="246",
            caller_target_pos=333,
            annotation_dict=annotation_dict_246
        )

        # Verify remove_failed_annotations was called correctly
        mock_add_transfer.assert_called_once_with(
            True,
            logger,
            transfer_dict,
            target_name_mock,
            333,
            'DISULFID | Interchain (with C-246); in linked form',
            mock_make_anno_total_disulfid_return['anno_total'],
            entry_mnemo_name_mock,
            entry_primary_accession_mock,
            True,
            'DISULFID | Interchain (with C-205); in linked form',
            mock_map_filter_disulfid_return[1]['anno_total']  # Access the desired anno_total
        )

        # Verify end state
        assert annot_pos_paireable_type_hit_bools["205"]["DISULFID"] is True
        assert annot_pos_paireable_type_hit_bools["246"]["DISULFID"] is True