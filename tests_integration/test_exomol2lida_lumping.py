"""
These are integration tests ensuring that the chunk_sizes and number of .trans files
do not have any effect on the states and transitions lumping outputs.
"""

from pathlib import Path

import pytest

from exomol2lida.read_inputs import MoleculeInput
from exomol2lida.process_dataset import DatasetProcessor

test_resources_dir = Path(__file__).parent / "resources"

mol_input = MoleculeInput(
    molecule_formula="FOO",
    **{
        "mol_slug": "HCN",
        "iso_slug": "1H-12C-14N",
        "dataset_name": "Harris",
        "states_header": [
            "i",
            "E",
            "g_tot",
            "J",
            "+/-",
            "kp",
            "iso",
            "v1",
            "v2",
            "l2",
            "v3",
        ],
        "resolve_vib": ["v1", "v2", "v3"],
        "energy_max": 9000,
        "only_with": {"iso": "1"},
    },
)
# test files are built from a real dataset, but were skimmed down to manageable size
# of 500,000 transitions (with split files of 100,000 transitions each),
# and only states involved in the transitions left are present in the states file.
states_path = test_resources_dir / "dummy_data.states.bz2"
trans_paths_full = [test_resources_dir / "dummy_data.trans.bz2"]
trans_paths_split = sorted(test_resources_dir.glob("dummy_data.trans_0*.bz2"))

shared_for_comparison = {
    "lumped_states": None,
    "states_map_lumped_to_original": None,
    "states_map_original_to_lumped": None,
    "lumped_states_lifetimes": [
        None,
    ],
    "lumped_transitions": None,
}


@pytest.mark.parametrize("chunk_size", (1_000_000, 100_000, 10_000, 5_000))
def test_states_lumping(monkeypatch, chunk_size):
    processor = DatasetProcessor(molecule=mol_input)
    processor.include_original_lifetimes = True
    monkeypatch.setattr(processor, "states_path", states_path)
    processor.states_chunk_size = chunk_size
    processor.lump_states()
    if shared_for_comparison["lumped_states"] is None:
        shared_for_comparison["lumped_states"] = processor.lumped_states.copy(deep=True)
        shared_for_comparison[
            "states_map_lumped_to_original"
        ] = processor.states_map_lumped_to_original.copy()
        shared_for_comparison[
            "states_map_original_to_lumped"
        ] = processor.states_map_original_to_lumped.copy()
    else:
        assert processor.lumped_states.equals(shared_for_comparison["lumped_states"])
        assert (
            processor.states_map_lumped_to_original
            == shared_for_comparison["states_map_lumped_to_original"]
        )
        assert (
            processor.states_map_original_to_lumped
            == shared_for_comparison["states_map_original_to_lumped"]
        )


@pytest.mark.parametrize(
    "trans_paths, chunk_size",
    (
        (trans_paths_full, 1_000_000),
        (trans_paths_full, 100_000),
        (trans_paths_full, 10_000),
        (trans_paths_split, 100_000),
        (trans_paths_split, 10_000),
        (trans_paths_split, 5_000),
    ),
)
def test_trans_lumping(monkeypatch, trans_paths, chunk_size):
    # prepare
    processor = DatasetProcessor(molecule=mol_input)
    processor.include_original_lifetimes = True
    monkeypatch.setattr(processor, "states_path", states_path)
    processor.states_chunk_size = 1_000_000
    processor.lump_states()
    # run tests
    monkeypatch.setattr(processor, "trans_paths", trans_paths)
    processor.trans_chunk_size = chunk_size
    processor.lump_transitions()
    if shared_for_comparison["lumped_states_lifetimes"] == [
        None,
    ]:
        shared_for_comparison["lumped_states_lifetimes"] = list(
            processor.lumped_states["tau"]
        )
        shared_for_comparison["lumped_transitions"] = processor.lumped_transitions.copy(
            deep=True
        )
    else:
        assert (
            list(processor.lumped_states["tau"])
            == shared_for_comparison["lumped_states_lifetimes"]
        )
        assert processor.lumped_transitions.equals(
            shared_for_comparison["lumped_transitions"]
        )
