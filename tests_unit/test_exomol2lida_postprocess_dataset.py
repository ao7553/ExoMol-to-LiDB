import pytest

from exomol2lida.postprocess_dataset import DatasetPostProcessor, CouldNotParseState


def test_default_electronic_state_parser(monkeypatch):
    monkeypatch.setattr(
        DatasetPostProcessor, "__init__", lambda self, mol_formula: None
    )
    dpp = DatasetPostProcessor("foo")
    assert dpp._parse_state_default(["g2Pi"]) == "g(2PI)"
    assert dpp._parse_state_default(["X4Sigma+"]) == "X(4SIGMA+)"
    assert dpp._parse_state_default(["X4Sigma-"]) == "X(4SIGMA-)"
    assert dpp._parse_state_default(["A1Delta"]) == "A(1DELTA)"
    assert dpp._parse_state_default(["Ap4Phi"]) == "A'(4PHI)"


def test_default_electronic_state_parser_fail(monkeypatch):
    monkeypatch.setattr(
        DatasetPostProcessor, "__init__", lambda self, mol_formula: None
    )
    dpp = DatasetPostProcessor("foo")
    with pytest.raises(CouldNotParseState):
        dpp._parse_state_default("g2Pi")
    with pytest.raises(CouldNotParseState):
        dpp._parse_state_default(["g", "2Pi"])
    with pytest.raises(CouldNotParseState):
        dpp._parse_state_default("foo")
