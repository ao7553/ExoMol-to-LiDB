"""
Module reading the dict from *input/molecules.py* input file and performing some
validation and sanitization on top of the input configuration data.

The ``molecules.data`` is expected to hold the following entries:
"molecule_formula": {
  "mol_slug": str,
  "iso_slug": str,
  "dataset_name": str,
  "states_header": list[str] (optional, bypasses .def file if specified),
  "resolve_el": list[str] (optional if `resolve_vib` not specified),
  "resolve_vib": list[str] (optional if `resolve_el` not specified),
  "energy_max": number (optional), in [eV],
  "only_with": dict[str, str] (optional),
  "only_without": dict[str, str] (optional),
}

`molecule_formula`: Identifier for the Lida database, does not have to correspond
    to the `mol_slug` or `iso_slug`. For example ``"HD+"`` might be used
    instead of ``"H2+"``. Must be, however, parsable by `PyValem` package.
`mol_slug`, `iso_slug`, `dataset_name`: Mandatory attributes for each
    `molecule_formula`, these three attributes will determine where to look for the
    data inside the `EXOMOL_DATA_DIR`.
`states_header`: The list of column names for the .states data table. If not supplied,
    the .def file for the given dataset is parsed and the .states column names are
    inferred from there. It should always be preferred to keep the .def files in order
    rather than patching inconsistent .def files by this higher-level fix.
`resolve_el`, `resolve_vib`: The lists of .states columns which need to be
    *resolved*. All the other columns will be lumped over and the resolved ones will
    stay in the output composite states and their transitions.
`energy_max`: Optional maximal energy in [eV]. States with higher energy will simply
    be filtered out and with them all their transitions (to and from).
`only_with`: This attribute allows for adjustable pre-filtering of the .states and
    .trans files. Keys of the `only_with` dict are .states column names and values
    are the only values allowed in these .states columns. As an example,
    ``"HCN": {..., "only_with": {"iso": "HNC"}, ...}`` will pre-filter the dataset and
    remove all the states (and their transitions) with ``"iso"`` column values
    (in .states file) not equaling ``"HNC"``.
`only_without`: Allows for similar pre-filtering as the `only_with`. As an example,
    the VO .states contains some values ``"0"`` in the ``States`` column, probably
    indicating unassigned states. These will be ignored by setting
    ``"only_without": {"State": "0"}``

A `MoleculeInput` class instantiated without exceptions signals data without
inconsistencies and ready to be processed into the Lida data product. All the possible
consistency checks are performed upon MoleculeInput instantiation.
The `MoleculeInput` instance will contain all the original ``molecules.data`` fields per
each `molecule_formula` saved as instance attributes, as well as couple of additional
attributes, such as `self.def_path`, `self.states_path` and `self.trans_paths`.
"""

from pathlib import Path

from exomole.read_def import DefParser, DefParseError
from exomole.utils import get_num_columns

from config.config import EXOMOL_DATA_DIR
from .exceptions import MoleculeInputError
from .utils import EV_IN_CM


class MoleculeInput:
    """Class representing Molecule Inputs.

    Loads all the attributes for a single `molecule_formula` from the
    ``input.molecules.molecules``, does whole lot of validation and preparation for the
    states lumping.

    Parameters
    ----------
    molecule_formula : str
    kwargs : dict, optional
        If not supplied, the class will load all the necessary arguments from the
        ``input.molecules.molecules``.
        If omitted, all the necessary arguments (see the structure of ``molecules.data``
        input file discussed in the docstring for `read_inputs` module) need to be
        passed explicitly.

    Attributes
    ----------
    formula : str
    iso_formula : str
    mol_slug : str
    iso_slug : str
    dataset_name : str
    states_header : list[str], optional
    resolve_el : list[str]
    resolve_vib : list[str]
    energy_max : float
        This is in [cm-1], converted from the input file.
    only_with : dict[str, str]
    def_path : Path
    states_path : Path
    trans_paths : list[Path]
    def_parser : DefParser
    def_parser_raised : Exception, optional
    version : int
    mass : float

    Raises
    ------
    MoleculeInputError
        If any inconsistencies are detected in the molecule input.
    DefParseError
        If the `MoleculeInput` relies on the information from the .def file (such as
        if the states_header is not explicitly provided in the input file) and the .def
        file cannot be parsed, this error is raised.
    """

    def __init__(self, molecule_formula, **kwargs):
        self.formula = molecule_formula
        self.raw_input = None
        self.def_parser = None
        self.def_parser_raised = None
        self.iso_formula = None
        self.version = None
        self.mass = None

        # stubs for all the arguments which might be expected:
        self.mol_slug = None
        self.iso_slug = None
        self.dataset_name = None
        self.states_header = None
        self.resolve_el = []
        self.resolve_vib = []
        self.energy_max = float("inf")
        self.only_with = {}
        self.only_without = {}

        # populate the attributes:
        if not len(kwargs):
            from input.molecules import molecules as molecules

            kwargs = molecules[molecule_formula]

        self.raw_input = kwargs.copy()
        for attr, val in kwargs.items():
            setattr(self, attr, val)
        self.energy_max = float(self.energy_max) * EV_IN_CM
        if not all([self.mol_slug, self.iso_slug, self.dataset_name]):
            raise MoleculeInputError(
                f"Input data for {molecule_formula} missing some of the mandatory "
                f'attributes ("mol_slug", "iso_slug", "dataset_name")'
            )
        if not any([len(self.resolve_el), len(self.resolve_vib)]):
            raise MoleculeInputError(
                f"Input data for {molecule_formula} missing one of the 'resolve_el', "
                f"'resolve_vib' attributes."
            )

        mol_root = EXOMOL_DATA_DIR / self.mol_slug
        if not mol_root.is_dir():
            raise MoleculeInputError(f"The molecule directory not found: {mol_root}")
        iso_root = mol_root / self.iso_slug
        if not iso_root.is_dir():
            raise MoleculeInputError(
                f"The isotopologue directory not found: {iso_root}"
            )
        ds_root = iso_root / self.dataset_name
        if not ds_root.is_dir():
            raise MoleculeInputError(f"The dataset directory not found: {ds_root}")

        file_name_stem = f"{self.iso_slug}__{self.dataset_name}"

        self.def_path = ds_root / f"{file_name_stem}.def"
        if not self.def_path.is_file():
            raise MoleculeInputError(f"The .def file not found under {self.def_path}")

        for states_path in [
            ds_root / f"{file_name_stem}.states.bz2",
            ds_root / f"{file_name_stem}.states",
        ]:
            if states_path.is_file():
                self.states_path = states_path
                break
        else:
            raise MoleculeInputError(f"No .states file found in {ds_root}")

        for trans_wc in [f"{file_name_stem}*.trans.bz2", f"{file_name_stem}*.trans"]:
            trans_paths = sorted(ds_root.glob(trans_wc))
            if len(trans_paths):
                self.trans_paths = trans_paths
                break
        else:
            raise MoleculeInputError(f"No .trans files found in {ds_root}")

        # try to parse the .def file as far as I can get.
        self.def_parser = DefParser(self.def_path)
        try:
            self.def_parser.parse(warn_on_comments=False)
        except DefParseError as e:
            self.def_parser_raised = e
        # We will need some data from the .def file for the outputs:
        if (
            self.def_parser.iso_formula is None
            or self.def_parser.version is None
            or self.def_parser.mass is None
        ):
            raise self.def_parser_raised
        else:
            self.iso_formula = self.def_parser.iso_formula
            self.version = self.def_parser.version
            self.mass = self.def_parser.mass

        # get .states column names:
        if self.states_header is not None:
            # some basic sanitation of the states header read from input:
            if self.states_header[:4] != ["i", "E", "g_tot", "J"]:
                raise MoleculeInputError(
                    f"Unexpected states_header for {molecule_formula}"
                )
        else:
            # states header is not explicitly specified in the input, get it from
            # the parsed .def file:
            if self.def_parser.quanta is None:
                # I need the quanta and few attributes before quanta to construct the
                # states header. If they are not there, that means the DefParser did
                # not finish cleanly.
                raise self.def_parser_raised
            self.states_header = self.def_parser.get_states_header()

        # resolve_el and resolve_vib must not share any states:
        if set(self.resolve_el).intersection(self.resolve_vib):
            raise MoleculeInputError(
                f"Common values found in 'resolve_el' and 'resolve_vib' for "
                f"{molecule_formula}."
            )
        resolved_quanta = set(self.resolve_el + self.resolve_vib)
        # only states might be in resolve_vib and resolve_el:
        if resolved_quanta.intersection(["i", "E", "g_tot", "J", "tau", "g_J"]):
            raise MoleculeInputError(
                f"Unsupported values found in 'resolve_el' or 'resolve_vib' for "
                f"{molecule_formula}."
            )
        quanta_available = set(self.states_header).difference(
            ["i", "E", "g_tot", "J", "tau", "g_J"]
        )
        if not resolved_quanta.issubset(quanta_available):
            raise MoleculeInputError(
                f"Unsupported values found in 'resolve_el' or 'resolve_vib' for "
                f"{molecule_formula}."
            )
        # "only_with" keys need to be subset of quanta nad J
        if not set(self.only_with).issubset(quanta_available | {"J"}):
            raise MoleculeInputError(
                f"Unrecognised 'only_with': {self.only_with} not among "
                f"quanta available in {self.dataset_name}."
            )
        # isomers failsafe: if "iso" among quanta, it needs to be either resolved, or
        # filtered via only_with parameter, so different isomers are not lumped
        if "iso" in quanta_available:
            if "iso" not in self.only_with and "iso" not in resolved_quanta:
                raise MoleculeInputError(
                    f"Cannot lump over 'iso' states. Isomers need to be either "
                    f"resolved or filtered with 'only_with' parameter for "
                    f"{molecule_formula}"
                )

        # check if the states header aligns with the .states file in number of
        # columns in .states:
        states_num_columns = get_num_columns(self.states_path)
        if len(self.states_header) != states_num_columns:
            msg = (
                f"{self.states_path.name} has {states_num_columns} "
                f"columns, while input or {self.def_path.name} specifies "
                f"{len(self.states_header)} columns."
            )
            raise MoleculeInputError(msg)

        # finally, check if the .trans file has the appropriate number of columns:
        num_columns_trans = get_num_columns(self.trans_paths[0])
        if num_columns_trans not in {3, 4}:
            msg = (
                f"{self.trans_paths[0].name} has {num_columns_trans} "
                f"columns, while 3 or 4 columns are expected!"
            )
            raise MoleculeInputError(msg)


def get_all_inputs(bypass_exceptions=False, verbose=True):
    """Get the `MoleculeInput` instances for all formulas specified in the input
    file.

    Parameters
    ----------
    bypass_exceptions : bool, optional
    verbose : bool, optional

    Returns
    -------
    dict[str, Optional[MoleculeInput]]
    """
    from input.molecules import molecules as inputs_dict

    all_inputs = {}
    num_exceptions_raised = 0
    for molecule_formula in inputs_dict:
        if not bypass_exceptions:
            mol_input = MoleculeInput(molecule_formula, **inputs_dict[molecule_formula])
        else:
            try:
                mol_input = MoleculeInput(
                    molecule_formula, **inputs_dict[molecule_formula]
                )
            except (MoleculeInputError, DefParseError) as e:
                num_exceptions_raised += 1
                mol_input = None
                if verbose:
                    print(f"{molecule_formula}: {e}")
        if mol_input is not None:
            all_inputs[molecule_formula] = mol_input
    if bypass_exceptions and num_exceptions_raised and verbose:
        print(
            f"{num_exceptions_raised}/{len(inputs_dict)} inconsistent inputs detected"
        )
    return all_inputs
