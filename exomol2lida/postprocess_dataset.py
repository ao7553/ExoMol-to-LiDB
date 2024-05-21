"""
This module contains the `DatasetPostProcessor` class and the top-level
`postprocess_molecule` function.

The idea of post-processing the data is to end up with the *exact* form needed by
the LIDA database.

The ``exomol2lida.process_dataset.DatasetProcessor`` will produce several files
containing the data for the LIDA database in a *raw* form.
Five files are produced by the `DatasetProcessor`, which are relevant for the LIDA db:
* `meta_data.json`
* `states_data.csv`
* `states_electronic_raw.csv`
* `states_vibrational.csv`
* `transitions_data.csv`
Most of those files are ready to go to LIDA, as in they are already in the form
compatible with the LIDA database populating logic.

The `DatasetPostProcessor` class potentially creates only the following single output
file:
* `states_electronic.csv` (if `states_electronic_raw.csv` exists).
The file after post-processing will be in the form of a two-column csv, with column
names ``"i"`` and ``"State"``, and with all the values in the ``"State"`` column being
valid pyvalem `MolecularTermSymbol` state strings.
The rules of post-processing are:
* Have a look at ``input.mapping_el.mapping_el["<mol_formula>"]`` dict, and if the
  mapping is found, apply it to all the rows of the original `states_electronic_raw`
  dataframe.
* If custom rules do not exist, apply the default parsing function.
* Check if all the new values under the State column are pyvalem-parseable.
* If some are not, raise the DatasetPostProcessorError with the original state strings
  which could not be parsed, prompting to implement custom rules into the
  ``input/mapping_el.py`` input file.

The existence of both `states_electronic_raw.csv` and `states_electronic.csv`
indicates that the post-processing happened already before, in which case an exception
is raised.
"""

import re

import pandas as pd
from pyvalem.states import MolecularTermSymbol
from pyvalem.states.molecular_term_symbol import MolecularTermSymbolError
from tqdm import tqdm

from config.config import OUTPUT_DIR
from .exceptions import DatasetPostProcessorError, CouldNotParseState


class DatasetPostProcessor:
    """Class post-processing the ``exomol2lida`` outputs into LIDA-ready files,
    expected by the LIDA's database populating logic.

    In general, this layer is needed because the ExoMol data files do not use
    pyvalem-compatible State strings. Lida database expects all state strings to
    be compatible with ``pyvalem.states.``

    Parameters
    ----------
    mol_formula : str
        Molecule formula string, as existing as a key in the input.molecules.

    Attributes
    ----------
    mol_formula : str
    output_dir : Path
    states_electronic_raw_path : Path
    states_electronic_raw : pd.DataFrame or None
    states_electronic_path : Path
    states_electronic : pd.DataFrame or None

    Raises
    ------
    DatasetPostProcessorError
        Under several circumstances:
        * If the formula is not among the outputs in the ``OUTPUT_DIR``.
        * If the formula has been post-processed already.
        * If the some of the states could not be converted into pyvalem-compatible
          strings.
    """

    # default electronic state pattern
    label_pattern = r"(?P<label>[a-zA-Z])"
    prime_pattern = r"(?P<prime>p)?"
    spin_pattern = r"(?P<spin>\d)"
    lambdas = ["Sigma", "Pi", "Delta", "Phi"]
    lambdas.extend([l.upper() for l in lambdas])
    lambda_pattern = rf"(?P<lambda>{'|'.join(l for l in lambdas)})"
    sym_pattern = r"(?P<sym>\+|\-)?"
    state_el_pattern = re.compile(
        f"^{label_pattern}{prime_pattern}{spin_pattern}{lambda_pattern}{sym_pattern}$"
    )

    def __init__(self, mol_formula):
        # first, verify that the dataset is among outputs and has not been processed
        # already:
        self.mol_formula = mol_formula
        self.output_dir = OUTPUT_DIR / mol_formula
        if not self.output_dir.joinpath("meta_data.json").is_file():
            raise DatasetPostProcessorError(
                f"The {mol_formula} data are not among the outputs in "
                f"{self.output_dir}!"
            )

        self.states_electronic_raw_path = self.output_dir / "states_electronic_raw.csv"
        if not self.states_electronic_raw_path.is_file():
            self.states_electronic_raw = None
        else:
            self.states_electronic_raw = pd.read_csv(
                self.states_electronic_raw_path, index_col=0, header=0
            )
        self.states_electronic_path = self.output_dir / "states_electronic.csv"
        if (
            self.states_electronic_raw is not None
            and self.states_electronic_path.is_file()
        ):
            raise DatasetPostProcessorError(
                f"The {mol_formula} output has already been post-processed, as "
                f"{self.states_electronic_path} exists!"
            )
        self.states_electronic = None

    @staticmethod
    def _parse_state_default(state_strs):
        """Electronic state string parsing function in the default case.

        This is the the default function called as a to parse the values
        of a single row of ``states_electronic.csv`` (excluding the index).
        It is intended solely for the case of a single column with values like
        "X1Sigma+", "a3Pi", "Ap1Pi", etc. These values do appear quite often in the
        .states files and they can be parsed into pyvalem-compatible states_str without
        specifying any custom parsing rules.
        If the parsing cannot be done (either if ``len(state_strs) != 1`` or if the
        structure is not like the examples shown above), `CouldNotParseState` is raised.

        Parameters
        ----------
        state_strs : iterable of str
            Values of a single row of the `states_electronic.csv` table, excluding the
            index. In most cases, this will be just a single value, but multiple
            values are theoretically supported (via the custom parsing rules), so the
            parameter takes the type of iterable[str].

        Returns
        -------
        str
            (Hopefully) PyValem-compatible state string (MolecularTermSymbol)

        Raises
        ------
        CouldNotParseState
            If the default parsing function does not apply for the `state_strs` passed.

        Examples
        --------
        >>> DatasetPostProcessor._parse_state_default(["g2Pi"])
        'g(2PI)'
        >>> DatasetPostProcessor._parse_state_default(["X4Sigma+"])
        'X(4SIGMA+)'
        >>> DatasetPostProcessor._parse_state_default(["Ap4Phi"])
        "A'(4PHI)"
        >>> DatasetPostProcessor._parse_state_default(["X", "4Sigma+"])
        Traceback (most recent call last):
          ...
        exomol2lida.exceptions.CouldNotParseState
        >>> DatasetPostProcessor._parse_state_default(["foo"])
        Traceback (most recent call last):
          ...
        exomol2lida.exceptions.CouldNotParseState
        """
        if len(state_strs) != 1:
            raise CouldNotParseState
        state_str = state_strs[0]
        match = DatasetPostProcessor.state_el_pattern.match(state_str)
        if match is None:
            raise CouldNotParseState
        label = match.group("label")
        if match.group("prime") is not None:
            label += "'"
        sym = "" if match.group("sym") is None else match.group("sym")
        return f"{label}({match.group('spin')}{match.group('lambda').upper()}{sym})"

    def _log_states_metadata(self):
        """Log the post-processed electronic states into the relevant table.

        The states are saved into the ``outputs/{mol_formula}/states_electronic.csv``,
        with two columns: state index ``i`` and state ``State``.

        Must be called only after `self.states_electronic` is populated.
        """
        assert self.states_electronic is not None, "Defense, should never happen"
        self.states_electronic.to_csv(self.states_electronic_path)

    def postprocess(self):
        """The main method for post-processing electronic states.

        Populates the self.states_electronic attribute with pd.DataFrame with a single
        column ``"State"`` and with values consisting of `pyvalem`-valid
        molecular term symbols.

        The dataframe gets also logged as ``states_electronic.csv`` into the relevant
        output directory.

        If the dataset under ``self.mol_formula`` does not resolve electronic states
        (does not include ``states_electronic_raw.csv``), this method does nothing.

        Raises
        ------
        DatasetPostProcessorError
            If the raw electronic states could no be parsed into pyvalem-compatible
            state strings by the default parsing method. This can be overcame by
            adding the pyvalem state strings into the `special_cases` ``dict``, as
            prompted by the error message.
        """
        if self.states_electronic_raw is None:
            return
        raw_states = [
            tuple(val) for val in self.states_electronic_raw.drop_duplicates().values
        ]
        map_raw_to_valid = {}
        failed_to_parse = []
        for raw_state in raw_states:
            # look into the special cases table:
            from input.mapping_el import mapping_el as special_cases

            valid_state = (
                special_cases.get(self.mol_formula, {}).get(raw_state, None)
            )
            try:
                if valid_state is None:
                    valid_state = self._parse_state_default(raw_state)
                MolecularTermSymbol(valid_state)
                map_raw_to_valid[raw_state] = valid_state
            except (CouldNotParseState, MolecularTermSymbolError):
                failed_to_parse.append(raw_state)
        if failed_to_parse:
            raise DatasetPostProcessorError(
                f"Add pyvalem-valid MolecularTermSymbol strings into "
                f"input.mapping_el.mapping_el['{self.mol_formula}'] under the "
                f"following keys: {str(failed_to_parse)[1:-1]}."
            )
        # now I have pyvalem-valid molecular term symbols, so just re-build the table
        self.states_electronic = pd.DataFrame(
            index=self.states_electronic_raw.index, columns=["State"]
        )
        for i in tqdm(
            self.states_electronic_raw.index, desc=f"{self.mol_formula} postprocessing"
        ):
            self.states_electronic.at[i, "State"] = map_raw_to_valid[
                tuple(self.states_electronic_raw.loc[i])
            ]
        # log the results and fuck off...
        self._log_states_metadata()


def postprocess_molecule(mol_formula, raise_exceptions=True):
    """A top-level function for post-processing exomol dataset which has already been
    processed (lumped and all).

    See the `DatasetPostProcessor` class for further documentation on errors etc.

    Parameters
    ----------
    mol_formula : str
        Molecular formula, must be among the subdirectories of ``config.OUTPUT_DIR``.
    raise_exceptions : bool
        If False, any exceptions raised by the `DataPostProcessor` constructor or its
        `postprocess` method will be caught and printed to stdout, instead of halting
        the program.

    Raises
    ------
    DatasetPostProcessorError
    """
    if raise_exceptions:
        mol_postprocessor = DatasetPostProcessor(mol_formula)
        mol_postprocessor.postprocess()
    else:
        try:
            mol_postprocessor = DatasetPostProcessor(mol_formula)
            mol_postprocessor.postprocess()
        except (DatasetPostProcessorError,) as e:
            print(f"{mol_formula}: POST-PROCESSING ABORTED: {type(e).__name__}: {e}")
    print()
