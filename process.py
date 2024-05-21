import sys
from functools import partial

from exomol2lida.process_dataset import process_molecule

if __name__ == "__main__":
    mol_formula = sys.argv[1]
    args = sys.argv[2:]
    allowed_args = {"--include-tau", "--postprocess"}
    assert set(args).issubset(allowed_args)

    proc_mol = partial(
        process_molecule,
        include_original_lifetimes=("--include-lifetimes" in args),
        postprocess=("--postprocess" in args),
        raise_exceptions=False,
    )

    if mol_formula.lower() == "all":
        from input.molecules import molecules as mol_formulas

        for mf in mol_formulas:
            proc_mol(mf)
    else:
        proc_mol(mol_formula)
