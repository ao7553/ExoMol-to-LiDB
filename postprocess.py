import sys

from exomol2lida.postprocess_dataset import postprocess_molecule

if __name__ == "__main__":
    mol_formula = sys.argv[1]

    if mol_formula.lower() == "all":
        from input.molecules import molecules as mol_formulas

        for mf in mol_formulas:
            postprocess_molecule(mf, raise_exceptions=False)
    else:
        postprocess_molecule(mol_formula, raise_exceptions=False)
