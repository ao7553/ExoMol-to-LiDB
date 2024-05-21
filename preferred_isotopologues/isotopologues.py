"""
A module with functionality to determine a "preferred isotopologue" for any molecular
formula.
"""

from pyvalem.formula import Formula

from .atomic_isotope_abundances import get_nist_isotopes_abundances


def get_preferred_isotopologue(molecule_formula):
    """Gives formula of the most abundand isotopologue.

    Decomposes the passed molecule formula into elements, finds the most abundant
    isotopes for each element, and reconstructs the formula back with the isotopes
    instead of the plain elements.

    Parameters
    ----------
    molecule_formula : str
        Molecular formula without isotope information, e.g. 'H2O', or 'C+'.

    Returns
    -------
    str
        The formula of what is assumed to be the most abundad isotopologue
        corresponding to the passed molecular formula. E.g. '(1H)2(16O)', or '(12C)+'.

    Examples
    --------
    >>> get_preferred_isotopologue('H2O')
    '(1H)2(16O)'
    >>> get_preferred_isotopologue('C+')
    '(12C)+'
    >>> get_preferred_isotopologue('trans-P2H2')
    'trans-(31P)2(1H)2'
    """

    formula = Formula(molecule_formula)
    iso_formula = ""
    if molecule_formula.startswith("cis-"):
        iso_formula += "cis-"
    elif molecule_formula.startswith("trans-"):
        iso_formula += "trans-"
    for atom, stoich in formula.atom_stoich.items():
        iso_formula += get_nist_isotopes_abundances()[atom][0][0]
        if stoich > 1:
            iso_formula += str(stoich)
    if formula.charge:
        assert formula.charge == 1
        iso_formula += "+"
    # validation of the name:
    Formula(iso_formula)

    return iso_formula
