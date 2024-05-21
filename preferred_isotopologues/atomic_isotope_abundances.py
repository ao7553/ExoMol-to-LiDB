"""
A module reading the relative atomic isotope abundances from NIST
(or from a json cache).
"""

import json
from functools import lru_cache
from pathlib import Path

import requests
from lxml import html


@lru_cache()
def get_nist_isotopes_abundances(overwrite_cache=False):
    """Get a dictionary of atomic isotopes abundances.

    Parses the NIST page of isotope abundances from the url
    https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all
    and returns the contents in a dictionary keyed by element symbols.

    Parameters
    ----------
    overwrite_cache : bool
        If True, loads the result from the json cache. If False, or the cache does not
        exist, parses the NIST page.

    Returns
    -------
    dict[str, list[list[str, float]]]
        Keys are element symbols and values are lists of list pairs of isotope
        symbols and relative abundances.

    Examples
    --------
    >>> get_nist_isotopes_abundances(overwrite_cache=False)['C']
    [['(12C)', 0.9893], ['(13C)', 0.0107]]
    """

    file_path = Path(__file__)
    file_dir, file_name = file_path.parent, file_path.name
    cache_path = file_dir.joinpath(f"{file_name[:-3]}.json")

    if overwrite_cache or not cache_path.is_file():
        page = requests.get(
            "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all"
        )
        tree = html.fromstring(page.content)

        rows = tree.xpath("//table/tbody/tr")[2:]

        isotopes_abundances = {}

        element = None
        for row in rows:
            if element is not None and len(row.xpath("td")) == 1:
                if element in isotopes_abundances:
                    isotopes_abundances[element].sort(key=lambda x: -x[1])
                element = None
            tds = row.xpath("td")
            if element is None:
                if len(tds) == 1:
                    continue
                else:
                    element = tds[1].text.strip()
                    if element == "Og":
                        break
                    offset = 2
            else:
                if element == "H":
                    offset = 1
                else:
                    offset = 0
            isotope = f"({int(tds[offset].text.strip())}{element})"
            abundance = tds[offset + 2].text.strip().replace("\xa0", "").split("(")[0]
            if not abundance:
                continue
            abundance = float(abundance)
            try:
                isotopes_abundances[element].append([isotope, abundance])
            except KeyError:
                isotopes_abundances[element] = [
                    [isotope, abundance],
                ]

        # save the result as a json
        with open(cache_path, "w") as fp:
            json.dump(isotopes_abundances, fp, indent=2)
        return isotopes_abundances
    else:
        with open(cache_path, "r") as fp:
            isotope_abundances = json.load(fp)
        return isotope_abundances
