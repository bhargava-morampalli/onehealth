"""Recognise a host by the names NCBI Taxonomy gives to animals and to people.

The lexicon's word lists name the hosts seen most often. This covers the rest: any
chordate by its scientific name at any rank ("Larus dominicanus", "Macaca", "Cervidae")
and any animal by its common name ("gull", "sika deer"). The names are in
data/animal_names.tsv.gz, made by tools/build_animal_names.py.

A value is recognised only as a whole name, never by a word inside it: "sea water" is
not the fish genus Sea, and "colon cancer" is neither of two animal genera.
"""

from __future__ import annotations

import gzip
from functools import lru_cache
from importlib.resources import files

from onehealth.text import tokenize

ANIMAL = "animal"
HUMAN = "human"
KINDS = (ANIMAL, HUMAN)

# After a genus: "Larus sp.", "Corvus spp.". After a species: a subspecies or breed.
_UNNAMED_SPECIES = frozenset({"sp", "spp", "species", "ssp", "subsp"})
_HUMAN_NAMES = frozenset({"human", "humans"})


@lru_cache(maxsize=1)
def _names() -> tuple[frozenset[str], frozenset[str]]:
    """(scientific names of chordates, common names of animals). Read once, on first use."""
    scientific: set[str] = set()
    common: set[str] = set()
    path = files("onehealth").joinpath("data/animal_names.tsv.gz")
    with path.open("rb") as raw, gzip.open(raw, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            kind, _, name = line.rstrip("\n").partition("\t")
            (scientific if kind == "s" else common).add(name)
    return frozenset(scientific), frozenset(common)


@lru_cache(maxsize=65536)
def taxon_kind(value: str | None) -> str:
    """"human", "animal", or "" when the value is not a name of either.

    People: the genus Homo ("Homo sapiens", "Homo sapiens sapiens", "Homo") and "human".
    Animals, in this order:
      * a common name, singular or plural ("gull" for NCBI's "gulls");
      * a scientific name of a chordate ("Macaca mulatta", "Macaca", "Aves");
      * a species followed by a subspecies or breed ("Sus scrofa domesticus"), or a
        genus followed by sp. ("Meleagris sp.").
    """
    words = tokenize(value)
    while words and words[-1] in _UNNAMED_SPECIES:
        words.pop()
    if not words:
        return ""
    name = " ".join(words)
    if words[0] == "homo" or name in _HUMAN_NAMES:
        return HUMAN
    scientific, common = _names()
    other_number = name[:-1] if name.endswith("s") else name + "s"
    if name in common or other_number in common or name in scientific:
        return ANIMAL
    if len(words) >= 3 and " ".join(words[:2]) in scientific:
        return ANIMAL
    return ""
