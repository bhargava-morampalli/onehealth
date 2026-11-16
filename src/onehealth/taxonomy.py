"""Recognise a host by the names NCBI Taxonomy gives to animals and to people.

The lexicon's word lists name the hosts seen most often. This covers the rest: any
chordate by its scientific name at any rank ("Larus dominicanus", "Macaca", "Cervidae"),
a mollusc, a crustacean or another animal with a common name by its species name
("Magallana gigas", "Musca domestica"), and any animal by its common name ("gull",
"sika deer"). The names are in data/animal_names.tsv.gz, made by
tools/build_animal_names.py.

taxon_kind reads a host value. The value must be a name, or end with one ("European
rook", "free-range white-tailed deer"): a word in the middle is not enough, so "sea
water" is not the fish genus Sea, and "colon cancer" is neither of two animal genera.

animal_name_within reads a text that is not a host value, such as the words around a
specimen ("stool (Odocoileus hemionus)", "pigeon stool"). It is stricter about single
words, because ordinary words are also animals' names.
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
# After the name in a host value: the animal's age. "white stork nestling" is a white stork.
_AGE_WORDS = frozenset({"nestling", "nestlings", "fledgling", "fledglings", "hatchling", "hatchlings",
                        "juvenile", "juveniles", "adult", "adults", "chick", "chicks"})
# A host value with one of these lists several things: "bologna, hot dogs, sausage-like products".
_LIST_MARKS = (",", ";", "/", "&", " and ", " or ")
# Vertebrates whose one-word English name is also an ordinary word in a sample description:
# the sole of a foot, an X-ray, a cardiac monitor, a bed rail.
_WORDS_NOT_NAMES = frozenset({"sole", "ray", "monitor", "crane", "rail", "swift", "cardinal", "bass", "char",
                              "kid", "permit", "ide", "ling", "dab", "tang", "grunt", "drum",
                              # "placentals" is NCBI's English name for a clade of mammals.
                              "placental", "placentals"})
_LONGEST_NAME = 4  # words


@lru_cache(maxsize=1)
def _tables() -> tuple[frozenset[str], frozenset[str], frozenset[str]]:
    """(scientific names, common names of chordates, common names of other animals). Read once, on first use."""
    tables: dict[str, set[str]] = {"s": set(), "v": set(), "c": set()}
    path = files("onehealth").joinpath("data/animal_names.tsv.gz")
    with path.open("rb") as raw, gzip.open(raw, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            kind, _, name = line.rstrip("\n").partition("\t")
            tables[kind].add(name)
    return frozenset(tables["s"]), frozenset(tables["v"]), frozenset(tables["c"])


def _names() -> tuple[frozenset[str], frozenset[str]]:
    """(scientific names, common names of any animal)."""
    scientific, vertebrate, other = _tables()
    return scientific, vertebrate | other


def _in_either_number(name: str, names: frozenset[str]) -> bool:
    """NCBI gives some common names in the plural ("gulls"); a host value may use either."""
    return name in names or (name[:-1] if name.endswith("s") else name + "s") in names


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
    if _is_animal_name(words):
        return ANIMAL
    # A name at the end of a longer value: "European rook", "forest musk deer", "others -
    # Daboia russelii". English puts the noun last, so the words before it describe the
    # animal. Not in a list of several things, where the last one need not be the host.
    if any(mark in (value or "").lower() for mark in _LIST_MARKS):
        return ""
    while len(words) > 1 and words[-1] in _AGE_WORDS:
        words.pop()
    vertebrate = _tables()[1]
    for start in range(1, len(words)):
        tail = words[start:]
        if len(tail) > _LONGEST_NAME:
            continue
        # One word at the end must name a vertebrate: a "kitchen sponge" is not an animal.
        if _is_animal_name(tail) if len(tail) > 1 else _in_either_number(tail[0], vertebrate):
            return ANIMAL
    return ""


def _is_animal_name(words: list[str]) -> bool:
    """The words, taken together, are a name NCBI Taxonomy gives to an animal."""
    scientific, common = _names()
    name = " ".join(words)
    if _in_either_number(name, common) or name in scientific:
        return True
    return len(words) >= 3 and " ".join(words[:2]) in scientific


@lru_cache(maxsize=65536)
def animal_name_within(text: str | None) -> str:
    """The first animal's name inside a text that is not a host value, or "".

    For the words around a specimen: "stool (Odocoileus hemionus)", "giant panda feces",
    "pigeon stool". A name of two or more words counts whatever the animal. A single word
    counts only when it is the common name of a vertebrate ("pigeon", "deer", "beaver")
    and not also an ordinary word ("sole", "monitor"). The English names of moths,
    butterflies and other invertebrates are ordinary words too often ("white", "tip",
    "sponge", "tick"), so one of those alone is not a name here.
    """
    words = tokenize(text)
    scientific, vertebrate, other = _tables()
    for length in range(min(_LONGEST_NAME, len(words)), 1, -1):
        for start in range(len(words) - length + 1):
            name = " ".join(words[start : start + length])
            if name in scientific or _in_either_number(name, vertebrate) or _in_either_number(name, other):
                return name
    for word in words:
        if word not in _WORDS_NOT_NAMES and len(word) > 2 and _in_either_number(word, vertebrate):
            return word
    return ""
