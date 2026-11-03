"""Normalization helpers for isolation-source matching."""

from __future__ import annotations

import re

from onehealth.constants import UNKNOWN_EXACT, UNKNOWN_GENERIC_TAILS

_SPLIT_RE = re.compile(r"[^a-z0-9]+")


def normalize_text(value: str | None) -> str:
    text = (value or "").strip().lower()
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text)
    return text


def tokenize(value: str | None) -> list[str]:
    normalized = normalize_text(value)
    if not normalized:
        return []
    return [part for part in _SPLIT_RE.split(normalized) if part]


def contains_phrase(haystack: list[str], phrase: list[str] | tuple[str, ...]) -> bool:
    needle = list(phrase)
    if not needle or len(needle) > len(haystack):
        return False
    n = len(needle)
    for i in range(len(haystack) - n + 1):
        if haystack[i : i + n] == needle:
            return True
    return False


def is_blank(value: str | None) -> bool:
    return normalize_text(value) == ""


def is_unknown_value(value: str | None) -> bool:
    """True when a field is empty, an unknown marker, or 'unknown <generic>'."""
    normalized = normalize_text(value)
    if normalized in UNKNOWN_EXACT:
        return True
    tokens = tokenize(normalized)
    if not tokens:
        return True
    # INSDC missing-value terms come as "missing" or "missing: <reason>".
    if tokens[0] in {"unknown", "unk", "missing"}:
        return True
    if tokens[-1] in {"unknown", "unk"} and set(tokens[:-1]).issubset(UNKNOWN_GENERIC_TAILS):
        return True
    return False


# Common pathogen abbreviations stuffed into isolation_source instead of a specimen.
_TAXON_ABBREVIATIONS = frozenset(
    {
        ("e", "coli"),
        ("s", "aureus"),
        ("s", "enterica"),
        ("s", "typhimurium"),
        ("k", "pneumoniae"),
        ("k", "oxytoca"),
        ("p", "aeruginosa"),
        ("c", "difficile"),
        ("l", "monocytogenes"),
        ("c", "jejuni"),
        ("e", "faecalis"),
        ("e", "faecium"),
        ("a", "baumannii"),
        ("e", "cloacae"),
        ("c", "freundii"),
        ("p", "mirabilis"),
    }
)

# English isolation-source tokens. If any appear, the field is not taxon-only.
_NON_TAXON_TOKENS = frozenset(
    {
        "blood",
        "serum",
        "plasma",
        "urine",
        "sputum",
        "stool",
        "feces",
        "faeces",
        "csf",
        "cerebrospinal",
        "wound",
        "pus",
        "abscess",
        "swab",
        "tissue",
        "biopsy",
        "bile",
        "catheter",
        "tracheal",
        "bronchial",
        "bronchoalveolar",
        "lavage",
        "bal",
        "peritoneal",
        "pleural",
        "hospital",
        "clinical",
        "healthcare",
        "ward",
        "icu",
        "environment",
        "environmental",
        "sink",
        "drain",
        "drainage",
        "faucet",
        "tap",
        "water",
        "river",
        "lake",
        "soil",
        "sewage",
        "wastewater",
        "sludge",
        "sediment",
        "air",
        "compost",
        "irrigation",
        "groundwater",
        "seawater",
        "marine",
        "food",
        "meat",
        "beef",
        "pork",
        "chicken",
        "poultry",
        "carcass",
        "milk",
        "dairy",
        "cheese",
        "lettuce",
        "salad",
        "flour",
        "spice",
        "produce",
        "vegetable",
        "fruit",
        "seafood",
        "retail",
        "animal",
        "livestock",
        "companion",
        "farm",
        "barn",
        "manure",
        "litter",
        "feed",
        "laboratory",
        "lab",
        "strain",
        "culture",
        "vitro",
        "reference",
        "type",
        "patient",
        "human",
        "carriage",
        "fresh",
        "freshwater",
        "surface",
        "well",
        "ground",
        "minced",
        "sample",
        "isolate",
        "source",
        "unknown",
        "missing",
        "cloaca",
        "cloacal",
        "dog",
        "cat",
        "pig",
        "cow",
        "calf",
        "sheep",
        "goat",
        "horse",
        "mouse",
        "mice",
        "rat",
        "turkey",
        "duck",
        "fish",
        "salmon",
        "bird",
        "cattle",
        "bovine",
        "swine",
        "porcine",
        "canine",
        "feline",
        "equine",
        "ovine",
        "caprine",
        "avian",
        "wildlife",
        "producing",
        "product",
        "surgical",
        "abdominal",
        "chest",
        "experimentally",
        "experimental",
        "experiment",
        "derived",
        "inbred",
        "knockout",
        "mutant",
        "alpha",
        "beta",
        "positive",
        "negative",
        "hemolytic",
        "resistant",
        "control",
        "blank",
        "pool",
        "pooled",
        "mixed",
        "from",
        "ventilator",
        "tubing",
        "circuit",
        "indoor",
        "outdoor",
        "ambient",
        "sea",
        "breast",
        "products",
        "specimen",
        "nasal",
        "rectal",
        "throat",
        "gut",
        "asymptomatic",
        "shower",
        "floor",
        "bathroom",
        "biofilm",
        "trap",
        "setting",
        "env",
        "collection",
        # Real sources once mistaken for binomials (IMP benchmark, NCBI Pathogen Detection).
        "influent",
        "effluent",
        "wwtp",
        "bedpan",
        "washer",
    }
)

# English word endings that no Latin binomial uses ("cleaning", "infection").
_ENGLISH_SUFFIXES = ("ing", "tion", "sion", "ment", "ness")


def is_taxon_only_source(value: str | None) -> bool:
    """True when isolation_source is only a scientific name, not a specimen.

    GenFS cleanup: taxonomy belongs in host. Do not treat `Escherichia coli`
    or `Bos taurus` as an isolation site. This is a word-shape guess (two or
    three plain words), so classify lets a matching lexicon rule overrule it.
    """
    tokens = tokenize(value)
    if not tokens:
        return False
    if any(not token.isalpha() for token in tokens):
        return False
    if tuple(tokens) in _TAXON_ABBREVIATIONS:
        return True
    if any(token in _NON_TAXON_TOKENS or token.endswith(_ENGLISH_SUFFIXES) for token in tokens):
        return False
    if len(tokens) not in {2, 3}:
        return False
    return all(len(token) >= 3 for token in tokens)
