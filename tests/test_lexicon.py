from pathlib import Path

import pytest

from onehealth.lexicon import load_lexicon
from onehealth.text import is_taxon_only_source

# Host binomials are meant to be caught as taxon-only when typed into isolation_source.
_BINOMIALS = {
    "gallus gallus",
    "bos taurus",
    "sus scrofa",
    "ovis aries",
    "mus musculus",
    "meleagris gallopavo",
    "anas platyrhynchos",
    "oryctolagus cuniculus",
    "homo sapiens",
    "homo sapien",
}


def test_seed_lexicon_phrases_are_not_mistaken_for_taxon_names(lexicon):
    """The taxon guard runs before matching; a new 2-3 word phrase must not trip it."""
    flagged = []
    for rule in lexicon.rules:
        for words in rule.phrases + rule.exact:
            phrase = " ".join(words)
            if phrase not in _BINOMIALS and is_taxon_only_source(phrase):
                flagged.append(f"{rule.rule_id}: {phrase}")
    assert flagged == []


def test_unknown_rule_field_is_rejected(tmp_path: Path):
    path = tmp_path / "lex.yaml"
    path.write_text(
        'version: "x"\nrules:\n  - id: r\n    fields: ["isolation source"]\n    phrases: [blood]\n    l1: Human\n',
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unknown field"):
        load_lexicon(path)


def test_seed_lexicon_version_matches_file(lexicon):
    assert lexicon.version == "2"
