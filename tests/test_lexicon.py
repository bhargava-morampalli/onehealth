from pathlib import Path

import pytest

from onehealth.constants import SEED_LEXICON_VERSION
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


def test_host_rule_phrases_are_not_mistaken_for_taxon_names(lexicon):
    """In a rule that reads `host`, a phrase that looks like a binomial is ignored in isolation_source.

    Rules for source text only are exempt: their matches overrule the taxon guess.
    """
    flagged = []
    for rule in lexicon.rules:
        if "host" not in rule.fields:
            continue
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
    assert lexicon.version == SEED_LEXICON_VERSION == "6"
