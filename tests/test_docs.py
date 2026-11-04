"""The documentation site (documentation/) must keep up with the code.

Each check lists what the code has and the page lacks, so a new option, column or
lexicon rule cannot ship undocumented.
"""

import argparse
import re
from pathlib import Path

import pytest
import yaml

from onehealth.cli import _build_parser
from onehealth.constants import (
    CLASSIFY_EXTRA_FIELDS,
    RESOLVE_FIELDS,
    SEED_LEXICON_VERSION,
    UNKNOWN_EXACT,
)
from onehealth.lexicon import shipped_lexicon_path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "documentation"

PAGES = {
    "resolve": "commands/resolve.md",
    "classify": "commands/classify.md",
    "summarize": "commands/summarize-plot.md",
    "plot": "commands/summarize-plot.md",
    "validate": "commands/validate.md",
    "report": "commands/report.md",
}


def _page(name: str) -> str:
    return (DOCS / name).read_text(encoding="utf-8")


def _headings(text: str):
    """(line number, level, title) of each markdown heading outside code blocks."""
    fenced = False
    for number, line in enumerate(text.splitlines()):
        if line.lstrip().startswith("```"):
            fenced = not fenced
        match = None if fenced else re.match(r"(#+) (.+)", line)
        if match:
            yield number, len(match.group(1)), match.group(2).strip()


def _section(text: str, heading: str) -> str:
    """The text under a markdown heading, up to the next heading of the same or a higher level."""
    lines = text.splitlines()
    headings = list(_headings(text))
    for i, (start, level, title) in enumerate(headings):
        if title == heading:
            end = next((n for n, lv, _ in headings[i + 1 :] if lv <= level), len(lines))
            return "\n".join(lines[start + 1 : end])
    raise AssertionError(f"no heading {heading!r}")


def _commands() -> dict[str, argparse.ArgumentParser]:
    parser = _build_parser()
    action = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))
    return dict(action.choices)


def test_every_command_has_a_page():
    assert set(_commands()) == set(PAGES)
    for page in set(PAGES.values()):
        assert (DOCS / page).is_file(), page


@pytest.mark.parametrize("command", sorted(PAGES))
def test_every_option_is_on_its_page(command):
    parser = _commands()[command]
    text = _page(PAGES[command])
    if PAGES[command] == PAGES["summarize"]:
        text = _section(text, command)  # summarize and plot share a page; each has its own options
    options = [s for a in parser._actions for s in a.option_strings if s.startswith("--") and s != "--help"]
    assert options, command
    missing = [option for option in options if f"`{option}`" not in text]
    assert not missing, f"{PAGES[command]} does not mention {missing}"


def test_output_columns_are_documented():
    resolve = _page("commands/resolve.md")
    assert not [f for f in RESOLVE_FIELDS if f"`{f}`" not in resolve]
    classify = _section(_page("commands/classify.md"), "Output")
    assert not [f for f in CLASSIFY_EXTRA_FIELDS if f"`{f}`" not in classify]


def test_every_missing_value_word_is_listed():
    text = _section(_page("guides/input-files.md"), "Missing values")
    missing = [word for word in sorted(UNKNOWN_EXACT) if word and f"`{word}`" not in text]
    assert not missing, f"guides/input-files.md does not list {missing}"


def test_every_rule_of_the_seed_lexicon_is_in_the_guide():
    text = _section(_page("guides/lexicon.md"), f"The rules in v{SEED_LEXICON_VERSION}")
    rules = yaml.safe_load(shipped_lexicon_path(SEED_LEXICON_VERSION).read_text(encoding="utf-8"))["rules"]
    rows = {m.group(2): int(m.group(1)) for m in re.finditer(r"^\| (\d+) \| `([^`]+)` \|", text, re.MULTILINE)}
    assert rows == {rule["id"]: int(rule.get("priority") or 0) for rule in rules}


def test_every_shipped_lexicon_version_is_listed():
    text = _section(_page("guides/lexicon.md"), "Versions")
    shipped = sorted(p.stem for p in shipped_lexicon_path(SEED_LEXICON_VERSION).parent.glob("v*.yaml"))
    assert not [v for v in shipped if f"| {v} |" not in text]


def test_every_page_is_in_the_navigation():
    nav = (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
    pages = sorted(p.relative_to(DOCS).as_posix() for p in DOCS.rglob("*.md"))
    assert not [p for p in pages if f": {p}" not in nav]
