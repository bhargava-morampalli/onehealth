"""What a release to PyPI depends on (.github/workflows/publish.yml)."""

import importlib.util
import re
from pathlib import Path

import onehealth

ROOT = Path(__file__).resolve().parent.parent


def _release_notes():
    spec = importlib.util.spec_from_file_location("release_notes", ROOT / "tools" / "release_notes.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _version() -> str:
    """The [project] version in pyproject.toml (tomllib is Python 3.11+, and the tests run on 3.10)."""
    found = re.search(r'^version = "([^"]+)"$', (ROOT / "pyproject.toml").read_text(encoding="utf-8"), re.MULTILINE)
    assert found, "pyproject.toml names no version"
    return found.group(1)


def test_the_package_names_one_version():
    assert onehealth.__version__ == _version()


def test_the_version_has_its_changelog_section_for_the_release_notes():
    notes = _release_notes().section((ROOT / "CHANGELOG.md").read_text(encoding="utf-8"), _version())
    assert len(notes) > 200 and not notes.startswith("## ")


def test_release_notes_are_the_versions_section_only():
    text = "# Changelog\n\nIntro.\n\n## 2 May: 0.5.0, later\n\nFive.\n\n## 1 May: 0.4.0, the first\n\nFour.\n\n## Older\n\nOld.\n"
    notes = _release_notes()
    assert notes.section(text, "0.4.0") == "Four.\n"
    assert notes.section(text, "0.5.0") == "Five.\n"
    try:
        notes.section(text, "0.4")  # not a prefix of another version
    except SystemExit as stop:
        assert "no section for version 0.4" in str(stop)
    else:
        raise AssertionError("a version with no section must stop the release")


def test_the_readme_has_no_relative_links():
    """PyPI shows README.md as the project page, where a link into the repository leads nowhere."""
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    relative = [target for target in re.findall(r"\]\(([^)\s]+)\)", readme) if not target.startswith(("https://", "http://", "#"))]
    assert relative == []


def test_the_install_instructions_name_the_current_version():
    for page in ("README.md", "documentation/install.md"):
        text = (ROOT / page).read_text(encoding="utf-8")
        assert "pip install onehealth\n" in text, page
        assert f"pip install onehealth=={_version()}" in text, page
