"""Print the CHANGELOG.md section of one version, for its GitHub Release.

    python tools/release_notes.py 0.4.0

The section is the one whose heading names the version ("## 17 November 2026: 0.4.0, ..."),
up to the next heading of the same level. A version with no section stops with an error,
so a release cannot go out without its changelog entry. Used by .github/workflows/publish.yml.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

CHANGELOG = Path(__file__).resolve().parent.parent / "CHANGELOG.md"


def section(text: str, version: str) -> str:
    """The body of the version's section, without its heading line."""
    lines = text.splitlines()
    named = re.compile(rf"^## .*[:\s]{re.escape(version)}[,\s]")
    for start, line in enumerate(lines):
        if named.match(line + " "):
            end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")), len(lines))
            return "\n".join(lines[start + 1 : end]).strip() + "\n"
    raise SystemExit(f"CHANGELOG.md has no section for version {version}: add one before releasing it")


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 1:
        raise SystemExit("usage: python tools/release_notes.py VERSION")
    sys.stdout.write(section(CHANGELOG.read_text(encoding="utf-8"), args[0]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
