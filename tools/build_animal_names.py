"""Rebuild src/onehealth/data/animal_names.tsv.gz from NCBI Taxonomy (standard library only).

    python tools/build_animal_names.py                      # downloads new_taxdump.tar.gz (160 MB)
    python tools/build_animal_names.py new_taxdump.tar.gz   # or reads a copy you already have

The lexicon recognises a host by a list of words (`chicken`, `bos taurus`, ...). No list
covers every species a pathogen is isolated from, so a host the list does not know
(a gull, a macaque, a water buffalo) went to review. This table lets the lexicon's
`taxonomy: animal` rule recognise any of them: it holds the names NCBI Taxonomy gives
to animals.

Two kinds of name are kept, one per line as `kind<TAB>name`, lower-case, punctuation
removed:

* `s`  scientific names of chordates (vertebrates and their relatives) at every rank:
       `larus dominicanus`, `macaca`, `cervidae`, `aves`.
* `c`  common names of any animal: `gull`, `sika deer`, `rhesus macaque`, `mud crab`.

Why only chordates for scientific names: the matcher needs a whole name, and a genus on
its own is a name. Invertebrate genera include Sepsis (a fly), Cancer (a crab), Colon,
Trachea and Tibia, and a host field can say "sepsis" or "cancer". Chordate names were
checked against the words used in 159,045 NCBI Pathogen Detection records. Almost every
host of a pathogen isolate is a vertebrate; other animals are matched by common name.

Left out: people (the genus Homo; the lexicon's own rules and `taxonomy: human` handle
them), names NCBI gives to groups that are not organisms ("environmental samples",
"unclassified Aves"), and unnamed entries ("Rattus sp. ABC-2011").

NCBI updates the taxonomy every day, so a rebuild is not byte-for-byte the same as the
file in the repository. The header of the file records which dump it was made from.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import re
import sys
import tarfile
import urllib.request
from pathlib import Path

URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
OUT = Path(__file__).resolve().parent.parent / "src" / "onehealth" / "data" / "animal_names.tsv.gz"

# Words in names of groups that are not organisms, and of entries without a real name.
NOT_A_NAME = frozenset(
    {"unclassified", "environmental", "samples", "sample", "incertae", "sedis", "unidentified", "uncultured",
     "mixed", "synthetic", "sp", "cf", "aff", "x"}
)
COMMON = ("common name", "genbank common name")


def normalise(name: str) -> str:
    """Lower-case words separated by single spaces: the form onehealth.text.tokenize gives."""
    return " ".join(part for part in re.split(r"[^a-z0-9]+", name.lower().replace("_", " ")) if part)


def table(archive: tarfile.TarFile, member: str):
    """Rows of a .dmp file: fields are separated by <TAB>|<TAB> and each line ends with <TAB>|."""
    handle = archive.extractfile(member)
    assert handle is not None, member
    for line in io.TextIOWrapper(handle, encoding="utf-8"):
        yield [cell.strip() for cell in line.rstrip("\n").rstrip("|").split("\t|")]


def usable(name: str, max_words: int) -> bool:
    words = name.split()
    return bool(words) and len(words) <= max_words and not (set(words) & NOT_A_NAME) and not any(c.isdigit() for c in name)


def build(data: bytes) -> tuple[set[str], set[str]]:
    scientific: set[str] = set()
    common: set[str] = set()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        # tax_id, name, species, genus, family, order, class, phylum, kingdom, domain
        lineage = {row[0]: (row[1], row[3], row[7], row[8]) for row in table(archive, "rankedlineage.dmp")}
        for tax_id, name, _unique, name_class, *_ in table(archive, "names.dmp"):
            if tax_id not in lineage:
                continue
            taxon, genus, phylum, kingdom = lineage[tax_id]
            if kingdom != "Metazoa" or genus == "Homo" or taxon == "Homo":
                continue
            text = normalise(name)
            if name_class == "scientific name" and phylum == "Chordata" and usable(text, 3):
                scientific.add(text)
            elif name_class in COMMON and usable(text, 4):
                common.add(text)
    return scientific, common


def main(argv: list[str]) -> int:
    if len(argv) > 1:
        data, source = Path(argv[1]).read_bytes(), Path(argv[1]).name
    else:
        print(f"downloading {URL}", file=sys.stderr)
        with urllib.request.urlopen(URL, timeout=900) as response:
            data, source = response.read(), URL
    scientific, common = build(data)
    lines = [
        "# Names NCBI Taxonomy gives to animals. Made by tools/build_animal_names.py; do not edit by hand.",
        f"# Source: {source}, md5 {hashlib.md5(data).hexdigest()}",
        "# s = scientific name of a chordate at any rank; c = common name of any animal. People are not listed.",
    ]
    lines += [f"s\t{name}" for name in sorted(scientific)]
    lines += [f"c\t{name}" for name in sorted(common)]
    OUT.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode("utf-8"), compresslevel=9, mtime=0))
    print(f"wrote {OUT}: {len(scientific)} scientific and {len(common)} common names, "
          f"{OUT.stat().st_size / 1e6:.2f} MB", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
