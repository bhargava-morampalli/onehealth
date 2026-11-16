"""Rebuild src/onehealth/data/animal_names.tsv.gz from NCBI Taxonomy (standard library only).

    python tools/build_animal_names.py                      # downloads new_taxdump.tar.gz (160 MB)
    python tools/build_animal_names.py new_taxdump.tar.gz   # or reads a copy you already have
    python tools/build_animal_names.py --out names.tsv.gz   # writes somewhere else

To see what a newer dump would change before replacing the table, use
tools/refresh_animal_names.py, which also runs monthly in .github/workflows/taxonomy.yml.

The lexicon recognises a host by a list of words (`chicken`, `bos taurus`, ...). No list
covers every species a pathogen is isolated from, so a host the list does not know
(a gull, a macaque, a water buffalo) went to review. This table lets the lexicon's
`taxonomy: animal` rule recognise any of them: it holds the names NCBI Taxonomy gives
to animals.

Three kinds of name are kept, one per line as `kind<TAB>name`, lower-case, punctuation
removed:

* `s`  scientific names. Of chordates (vertebrates and their relatives) at every rank:
       `larus dominicanus`, `macaca`, `cervidae`, `aves`. Of other animals, only names
       of two or three words, which are species, and only of molluscs, of crustaceans,
       and of species NCBI also gives a common name (the ones people sample):
       `magallana gigas` (an oyster), `penaeus vannamei` (a shrimp), `musca domestica`.
       Synonyms of two or three words are kept too (`litopenaeus vannamei`). All
       quarter of a million insect species would make the table four times larger.
* `v`  common names of chordates: `gull`, `sika deer`, `rhesus macaque`.
* `c`  common names of other animals: `mud crab`, `house fly`.

Why a genus on its own only for chordates: the matcher needs a whole name, and
invertebrate genera include Sepsis (a fly), Cancer (a crab), Colon, Trachea and Tibia;
a host field can say "sepsis" or "cancer". A two-word name is safe. Chordate names were
checked against the words used in 159,045 NCBI Pathogen Detection records.

Why common names are split: the English names of moths and butterflies are ordinary
words ("the nun", "emperor", "white", "blue", "the uncertain"), and so are "sponge",
"tick" and "bug". As a whole host value they are still names. Inside a longer text
(onehealth.taxonomy.animal_name_within) a single word counts only when it names a
vertebrate.

Left out: people (the genus Homo; the lexicon's own rules and `taxonomy: human` handle
them), names NCBI gives to groups that are not organisms ("environmental samples",
"unclassified Aves"), and unnamed entries ("Rattus sp. ABC-2011").

NCBI updates the taxonomy every day, so a rebuild is not byte-for-byte the same as the
file in the repository. The header of the file records which dump it was made from: its
checksum, and the date NCBI wrote it (`# NCBI Taxonomy of 2026-11-14`), which `classify`
prints, so a methods section can name the names it used.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import re
import sys
import tarfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
OUT = Path(__file__).resolve().parent.parent / "src" / "onehealth" / "data" / "animal_names.tsv.gz"

# Words in names of groups that are not organisms, and of entries without a real name.
NOT_A_NAME = frozenset(
    {"unclassified", "environmental", "samples", "sample", "incertae", "sedis", "unidentified", "uncultured",
     "mixed", "synthetic", "sp", "cf", "aff", "x"}
)
COMMON = ("common name", "genbank common name")
# Classes of crustaceans (NCBI has no rank for Crustacea itself in rankedlineage.dmp).
CRUSTACEANS = frozenset({"Malacostraca", "Branchiopoda", "Hexanauplia", "Thecostraca", "Ostracoda", "Copepoda", "Maxillopoda"})
SYNONYM = ("synonym", "equivalent name", "genbank synonym")


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


def build(data: bytes) -> tuple[set[str], set[str], set[str]]:
    """(scientific names, common names of chordates, common names of other animals)."""
    scientific: set[str] = set()
    vertebrate: set[str] = set()
    common: set[str] = set()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        # tax_id, name, species, genus, family, order, class, phylum, kingdom, domain
        lineage = {row[0]: (row[1], row[3], row[6], row[7])
                   for row in table(archive, "rankedlineage.dmp")
                   if row[8] == "Metazoa" and row[3] != "Homo" and row[1] != "Homo"}
        names = [row[:4] for row in table(archive, "names.dmp") if row[0] in lineage]
        has_common_name = {tax_id for tax_id, _name, _unique, name_class in names if name_class in COMMON}
        for tax_id, name, _unique, name_class in names:
            _taxon, _genus, taxon_class, phylum = lineage[tax_id]
            text = normalise(name)
            words = len(text.split())
            chordate = phylum == "Chordata"
            sampled = chordate or phylum == "Mollusca" or taxon_class in CRUSTACEANS or tax_id in has_common_name
            if name_class == "scientific name" and usable(text, 3) and (chordate or (words >= 2 and sampled)):
                scientific.add(text)
            elif name_class in SYNONYM and usable(text, 3) and words >= 2 and sampled:
                scientific.add(text)
            elif name_class in COMMON and usable(text, 4):
                (vertebrate if chordate else common).add(text)
    return scientific, vertebrate, common - vertebrate


def read_dump(path: str | None) -> tuple[bytes, str]:
    """The dump's bytes and what to call it in the header: a file name, or the URL it came from."""
    if path:
        return Path(path).read_bytes(), Path(path).name
    print(f"downloading {URL}", file=sys.stderr)
    with urllib.request.urlopen(URL, timeout=900) as response:
        return response.read(), URL


def dump_date(data: bytes) -> str:
    """The day NCBI wrote the dump: the date of its names.dmp, in UTC."""
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        stamp = archive.getmember("names.dmp").mtime
    return datetime.fromtimestamp(stamp, tz=timezone.utc).date().isoformat()


def render(names: tuple[set[str], set[str], set[str]], source: str, md5: str, date: str) -> bytes:
    """The table as it is stored: gzip with no time stamp, so the same names give the same bytes."""
    scientific, vertebrate, common = names
    lines = [
        "# Names NCBI Taxonomy gives to animals. Made by tools/build_animal_names.py; do not edit by hand.",
        f"# Source: {source}, md5 {md5}",
        f"# NCBI Taxonomy of {date} (the date of names.dmp in the dump)",
        "# s = scientific name (a chordate at any rank; any other animal by a name of two or three words);",
        "# v = common name of a chordate; c = common name of another animal. People are not listed.",
    ]
    lines += [f"s\t{name}" for name in sorted(scientific)]
    lines += [f"v\t{name}" for name in sorted(vertebrate)]
    lines += [f"c\t{name}" for name in sorted(common)]
    return gzip.compress(("\n".join(lines) + "\n").encode("utf-8"), compresslevel=9, mtime=0)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Rebuild onehealth's table of animal names from NCBI Taxonomy.")
    parser.add_argument("dump", nargs="?", help="new_taxdump.tar.gz (default: download it from NCBI)")
    parser.add_argument("--out", type=Path, default=OUT, help=f"Where to write the table (default: {OUT})")
    args = parser.parse_args(argv)
    data, source = read_dump(args.dump)
    names = build(data)
    args.out.write_bytes(render(names, source, hashlib.md5(data).hexdigest(), dump_date(data)))
    scientific, vertebrate, common = names
    print(f"wrote {args.out}: {len(scientific)} scientific names, {len(vertebrate)} common names of chordates and "
          f"{len(common)} of other animals, {args.out.stat().st_size / 1e6:.2f} MB", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
