# Install

`onehealth` needs Python 3.10 or newer. It depends on numpy, matplotlib and PyYAML, which are installed with it.

## Install

No clone is needed. In a virtual environment, or any environment you use for analysis:

```bash
pip install "git+https://github.com/bhargava-morampalli/onehealth.git"
```

This installs the newest version and puts the `onehealth` command on your path.

### A fixed version

For a study, install one release and name it in your methods, so that the same command gives the same classes later. Releases are tagged `v` and the version number; the [tags page](https://github.com/bhargava-morampalli/onehealth/tags) lists them. Add `@` and the tag:

```bash
pip install "git+https://github.com/bhargava-morampalli/onehealth.git@v0.2.0"
```

`onehealth --version` prints the installed version, and every classified row records the lexicon version that made the call (`lexicon_version`). The [changelog](changelog.md) says what each version changed; versions that change what a class means say so at the top of their entry.

### As a command of its own

[uv](https://docs.astral.sh/uv/) and [pipx](https://pipx.pypa.io/) install a command in an environment of its own, so that it does not meet the packages of your other projects:

```bash
uv tool install "git+https://github.com/bhargava-morampalli/onehealth.git"
```

```bash
pipx install "git+https://github.com/bhargava-morampalli/onehealth.git"
```

Both take `@v0.2.0` at the end of the address for a fixed version.

### From a clone, to change the code or run the tests

```bash
git clone https://github.com/bhargava-morampalli/onehealth.git
cd onehealth
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Or with uv:

```bash
git clone https://github.com/bhargava-morampalli/onehealth.git
cd onehealth
uv sync --extra dev
uv run onehealth --help
```

The examples in the [tutorial](tutorial.md) use files from the repository (`examples/`), so they need a clone too.

## Check it works

```bash
onehealth --version
onehealth --help
```

Every command has its own help, for example `onehealth classify --help`.

## Optional extras

| Extra | Installs | For |
|---|---|---|
| `dev` | pytest, scipy, ruff, mypy | Running the test suite and the checks: `pip install -e ".[dev]"`, then `pytest`, `ruff check` and `mypy`. scipy is used only by the tests, to check the statistics. |
| `docs` | MkDocs Material | Building this site: `pip install -e ".[docs]"`, then `mkdocs serve` |

## Network access

Only `resolve` uses the network, to look IDs up at NCBI. It reads metadata only and never downloads sequences. Every other command works offline, and `resolve --no-network` uses only the local metadata you give it.
