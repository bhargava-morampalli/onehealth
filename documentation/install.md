# Install

`onehealth` needs Python 3.10 or newer. It depends on numpy, matplotlib and PyYAML, which are installed with it.

## With pip

```bash
git clone https://github.com/bhargava-morampalli/onehealth.git
cd onehealth
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

## With uv

```bash
git clone https://github.com/bhargava-morampalli/onehealth.git
cd onehealth
uv sync
uv run onehealth --help
```

## Check it works

```bash
onehealth --version
onehealth --help
```

Every command has its own help, for example `onehealth classify --help`.

## Optional extras

| Extra | Installs | For |
|---|---|---|
| `dev` | pytest | Running the test suite: `pip install -e ".[dev]"`, then `pytest` |
| `docs` | MkDocs Material | Building this site: `pip install -e ".[docs]"`, then `mkdocs serve` |

## Network access

Only `resolve` uses the network, to look IDs up at NCBI. It reads metadata only and never downloads sequences. Every other command works offline, and `resolve --no-network` uses only the local metadata you give it.
