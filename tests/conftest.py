from pathlib import Path

import pytest

from onehealth.lexicon import load_lexicon


@pytest.fixture(scope="session")
def lexicon():
    return load_lexicon()


@pytest.fixture
def fixtures() -> Path:
    return Path(__file__).parent / "fixtures"
