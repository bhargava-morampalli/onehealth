"""Errors in what the user gave onehealth: a file, a table or an option."""

from __future__ import annotations


class InputError(ValueError):
    """A problem with the input that the command line reports in one line, without a traceback.

    Anything else that goes wrong is a bug in onehealth and keeps its full traceback.
    """
