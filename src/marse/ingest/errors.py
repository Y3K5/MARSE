"""Errors raised while a table is read or checked."""


class IngestError(ValueError):
    """A table is missing a column, or a value fails a check."""
