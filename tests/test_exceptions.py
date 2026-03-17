"""Tests for pdbe_sifts.base.exceptions."""

import pytest
from pdbe_sifts.base.exceptions import (
    AccessionNotFound,
    BatchRunException,
    EntryFailedException,
    EntryTimedOutException,
    NoSegmentsError,
    NotAPolyPeptide,
    ObsoleteUniProtError,
    ProcessFailedError,
    ReleaseCheckFailedException,
    SplitAccessionError,
    TooManyFailedEntries,
)


class TestSimpleExceptions:
    """Exceptions with no custom logic — just subclass Exception."""

    @pytest.mark.parametrize(
        "exc_class",
        [
            BatchRunException,
            ProcessFailedError,
            EntryFailedException,
            EntryTimedOutException,
            TooManyFailedEntries,
            ObsoleteUniProtError,
            AccessionNotFound,
            SplitAccessionError,
            NoSegmentsError,
            NotAPolyPeptide,
        ],
    )
    def test_is_exception(self, exc_class):
        assert issubclass(exc_class, Exception)

    @pytest.mark.parametrize(
        "exc_class",
        [
            BatchRunException,
            ProcessFailedError,
            EntryFailedException,
            EntryTimedOutException,
            TooManyFailedEntries,
            ObsoleteUniProtError,
            AccessionNotFound,
            SplitAccessionError,
            NoSegmentsError,
            NotAPolyPeptide,
        ],
    )
    def test_can_raise_and_catch(self, exc_class):
        msg = "test message"
        with pytest.raises(exc_class, match=msg):
            raise exc_class(msg)


class TestReleaseCheckFailedException:
    def test_attributes(self):
        exc = ReleaseCheckFailedException(
            entry_id="1ABC",
            message="check failed",
            context={"key": "val"},
        )
        assert exc.entry_id == "1ABC"
        assert exc.message == "check failed"
        assert exc.context == {"key": "val"}

    def test_str(self):
        exc = ReleaseCheckFailedException(entry_id="1ABC", message="bad entry")
        assert "bad entry" in str(exc)

    def test_context_defaults_to_none(self):
        exc = ReleaseCheckFailedException(entry_id="X", message="m")
        assert exc.context is None

    def test_is_exception(self):
        assert issubclass(ReleaseCheckFailedException, Exception)
