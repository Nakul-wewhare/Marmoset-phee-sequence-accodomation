from pathlib import Path

import pandas as pd
import pytest

from marmoset_convergence import (
    AudioPreflightError,
    DataValidationError,
    DuplicateCallConflictError,
    clean_call_path,
    deduplicate_calls,
    preflight_audio,
)


def test_clean_call_path_rebases_windows_path_to_portable_relative_path():
    assert clean_call_path(
        r"C:\\old\\machine\\calls\\a.wav",
        filename="a.wav",
        audio_dir="extracted calls",
    ) == "extracted calls/a.wav"


def test_deduplicate_calls_keeps_identical_copy_and_cleans_path():
    calls = pd.DataFrame(
        {
            "filename": [r"old\\a.wav", "a.wav", "b.wav"],
            "individual_id": ["A", "A", "B"],
            "value": [1.0, 1.0, 2.0],
            "call_audio_path": [r"X:\\audio\\a.wav", "/tmp/a.wav", "/tmp/b.wav"],
        }
    )

    result = deduplicate_calls(
        calls,
        path_col="call_audio_path",
        audio_dir="extracted calls",
    )

    assert result["call_id"].tolist() == ["a.wav", "b.wav"]
    assert result["call_audio_path"].tolist() == [
        "extracted calls/a.wav",
        "extracted calls/b.wav",
    ]


def test_deduplicate_calls_rejects_conflicting_duplicate_metadata():
    calls = pd.DataFrame(
        {
            "filename": ["a.wav", "a.wav"],
            "individual_id": ["A", "B"],
        }
    )
    with pytest.raises(DuplicateCallConflictError, match="individual_id"):
        deduplicate_calls(calls, path_col=None)


def test_deduplicate_calls_rejects_case_only_collision():
    calls = pd.DataFrame({"filename": ["a.wav", "A.wav"], "value": [1, 1]})
    with pytest.raises(DataValidationError, match="differ only by case"):
        deduplicate_calls(calls, path_col=None)


def test_audio_preflight_reports_missing_without_reading_wavs(tmp_path: Path):
    audio = tmp_path / "audio"
    audio.mkdir()
    (audio / "a.wav").write_bytes(b"not decoded")
    calls = pd.DataFrame(
        {"call_id": ["a.wav", "b.wav"], "filename": ["a.wav", "b.wav"]}
    )

    report = preflight_audio(calls, audio)

    assert report.missing == ("b.wav",)
    with pytest.raises(AudioPreflightError, match="missing 1 WAV"):
        report.require_complete()


def test_audio_preflight_accepts_one_file_per_canonical_call(tmp_path: Path):
    audio = tmp_path / "audio"
    audio.mkdir()
    (audio / "a.wav").write_bytes(b"a")
    (audio / "b.WAV").write_bytes(b"b")
    calls = pd.DataFrame(
        {"call_id": ["a.wav", "b.WAV"], "filename": ["a.wav", "b.WAV"]}
    )

    report = preflight_audio(calls, audio)

    assert report.complete
    report.require_complete()
