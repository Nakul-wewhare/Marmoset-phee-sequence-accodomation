"""Canonical call-table cleaning and full-mode audio preflight checks."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Optional, Sequence, Tuple, Union

import pandas as pd

from .exceptions import AudioPreflightError, DataValidationError, DuplicateCallConflictError

PathLike = Union[str, Path]


def call_basename(value: Any) -> str:
    """Extract a filename from either POSIX or Windows-style path text."""

    if value is None or pd.isna(value):
        raise DataValidationError("Call filename cannot be missing")
    cleaned = str(value).strip().strip('"').strip("'").replace("\\", "/")
    name = PurePosixPath(cleaned).name
    if not name or name in {".", ".."}:
        raise DataValidationError(f"Invalid call filename/path: {value!r}")
    return name


def clean_call_path(
    value: Any,
    *,
    filename: Optional[Any] = None,
    audio_dir: Optional[PathLike] = "extracted calls",
) -> str:
    """Return a portable call path, optionally rebased under ``audio_dir``.

    Rebasing strips machine-specific absolute prefixes while retaining the
    validated basename.  Set ``audio_dir=None`` to keep the normalized path.
    """

    normalized = str(value).strip().strip('"').strip("'").replace("\\", "/")
    path_name = call_basename(normalized)
    expected_name = call_basename(filename) if filename is not None else path_name
    if path_name != expected_name:
        raise DataValidationError(
            f"Audio path basename {path_name!r} does not match filename {expected_name!r}"
        )
    if audio_dir is None:
        return PurePosixPath(normalized).as_posix()
    directory = str(audio_dir).strip().replace("\\", "/").rstrip("/")
    if not directory:
        return expected_name
    return (PurePosixPath(directory) / expected_name).as_posix()


def _same_scalar(left: Any, right: Any) -> bool:
    left_missing = pd.isna(left)
    right_missing = pd.isna(right)
    if bool(left_missing) or bool(right_missing):
        return bool(left_missing) and bool(right_missing)
    return bool(left == right)


def deduplicate_calls(
    calls: pd.DataFrame,
    *,
    filename_col: str = "filename",
    call_id_col: str = "call_id",
    path_col: Optional[str] = "audio_path",
    audio_dir: Optional[PathLike] = "extracted calls",
    ignore_columns: Sequence[str] = (),
) -> pd.DataFrame:
    """Deduplicate a call table by normalized filename.

    All non-path values in duplicate rows must be identical.  A disagreement is
    treated as a data error rather than silently selecting one copy.  The first
    occurrence determines row order; ``call_id`` is the normalized filename.
    """

    if filename_col not in calls:
        raise DataValidationError(f"Missing required column {filename_col!r}")
    result = calls.copy()
    result[filename_col] = result[filename_col].map(call_basename)

    folded = result[filename_col].str.casefold()
    case_collisions = (
        result.assign(_folded=folded)
        .groupby("_folded", sort=False)[filename_col]
        .nunique(dropna=False)
    )
    if (case_collisions > 1).any():
        names = case_collisions[case_collisions > 1].index.tolist()
        raise DataValidationError(
            "Call filenames differ only by case, which is not portable: "
            + ", ".join(names)
        )

    excluded = {filename_col, call_id_col, *ignore_columns}
    if path_col is not None:
        excluded.add(path_col)
    validation_columns = [column for column in result.columns if column not in excluded]

    conflicts = []
    for filename, group in result.groupby(filename_col, sort=False, dropna=False):
        if len(group) == 1:
            continue
        first = group.iloc[0]
        differing = [
            column
            for column in validation_columns
            if any(not _same_scalar(first[column], value) for value in group[column].iloc[1:])
        ]
        if differing:
            conflicts.append(f"{filename}: {', '.join(differing)}")
    if conflicts:
        raise DuplicateCallConflictError(
            "Rows sharing a call filename disagree in non-path data: "
            + "; ".join(conflicts)
        )

    result = result.drop_duplicates(filename_col, keep="first").copy()
    if call_id_col in result:
        result = result.drop(columns=[call_id_col])
    result.insert(0, call_id_col, result[filename_col])

    if path_col is not None:
        if path_col in result:
            cleaned_paths = []
            for filename, value in zip(result[filename_col], result[path_col]):
                if pd.isna(value) or not str(value).strip():
                    value = filename
                cleaned_paths.append(
                    clean_call_path(value, filename=filename, audio_dir=audio_dir)
                )
            result[path_col] = cleaned_paths
        elif audio_dir is not None:
            result[path_col] = [
                clean_call_path(filename, filename=filename, audio_dir=audio_dir)
                for filename in result[filename_col]
            ]
    return result.reset_index(drop=True)


@dataclass(frozen=True)
class AudioPreflight:
    """Read-only inventory comparison performed before full recomputation."""

    expected_count: int
    available_count: int
    missing: Tuple[str, ...]
    ambiguous: Tuple[str, ...]
    unexpected: Tuple[str, ...]

    @property
    def complete(self) -> bool:
        return not self.missing and not self.ambiguous

    def require_complete(self) -> None:
        if self.complete:
            return
        details = []
        if self.missing:
            preview = ", ".join(self.missing[:5])
            suffix = " ..." if len(self.missing) > 5 else ""
            details.append(f"missing {len(self.missing)} WAV(s): {preview}{suffix}")
        if self.ambiguous:
            details.append(
                "duplicate WAV basenames in audio tree: " + ", ".join(self.ambiguous[:5])
            )
        raise AudioPreflightError("Full-mode audio preflight failed; " + "; ".join(details))


def preflight_audio(
    calls: pd.DataFrame,
    audio_dir: PathLike,
    *,
    filename_col: str = "filename",
    call_id_col: str = "call_id",
) -> AudioPreflight:
    """Compare canonical call rows with available WAVs without reading audio.

    Cached mode should not call this function.  Full mode should call
    ``require_complete`` on the result before spectrogram or DTW work begins.
    """

    for column in (filename_col, call_id_col):
        if column not in calls:
            raise DataValidationError(f"Missing required column {column!r}")
    filenames = tuple(calls[filename_col].map(call_basename))
    call_ids = calls[call_id_col].astype(str)
    if call_ids.duplicated().any():
        duplicated = sorted(call_ids[call_ids.duplicated(keep=False)].unique())
        raise DataValidationError("Duplicate call_id values: " + ", ".join(duplicated[:5]))
    if len(set(filenames)) != len(filenames):
        raise DataValidationError("Canonical call table contains duplicate filenames")

    directory = Path(audio_dir).expanduser().resolve()
    if not directory.is_dir():
        return AudioPreflight(
            expected_count=len(filenames),
            available_count=0,
            missing=tuple(sorted(filenames)),
            ambiguous=(),
            unexpected=(),
        )

    by_name = {}
    for path in directory.rglob("*"):
        if path.is_file() and path.suffix.casefold() == ".wav":
            by_name.setdefault(path.name, []).append(path)
    expected = set(filenames)
    available = set(by_name)
    return AudioPreflight(
        expected_count=len(expected),
        available_count=sum(len(paths) for paths in by_name.values()),
        missing=tuple(sorted(expected - available)),
        ambiguous=tuple(sorted(name for name, paths in by_name.items() if len(paths) > 1)),
        unexpected=tuple(sorted(available - expected)),
    )
