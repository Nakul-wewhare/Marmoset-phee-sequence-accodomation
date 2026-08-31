#!/usr/bin/env python3
"""Convert the historical wide CSV files to the public canonical schema.

This is an inexpensive, deterministic migration step. It does not read audio,
calculate distances, fit an embedding, train the VAE, or fit a statistical
model. The historical files remain the provenance inputs; downstream analysis
uses only the canonical tables written under ``data/processed``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import quote

import pandas as pd

from marmoset_convergence.calls import deduplicate_calls
from marmoset_convergence.identifiers import canonical_session_id, normalize_session_number


BONDED_PAIRS = {
    "A": ("Tabor", "Lola"),
    "B": ("Wuschel", "Olympia"),
    "C": ("Odin", "Nougatti"),
}
PAIR_BY_INDIVIDUAL = {
    individual: pair_id
    for pair_id, individuals in BONDED_PAIRS.items()
    for individual in individuals
}
CALL_TYPES = frozenset("AKPRST")
STP_COLUMNS = {
    "Dur 90% (s)": "duration_90_s",
    "Dur 50% (s)": "duration_50_s",
    "Center Freq (Hz)": "center_frequency_hz",
    "Freq 5% (Hz)": "frequency_5_hz",
    "Freq 25% (Hz)": "frequency_25_hz",
    "Freq 75% (Hz)": "frequency_75_hz",
    "Freq 95% (Hz)": "frequency_95_hz",
    "BW 50% (Hz)": "bandwidth_50_hz",
    "BW 90% (Hz)": "bandwidth_90_hz",
    "Avg Entropy (bits)": "average_entropy_bits",
    "Agg Entropy (bits)": "aggregate_entropy_bits",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def context_label(value: object) -> str:
    normalized = str(value).strip().lower().replace("_", "-")
    aliases = {
        "partner": "partner",
        "paired": "partner",
        "non-partner": "non-partner",
        "nonpartner": "non-partner",
        "non-paired": "non-partner",
        "stranger": "non-partner",
    }
    if normalized not in aliases:
        raise ValueError(f"Unknown audience-context label: {value!r}")
    return aliases[normalized]


def pair_id_for(individual: object) -> str:
    name = str(individual).strip()
    if name not in PAIR_BY_INDIVIDUAL:
        raise ValueError(f"No eventual pair is registered for {name!r}")
    return PAIR_BY_INDIVIDUAL[name]


def clean_sequence(value: object) -> str:
    """Remove one known legacy concatenation sentinel and validate symbols."""

    sequence = str(value).strip().replace("nan", "").replace("x", "").replace("y", "")
    if len(sequence) < 2:
        raise ValueError(f"Sequence has fewer than two retained elements: {value!r}")
    unknown = sorted(set(sequence) - CALL_TYPES)
    if unknown:
        raise ValueError(f"Sequence contains unknown call symbols {unknown}: {value!r}")
    return sequence


def recording_id(row: pd.Series) -> str:
    participants = sorted((str(row["individual_id"]), str(row["receiver_id"])))
    components = (
        f"participants={quote('+'.join(participants), safe='-._~')}",
        f"stage={quote(str(row['stage']), safe='-._~')}",
        f"session={quote(str(row['session_number']), safe='-._~')}",
    )
    return "|".join(components)


def add_shared_identifiers(table: pd.DataFrame) -> pd.DataFrame:
    result = table.copy()
    result["individual_id"] = result["focal ID"].astype(str).str.strip()
    result["receiver_id"] = result["conspecific_ID"].astype(str).str.strip()
    result["stage"] = result["stage"].astype(str).str.strip().str.lower()
    result["context"] = result["paired_status"].map(context_label)
    result["session_number"] = result["session_number"].map(normalize_session_number)
    result["pair_id"] = result["individual_id"].map(pair_id_for)
    result["repertoire_id"] = [
        canonical_session_id(focal, receiver, stage, session)
        for focal, receiver, stage, session in zip(
            result["individual_id"],
            result["receiver_id"],
            result["stage"],
            result["session_number"],
        )
    ]
    result["recording_id"] = result.apply(recording_id, axis=1)
    datetime_text = (
        result["date_x"].astype(str).str.strip()
        + " "
        + result["time_x"].astype(str).str.strip().str.replace("-", ":", regex=False)
    )
    parsed = pd.to_datetime(datetime_text, errors="coerce")
    if parsed.isna().any():
        raise ValueError("One or more recording date/time values could not be parsed")
    result["recording_datetime"] = parsed.dt.strftime("%Y-%m-%dT%H:%M:%S")
    return result


def migrate_sequences(source: Path) -> pd.DataFrame:
    raw = pd.read_csv(source)
    required = {
        "0",
        "true_sequence",
        "seq_num",
        "stage",
        "date_x",
        "time_x",
        "focal ID",
        "session_number",
        "conspecific_ID",
        "paired_status",
    }
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"Historical sequence table is missing columns: {missing}")

    table = add_shared_identifiers(raw)
    table["sequence"] = table["true_sequence"].map(clean_sequence)
    table["n_elements"] = table["sequence"].str.len()
    table["source_file"] = table["0"].map(lambda value: Path(str(value)).name)
    sequence_number = table["seq_num"].map(normalize_session_number)
    table["sequence_id"] = table["source_file"] + "::sequence=" + sequence_number
    if table["sequence_id"].duplicated().any():
        raise ValueError("Canonical sequence identifiers are not unique")

    columns = [
        "sequence_id",
        "repertoire_id",
        "recording_id",
        "individual_id",
        "receiver_id",
        "pair_id",
        "stage",
        "context",
        "sequence",
        "n_elements",
        "session_number",
        "recording_datetime",
        "source_file",
    ]
    result = table[columns].copy()
    if len(result) != 1619:
        raise ValueError(f"Expected 1,619 sequences, observed {len(result):,}")
    return result


def migrate_calls(source: Path) -> pd.DataFrame:
    raw = pd.read_csv(source)
    required = {
        "filename",
        "call_audio_path",
        "stage",
        "date_x",
        "time_x",
        "focal ID",
        "session_number",
        "conspecific_ID",
        "paired_status",
        *STP_COLUMNS,
        *(str(index) for index in range(1, 133)),
    }
    missing = sorted(required - set(raw.columns))
    if missing:
        raise ValueError(f"Historical call table is missing columns: {missing}")

    table = raw.rename(columns={"call_audio_path": "audio_path", **STP_COLUMNS})
    table = deduplicate_calls(
        table,
        filename_col="filename",
        call_id_col="call_id",
        path_col="audio_path",
        audio_dir="extracted calls",
    )
    table = add_shared_identifiers(table)
    mfcc_rename = {
        str(index): f"mfcc_{index:03d}" for index in range(1, 133)
    }
    table = table.rename(columns=mfcc_rename)
    columns = [
        "call_id",
        "filename",
        "repertoire_id",
        "recording_id",
        "individual_id",
        "receiver_id",
        "pair_id",
        "stage",
        "context",
        "session_number",
        "recording_datetime",
        "audio_path",
        *STP_COLUMNS.values(),
        *mfcc_rename.values(),
    ]
    result = table[columns].copy()
    if len(result) != 3612:
        raise ValueError(f"Expected 3,612 unique calls, observed {len(result):,}")
    if result["call_id"].duplicated().any():
        raise ValueError("Canonical call identifiers are not unique")
    return result


def write_manifest(
    path: Path,
    sequence_source: Path,
    call_source: Path,
    sequence_output: Path,
    call_output: Path,
) -> None:
    project_root = path.resolve().parents[2]

    def repository_path(file_path: Path) -> str:
        try:
            return file_path.resolve().relative_to(project_root).as_posix()
        except ValueError as exc:
            raise ValueError(
                f"Manifest input is outside the repository: {file_path}"
            ) from exc

    payload = {
        "schema_version": "1.0",
        "kind": "processed-input-manifest",
        "migration": "scripts/import_legacy_processed_data.py",
        "parameters": {
            "bonded_pairs": BONDED_PAIRS,
            "sequence_alphabet": sorted(CALL_TYPES),
            "sequence_sentinel_rule": "remove literal nan/x/y concatenation sentinels",
            "mfcc_columns": 132,
            "duplicate_call_rule": "identical rows by normalized filename; keep first",
            "audio_path_base": "extracted calls",
            "legacy_call_audio_path_rule": (
                "author-machine prefixes removed by "
                "scripts/sanitize_legacy_source_paths.py"
            ),
        },
        "sources": [
            {
                "path": repository_path(sequence_source),
                "sha256": sha256_file(sequence_source),
            },
            {
                "path": repository_path(call_source),
                "sha256": sha256_file(call_source),
            },
        ],
        "artifacts": [
            {
                "path": repository_path(sequence_output),
                "rows": 1619,
                "sha256": sha256_file(sequence_output),
            },
            {
                "path": repository_path(call_output),
                "rows": 3612,
                "sha256": sha256_file(call_output),
            },
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--sequence-source",
        type=Path,
        default=Path("data/source/legacy_processed_sequences.csv"),
    )
    result.add_argument(
        "--call-source",
        type=Path,
        default=Path("data/source/legacy_processed_calls.csv"),
    )
    result.add_argument(
        "--output-dir", type=Path, default=Path("data/processed")
    )
    return result


def main() -> None:
    args = parser().parse_args()
    sequence_source = args.sequence_source.resolve()
    call_source = args.call_source.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    sequence_output = output_dir / "sequences.csv"
    call_output = output_dir / "calls.csv"

    sequences = migrate_sequences(sequence_source)
    calls = migrate_calls(call_source)
    sequences.to_csv(sequence_output, index=False)
    calls.to_csv(call_output, index=False)
    write_manifest(
        output_dir / "manifest.json",
        sequence_source,
        call_source,
        sequence_output,
        call_output,
    )
    print(f"Wrote {len(sequences):,} sequences to {sequence_output}")
    print(f"Wrote {len(calls):,} unique calls to {call_output}")


if __name__ == "__main__":
    main()
