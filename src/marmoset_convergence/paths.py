"""Canonical, repository-relative paths used by notebooks and scripts.

This module performs no I/O at import time.  Directory creation is deliberately
an explicit action so that cached-mode analyses remain read-only by default.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Union

PathLike = Union[str, Path]


@dataclass(frozen=True)
class ProjectPaths:
    """Resolved paths for the public reproducibility layout."""

    root: Path

    @classmethod
    def from_root(cls, root: PathLike) -> "ProjectPaths":
        return cls(Path(root).expanduser().resolve())

    @classmethod
    def discover(cls, start: Optional[PathLike] = None) -> "ProjectPaths":
        """Find the nearest parent containing ``pyproject.toml`` and ``data``.

        Discovery is intended for interactive use.  Production scripts should
        pass ``--root`` so their behavior does not depend on the working
        directory.
        """

        cursor = Path(start or Path.cwd()).expanduser().resolve()
        for candidate in (cursor, *cursor.parents):
            if (candidate / "pyproject.toml").is_file():
                return cls(candidate)
        raise FileNotFoundError(
            f"Could not find a project root above {cursor}; pass an explicit root"
        )

    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def processed(self) -> Path:
        return self.data / "processed"

    @property
    def cache(self) -> Path:
        return self.data / "cache"

    @property
    def derived(self) -> Path:
        return self.data / "derived"

    @property
    def sequences(self) -> Path:
        return self.processed / "sequences.csv"

    @property
    def calls(self) -> Path:
        return self.processed / "calls.csv"

    @property
    def sequence_inventory(self) -> Path:
        return self.cache / "sequence_session_inventory.csv"

    @property
    def call_inventory(self) -> Path:
        return self.cache / "call_session_inventory.csv"

    @property
    def sequence_representations(self) -> Path:
        return self.cache / "sequence_representations.npz"

    @property
    def sequence_distances(self) -> Path:
        return self.cache / "sequence_distances.npz"

    @property
    def dtw_distances(self) -> Path:
        return self.cache / "dtw_distances_condensed.npz"

    @property
    def stp_pca_scores(self) -> Path:
        return self.cache / "stp_pca_scores.csv"

    @property
    def mfcc_pca_scores(self) -> Path:
        return self.cache / "mfcc_pca_scores.csv"

    @property
    def vae_latent_means(self) -> Path:
        return self.cache / "vae_latent_means.csv"

    @property
    def embeddings(self) -> Path:
        return self.cache / "embeddings"

    @property
    def manifest(self) -> Path:
        return self.cache / "manifest.json"

    @property
    def sequence_repertoire_distances(self) -> Path:
        return self.derived / "sequence_repertoire_distances.csv"

    @property
    def call_repertoire_distances(self) -> Path:
        return self.derived / "call_repertoire_distances.csv"

    @property
    def model_sequence(self) -> Path:
        return self.derived / "model_sequence.csv"

    @property
    def model_call(self) -> Path:
        return self.derived / "model_call.csv"

    @property
    def audio(self) -> Path:
        """Location of the WAV files in the current public repository."""

        return self.root / "extracted calls"

    def create_output_directories(self) -> None:
        """Create only generated-data directories, when explicitly requested."""

        for directory in (self.processed, self.cache, self.derived, self.embeddings):
            directory.mkdir(parents=True, exist_ok=True)

    def canonical_artifacts(self) -> Iterable[Path]:
        """Return the fixed artifacts used by the reader-facing notebooks."""

        return (
            self.sequences,
            self.calls,
            self.sequence_inventory,
            self.call_inventory,
            self.sequence_representations,
            self.sequence_distances,
            self.dtw_distances,
            self.stp_pca_scores,
            self.mfcc_pca_scores,
            self.vae_latent_means,
            self.sequence_repertoire_distances,
            self.call_repertoire_distances,
            self.model_sequence,
            self.model_call,
        )
