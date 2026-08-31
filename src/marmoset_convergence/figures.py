"""Cheap analytical figure regeneration from validated tables and coordinates.

This module never reads audio, calculates a distance, or fits an embedding.
Every export path is derived from a locked manuscript path by adding the
``_regenerated`` suffix, so the provenance-locked snapshot cannot be replaced.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log, pi
from pathlib import Path
from typing import Mapping, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, FancyArrowPatch, Patch

from .exceptions import DataValidationError
from .sessions import normalize_context
from .identifiers import normalize_stage

PathLike = Union[str, Path]

CALL_TYPE_LABELS: Mapping[str, str] = {
    "A": "phee",
    "P": "peep",
    "R": "trillphee",
    "K": "egg",
    "S": "sd-peep",
    "T": "tsk",
}
CALL_TYPE_ORDER = ("A", "P", "R", "K", "S", "T")
EMBEDDING_COORDINATE_COLUMNS = ("embedding_x", "embedding_y")
EMBEDDING_GROUP_COLUMNS = ("individual_id", "pair_id", "stage", "context")

INDIVIDUAL_COLORS: Mapping[str, str] = {
    "Tabor": "#1f77b4",
    "Lola": "#8ecae6",
    "Wuschel": "#7f4a3b",
    "Olympia": "#c9a58f",
    "Odin": "#e66101",
    "Nougatti": "#fdb863",
}


@dataclass(frozen=True)
class EmbeddingFigureSpec:
    code: str
    metric: str
    title: str
    level: str
    id_column: str
    locked_relative_path: str

    def __post_init__(self) -> None:
        if self.level not in {"sequence", "call"}:
            raise DataValidationError("Embedding figure level must be sequence or call")

    @property
    def annotation_note(self) -> str:
        if self.level == "call":
            return (
                "Linked spectrogram cards in the locked manuscript snapshot are "
                "omitted: the complete WAV collection is not part of the cache-only "
                "figure contract. This export reads pooled cached coordinates only."
            )
        return (
            "Linked repertoire cards in the locked manuscript snapshot are omitted "
            "from this compact analytical export. All panels reuse the pooled cached "
            "coordinates and fixed limits."
        )


EMBEDDING_FIGURES: Mapping[str, EmbeddingFigureSpec] = {
    "4": EmbeddingFigureSpec(
        "4",
        "transition_probability",
        "Transition probability space",
        "sequence",
        "repertoire_id",
        "results/figures/main/figure_4.png",
    ),
    "s2": EmbeddingFigureSpec(
        "s2",
        "bigram",
        "Bigram frequency space",
        "sequence",
        "repertoire_id",
        "results/figures/supplement/figure_s2_bigram.png",
    ),
    "s3": EmbeddingFigureSpec(
        "s3",
        "local_alignment",
        "Local-alignment sequence space",
        "sequence",
        "repertoire_id",
        "results/figures/supplement/figure_s3_local_alignment.png",
    ),
    "s4": EmbeddingFigureSpec(
        "s4",
        "phee_repeat",
        "Phee-repeat sequence space",
        "sequence",
        "repertoire_id",
        "results/figures/supplement/figure_s4_phee_repeat.png",
    ),
    "s5": EmbeddingFigureSpec(
        "s5",
        "mfcc",
        "MFCC acoustic space",
        "call",
        "call_id",
        "results/figures/supplement/figure_s5_mfcc.png",
    ),
    "s6": EmbeddingFigureSpec(
        "s6",
        "stp",
        "Spectro-temporal acoustic space",
        "call",
        "call_id",
        "results/figures/supplement/figure_s6_stp.png",
    ),
    "s7": EmbeddingFigureSpec(
        "s7",
        "dtw",
        "Dynamic-time-warping acoustic space",
        "call",
        "call_id",
        "results/figures/supplement/figure_s7_dtw.png",
    ),
    "s8": EmbeddingFigureSpec(
        "s8",
        "vae",
        "VAE acoustic space",
        "call",
        "call_id",
        "results/figures/supplement/figure_s8_vae.png",
    ),
}


def regenerated_figure_path(locked_path: PathLike) -> Path:
    """Derive an unmistakably separate output beside a locked snapshot."""

    locked = Path(locked_path)
    if locked.stem.endswith("_regenerated"):
        raise DataValidationError(
            "Pass the locked manuscript path; the regenerated suffix is added automatically"
        )
    if locked.suffix.casefold() != ".png":
        raise DataValidationError("Reader-facing regenerated exports must be PNG files")
    return locked.with_name(f"{locked.stem}_regenerated{locked.suffix}")


def save_regenerated_figure(
    figure: plt.Figure,
    locked_path: PathLike,
    *,
    overwrite_regenerated: bool = False,
    dpi: int = 200,
) -> Path:
    """Save only to the derived ``_regenerated`` path, never ``locked_path``."""

    locked = Path(locked_path)
    destination = regenerated_figure_path(locked)
    if destination.exists() and not overwrite_regenerated:
        raise FileExistsError(
            f"Regenerated export already exists: {destination}. Pass explicit "
            "overwrite_regenerated=True to replace only that derived file."
        )
    if locked.exists() and destination.exists():
        try:
            if locked.resolve() == destination.resolve():
                raise DataValidationError("Regenerated output resolves to the locked snapshot")
        except OSError:
            pass
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=dpi, bbox_inches="tight", facecolor="white")
    return destination


def _sequence_values(sequences: pd.DataFrame) -> Tuple[str, ...]:
    if "sequence" not in sequences:
        raise DataValidationError("Processed sequence table lacks the 'sequence' column")
    values = tuple(sequences["sequence"].astype(str).str.strip())
    if not values or any(not value for value in values):
        raise DataValidationError("Processed sequences must be nonempty")
    unknown = sorted(set("".join(values)).difference(CALL_TYPE_LABELS))
    if unknown:
        raise DataValidationError(f"Processed sequences contain unknown symbols: {unknown}")
    return values


def _panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.14,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=18,
        fontweight="bold",
        va="top",
    )


def _draw_transition_network(
    axis: plt.Axes,
    transition: np.ndarray,
    counts: np.ndarray,
    symbols: Sequence[str],
) -> None:
    angles = np.linspace(0, 2 * pi, len(symbols), endpoint=False)
    positions = {
        symbol: np.array([np.cos(angle), np.sin(angle)])
        for symbol, angle in zip(symbols, angles)
    }
    maximum = float(transition.max()) if transition.size else 1.0
    for i, source in enumerate(symbols):
        for j, target in enumerate(symbols):
            probability = transition[i, j]
            if probability <= 0:
                continue
            start = positions[source]
            end = positions[target]
            width = 0.4 + 4.0 * probability / maximum
            if source == target:
                start = start + np.array([-0.08, 0.04])
                end = end + np.array([0.08, 0.04])
                connection = "arc3,rad=-1.4"
            else:
                connection = "arc3,rad=0.08"
            arrow = FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=8,
                connectionstyle=connection,
                color="#444444",
                linewidth=width,
                alpha=0.85,
                shrinkA=10,
                shrinkB=10,
                zorder=1,
            )
            axis.add_patch(arrow)
    node_sizes = 180 + 850 * counts / counts.max()
    colors = plt.colormaps["viridis"](np.linspace(0.05, 0.95, len(symbols)))
    for symbol, size, color in zip(symbols, node_sizes, colors):
        x, y = positions[symbol]
        axis.scatter(x, y, s=size, color=color, edgecolor="black", zorder=3)
        axis.text(
            x * 1.20,
            y * 1.20,
            CALL_TYPE_LABELS[symbol],
            ha="center",
            va="center",
            fontsize=10,
        )
    axis.set_xlim(-1.45, 1.45)
    axis.set_ylim(-1.35, 1.35)
    axis.set_aspect("equal")
    axis.axis("off")


def plot_sequence_summary(sequences: pd.DataFrame) -> plt.Figure:
    """Build the four-panel analytical Figure 2 from processed sequences only."""

    values = _sequence_values(sequences)
    symbols = tuple(symbol for symbol in CALL_TYPE_ORDER if symbol in set("".join(values)))
    index = {symbol: position for position, symbol in enumerate(symbols)}
    call_counts = np.zeros(len(symbols), dtype=int)
    bigram_counts = np.zeros((len(symbols), len(symbols)), dtype=int)
    for sequence in values:
        for symbol in sequence:
            call_counts[index[symbol]] += 1
        for preceding, following in zip(sequence[:-1], sequence[1:]):
            bigram_counts[index[preceding], index[following]] += 1
    row_totals = bigram_counts.sum(axis=1, keepdims=True)
    transition = np.divide(
        bigram_counts,
        row_totals,
        out=np.zeros_like(bigram_counts, dtype=float),
        where=row_totals != 0,
    )
    lengths = pd.Series([len(value) for value in values]).value_counts().sort_index()

    figure, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_composition, ax_length, ax_network, ax_heatmap = axes.ravel()

    colors = plt.colormaps["viridis"](np.linspace(0.05, 0.98, len(symbols)))
    x = np.arange(len(symbols))
    bars = ax_composition.bar(x, call_counts, color=colors, width=0.78)
    ax_composition.set_yscale("log")
    ax_composition.set_xticks(x, [CALL_TYPE_LABELS[symbol] for symbol in symbols], rotation=28)
    ax_composition.set_ylabel("Count (log scale)")
    ax_composition.set_xlabel("Call type")
    ax_composition.set_title("Call-type composition")
    total_calls = call_counts.sum()
    for bar, count in zip(bars, call_counts):
        ax_composition.text(
            bar.get_x() + bar.get_width() / 2,
            count * 1.15,
            f"{100 * count / total_calls:.1f}%",
            ha="center",
            fontsize=10,
        )
    _panel_label(ax_composition, "A")

    ax_length.plot(
        lengths.index,
        lengths.values,
        color="black",
        linewidth=1.5,
        marker="o",
        markerfacecolor="#9acd32",
        markeredgecolor="black",
    )
    ax_length.set_yscale("log")
    ax_length.set_xlabel("Sequence length")
    ax_length.set_ylabel("Count (log scale)")
    ax_length.set_title("Sequence length distribution")
    for length, count in lengths.items():
        ax_length.annotate(
            str(int(count)),
            (length, count),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            fontsize=10,
        )
    _panel_label(ax_length, "B")

    _draw_transition_network(ax_network, transition, call_counts, symbols)
    ax_network.set_title("Transition network")
    _panel_label(ax_network, "C")

    display_counts = bigram_counts.T + 1
    image = ax_heatmap.imshow(
        display_counts,
        origin="lower",
        aspect="equal",
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max(1, int(display_counts.max()))),
    )
    labels = [CALL_TYPE_LABELS[symbol] for symbol in symbols]
    ax_heatmap.set_xticks(range(len(symbols)), labels, rotation=38, ha="right")
    ax_heatmap.set_yticks(range(len(symbols)), labels)
    ax_heatmap.set_xlabel("First element")
    ax_heatmap.set_ylabel("Second element")
    ax_heatmap.set_title("Bigram count (+1 for display)")
    figure.colorbar(image, ax=ax_heatmap, fraction=0.046, pad=0.04, label="Count + 1")
    _panel_label(ax_heatmap, "D")

    figure.suptitle("Structural features of phee sequences", fontsize=18, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    return figure


def regenerate_sequence_summary(
    sequences: pd.DataFrame,
    locked_path: PathLike,
    *,
    overwrite_regenerated: bool = False,
) -> Path:
    figure = plot_sequence_summary(sequences)
    try:
        return save_regenerated_figure(
            figure,
            locked_path,
            overwrite_regenerated=overwrite_regenerated,
        )
    finally:
        plt.close(figure)


def resolve_embedding_columns(
    embedding: pd.DataFrame,
    *,
    x_column: Optional[str] = None,
    y_column: Optional[str] = None,
) -> Tuple[str, str]:
    """Resolve common cached coordinate names without guessing from metadata."""

    if (x_column is None) != (y_column is None):
        raise DataValidationError("Specify both x_column and y_column, or neither")
    if x_column is not None:
        missing = [column for column in (x_column, y_column) if column not in embedding]
        if missing:
            raise DataValidationError("Embedding lacks coordinate columns: " + ", ".join(missing))
        return x_column, str(y_column)

    normalized = {
        "".join(character for character in str(column).casefold() if character.isalnum()): column
        for column in embedding.columns
    }
    candidates = (
        ("embeddingx", "embeddingy"),
        ("embedding1", "embedding2"),
        ("pacmap1", "pacmap2"),
        ("umap1", "umap2"),
        ("dimension1", "dimension2"),
        ("dim1", "dim2"),
        ("x", "y"),
    )
    matches = [
        (normalized[left], normalized[right])
        for left, right in candidates
        if left in normalized and right in normalized
    ]
    if len(matches) != 1:
        raise DataValidationError(
            "Could not identify exactly one coordinate pair; use explicit x_column/y_column"
        )
    return matches[0]


def prepare_embedding_table(
    embedding: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    id_column: str,
    x_column: Optional[str] = None,
    y_column: Optional[str] = None,
) -> Tuple[pd.DataFrame, str, str]:
    """Validate a canonical pooled-coordinate cache against its source rows.

    Public embedding CSVs deliberately have one exact seven-column schema: the
    observation identifier, ``embedding_x``/``embedding_y``, and the four
    grouping columns used in the panels.  Requiring row order and grouping
    values to match the canonical metadata prevents a set-equal but misaligned
    coordinate cache from producing a plausible-looking figure.
    """

    required_metadata = {id_column, "individual_id", "pair_id", "stage", "context"}
    missing_metadata = sorted(required_metadata.difference(metadata.columns))
    if missing_metadata:
        raise DataValidationError(
            "Metadata lacks columns required for embedding figures: "
            + ", ".join(missing_metadata)
        )
    if x_column is not None or y_column is not None:
        if (x_column, y_column) != EMBEDDING_COORDINATE_COLUMNS:
            raise DataValidationError(
                "Canonical embedding coordinates must be named embedding_x and embedding_y"
            )
    x_name, y_name = EMBEDDING_COORDINATE_COLUMNS
    expected_columns = {
        id_column,
        *EMBEDDING_COORDINATE_COLUMNS,
        *EMBEDDING_GROUP_COLUMNS,
    }
    observed_columns = set(embedding.columns)
    missing_columns = sorted(expected_columns.difference(observed_columns))
    extra_columns = sorted(observed_columns.difference(expected_columns))
    if missing_columns or extra_columns:
        details = []
        if missing_columns:
            details.append("missing " + ", ".join(missing_columns))
        if extra_columns:
            details.append("unexpected " + ", ".join(extra_columns))
        raise DataValidationError(
            "Embedding must use the exact canonical seven-column schema: "
            + "; ".join(details)
        )

    table = embedding[
        [id_column, x_name, y_name, *EMBEDDING_GROUP_COLUMNS]
    ].copy()
    canonical_metadata = metadata[
        [id_column, *EMBEDDING_GROUP_COLUMNS]
    ].copy()
    if table[[id_column, *EMBEDDING_GROUP_COLUMNS]].isna().any().any():
        raise DataValidationError("Embedding identifiers and grouping values must be nonmissing")
    if canonical_metadata.isna().any().any():
        raise DataValidationError("Canonical metadata contains missing grouping values")
    table[id_column] = table[id_column].astype(str)
    canonical_metadata[id_column] = canonical_metadata[id_column].astype(str)
    if table[id_column].duplicated().any() or canonical_metadata[id_column].duplicated().any():
        raise DataValidationError(f"{id_column} values must be unique in both tables")
    observed_ids = tuple(table[id_column])
    expected_ids = tuple(canonical_metadata[id_column])
    if observed_ids != expected_ids:
        mismatch = next(
            (
                position
                for position, (observed, expected) in enumerate(
                    zip(observed_ids, expected_ids)
                )
                if observed != expected
            ),
            min(len(observed_ids), len(expected_ids)),
        )
        raise DataValidationError(
            f"Embedding {id_column} order differs from canonical metadata at row {mismatch}"
        )
    table[x_name] = pd.to_numeric(table[x_name], errors="coerce")
    table[y_name] = pd.to_numeric(table[y_name], errors="coerce")
    if not np.isfinite(table[[x_name, y_name]].to_numpy()).all():
        raise DataValidationError("Embedding coordinates must be finite numeric values")
    table["stage"] = table["stage"].map(normalize_stage)
    table["context"] = table["context"].map(normalize_context)
    table["individual_id"] = table["individual_id"].astype(str)
    table["pair_id"] = table["pair_id"].astype(str)
    canonical_metadata["stage"] = canonical_metadata["stage"].map(normalize_stage)
    canonical_metadata["context"] = canonical_metadata["context"].map(normalize_context)
    canonical_metadata["individual_id"] = canonical_metadata["individual_id"].astype(str)
    canonical_metadata["pair_id"] = canonical_metadata["pair_id"].astype(str)
    for column in EMBEDDING_GROUP_COLUMNS:
        observed = table[column].to_numpy(dtype=str)
        expected = canonical_metadata[column].to_numpy(dtype=str)
        if not np.array_equal(observed, expected):
            mismatch = int(np.flatnonzero(observed != expected)[0])
            raise DataValidationError(
                f"Embedding grouping column {column!r} differs from canonical "
                f"metadata at row {mismatch}"
            )
    return table, x_name, y_name


def pooled_axis_limits(values: Sequence[float], *, padding_fraction: float = 0.05) -> Tuple[float, float]:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or not len(array) or not np.isfinite(array).all():
        raise DataValidationError("Axis coordinates must be one finite nonempty vector")
    minimum = float(array.min())
    maximum = float(array.max())
    span = maximum - minimum
    padding = (span if span > 0 else 1.0) * padding_fraction
    return minimum - padding, maximum + padding


def _color_map(individuals: Sequence[str]) -> Mapping[str, object]:
    ordered = tuple(sorted(set(individuals)))
    fallback = plt.colormaps["tab20"](np.linspace(0, 1, max(1, len(ordered))))
    return {
        individual: INDIVIDUAL_COLORS.get(individual, fallback[position])
        for position, individual in enumerate(ordered)
    }


def _add_gaussian_covariance_ellipse(
    axis: plt.Axes,
    points: pd.DataFrame,
    x_column: str,
    y_column: str,
    color: object,
) -> None:
    if len(points) < 3:
        return
    coordinates = points[[x_column, y_column]].to_numpy(dtype=float)
    covariance = np.cov(coordinates, rowvar=False, ddof=1)
    if covariance.shape != (2, 2) or not np.isfinite(covariance).all():
        return
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    eigenvalues = np.maximum(eigenvalues, 0)
    if eigenvalues.max() <= 1e-15:
        return
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    central_50_scale = np.sqrt(2 * log(2))
    width, height = 2 * central_50_scale * np.sqrt(eigenvalues)
    center = coordinates.mean(axis=0)
    ellipse = Ellipse(
        center,
        width,
        height,
        angle=angle,
        facecolor=color,
        edgecolor=color,
        linewidth=1.1,
        alpha=0.12,
        zorder=1,
    )
    axis.add_patch(ellipse)


def _plot_embedding_points(
    axis: plt.Axes,
    table: pd.DataFrame,
    x_column: str,
    y_column: str,
    colors: Mapping[str, object],
    *,
    point_size: float,
    show_summaries: bool,
) -> None:
    marker_by_context = {"partner": "^", "non-partner": "s"}
    for (individual, context), group in table.groupby(
        ["individual_id", "context"], sort=True, observed=True
    ):
        axis.scatter(
            group[x_column],
            group[y_column],
            s=point_size,
            marker=marker_by_context[context],
            color=colors[individual],
            alpha=0.62,
            edgecolors="none",
            zorder=2,
        )
    if not show_summaries or table.empty:
        return
    centroids = (
        table.groupby(["pair_id", "individual_id"], sort=True, observed=True)[
            [x_column, y_column]
        ]
        .mean()
        .reset_index()
    )
    for _, group in centroids.groupby("pair_id", sort=True):
        if len(group) == 2:
            axis.plot(
                group[x_column],
                group[y_column],
                color="#555555",
                linewidth=1.4,
                zorder=3,
            )
    for _, centroid in centroids.iterrows():
        individual = centroid["individual_id"]
        axis.scatter(
            centroid[x_column],
            centroid[y_column],
            s=125,
            marker="o",
            color=colors[individual],
            edgecolor="black",
            linewidth=1.0,
            zorder=4,
        )
    for individual, group in table.groupby("individual_id", sort=True, observed=True):
        _add_gaussian_covariance_ellipse(
            axis, group, x_column, y_column, colors[individual]
        )


def plot_embedding_space(
    embedding: pd.DataFrame,
    metadata: pd.DataFrame,
    spec: EmbeddingFigureSpec,
    *,
    x_column: Optional[str] = None,
    y_column: Optional[str] = None,
) -> plt.Figure:
    """Plot pooled and faceted views from one already-fitted coordinate cache."""

    table, x_name, y_name = prepare_embedding_table(
        embedding,
        metadata,
        id_column=spec.id_column,
        x_column=x_column,
        y_column=y_column,
    )
    x_limits = pooled_axis_limits(table[x_name])
    y_limits = pooled_axis_limits(table[y_name])
    colors = _color_map(table["individual_id"])
    point_size = 11 if spec.level == "call" else 30

    figure = plt.figure(figsize=(16, 9))
    grid = figure.add_gridspec(2, 3, width_ratios=(1.18, 1, 1), wspace=0.22, hspace=0.27)
    pooled_axis = figure.add_subplot(grid[:, 0])
    facet_axes = (
        figure.add_subplot(grid[0, 1]),
        figure.add_subplot(grid[0, 2]),
        figure.add_subplot(grid[1, 1]),
        figure.add_subplot(grid[1, 2]),
    )

    _plot_embedding_points(
        pooled_axis,
        table,
        x_name,
        y_name,
        colors,
        point_size=point_size,
        show_summaries=True,
    )
    pooled_axis.set_title(f"A   {spec.title}", loc="left", fontweight="bold", fontsize=17)
    pooled_axis.set_xlabel("Cached embedding dimension 1")
    pooled_axis.set_ylabel("Cached embedding dimension 2")

    facets = (
        ("partner", "before", "Partner · Before"),
        ("partner", "after", "Partner · After"),
        ("non-partner", "before", "Non-partner · Before"),
        ("non-partner", "after", "Non-partner · After"),
    )
    for axis, (context, stage, title) in zip(facet_axes, facets):
        subset = table.loc[(table["context"] == context) & (table["stage"] == stage)]
        _plot_embedding_points(
            axis,
            subset,
            x_name,
            y_name,
            colors,
            point_size=point_size,
            show_summaries=True,
        )
        axis.set_title(title)
        axis.set_xlabel("Cached dimension 1" if stage == "after" else "")
        axis.set_ylabel("Cached dimension 2" if context == "non-partner" else "")

    for axis in (pooled_axis, *facet_axes):
        axis.set_xlim(x_limits)
        axis.set_ylim(y_limits)
        axis.grid(color="#dddddd", linewidth=0.5, alpha=0.5)
    facet_axes[0].text(
        -0.22,
        1.15,
        "B   Spread by stage and context",
        transform=facet_axes[0].transAxes,
        fontsize=17,
        fontweight="bold",
    )

    individual_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=colors[individual],
            markeredgecolor="black",
            label=individual,
            markersize=8,
        )
        for individual in sorted(colors)
    ]
    semantic_handles = [
        Line2D([0], [0], marker="^", color="none", markerfacecolor="#888888", label="Partner"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#888888", label="Non-partner"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#aaaaaa", markeredgecolor="black", label="Centroid"),
        Line2D([0, 1], [0, 0], color="#555555", label="Bonded-pair link"),
        Patch(
            facecolor="#aaaaaa",
            edgecolor="#777777",
            alpha=0.2,
            label="50% Gaussian covariance ellipse",
        ),
    ]
    figure.legend(
        handles=[*individual_handles, *semantic_handles],
        loc="lower center",
        ncol=6,
        frameon=False,
        bbox_to_anchor=(0.5, 0.075),
    )
    figure.text(0.5, 0.025, spec.annotation_note, ha="center", fontsize=9, color="#444444")
    figure.subplots_adjust(bottom=0.18, top=0.91)
    return figure


def regenerate_embedding_space(
    embedding: pd.DataFrame,
    metadata: pd.DataFrame,
    spec: EmbeddingFigureSpec,
    locked_path: PathLike,
    *,
    x_column: Optional[str] = None,
    y_column: Optional[str] = None,
    overwrite_regenerated: bool = False,
) -> Path:
    figure = plot_embedding_space(
        embedding,
        metadata,
        spec,
        x_column=x_column,
        y_column=y_column,
    )
    try:
        return save_regenerated_figure(
            figure,
            locked_path,
            overwrite_regenerated=overwrite_regenerated,
        )
    finally:
        plt.close(figure)
