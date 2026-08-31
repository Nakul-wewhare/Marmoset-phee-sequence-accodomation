"""Render the final Script 3 sequence-space distribution figures.

This focused helper contains only the plotting logic needed by
``script_3_sequence_spaces.ipynb``. It reads the fixed 107-repertoire
coordinates and figure selections bundled with the repository; it never fits
an embedding or recalculates a sequence-distance matrix.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field, replace
from itertools import product
import json
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Ellipse, Rectangle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


SUBJECT_ID_CROSSWALK = {
    "Tabor": {"previous_subject_id": "A", "new_subject_id": "Male A"},
    "Lola": {"previous_subject_id": "D", "new_subject_id": "Female A"},
    "Wuschel": {"previous_subject_id": "B", "new_subject_id": "Male B"},
    "Olympia": {"previous_subject_id": "F", "new_subject_id": "Female B"},
    "Odin": {"previous_subject_id": "C", "new_subject_id": "Male C"},
    "Nougatti": {"previous_subject_id": "E", "new_subject_id": "Female C"},
}
INDIVIDUAL_LEGEND_ORDER = (
    "Tabor", "Wuschel", "Odin", "Lola", "Olympia", "Nougatti"
)
INDIVIDUAL_DISPLAY_LABELS = {
    name: identifiers["new_subject_id"]
    for name, identifiers in SUBJECT_ID_CROSSWALK.items()
}
CALL_TOKEN_FULL_NAMES = {
    "A": "Phee", "K": "Egg", "P": "Peep",
    "R": "Trillphee", "S": "SD-peep", "T": "Tsk",
}
CALL_TOKEN_DISPLAY_LABELS = {
    token: f"{token} – {full_name}"
    for token, full_name in CALL_TOKEN_FULL_NAMES.items()
}
PAIR_MEMBERS = {
    "Tabor-Lola": ("Tabor", "Lola"),
    "Odin-Nougatti": ("Odin", "Nougatti"),
    "Wuschel-Olympia": ("Wuschel", "Olympia"),
}
INDIVIDUAL_ORDER = ("Tabor", "Lola", "Odin", "Nougatti", "Wuschel", "Olympia")


@dataclass(frozen=True)
class UmapConfig:
    n_neighbors: int = 30
    min_dist: float = 0.10
    random_state: int = 42


@dataclass(frozen=True)
class PacmapConfig:
    n_neighbors: int = 30
    mn_ratio: float = 0.5
    fp_ratio: float = 2.0
    random_state: int = 42
    distance: str = "euclidean"


FINAL_INDIVIDUAL_COLORS = {
    "Tabor": "#1F78B4",
    "Lola": "#8EC7E8",
    "Odin": "#D95F02",
    "Nougatti": "#FDB863",
    "Wuschel": "#6F3B2E",
    "Olympia": "#C49A87",
}
FINAL_PAIR_LINE_COLORS = {
    "Tabor-Lola": "#1F78B4",
    "Odin-Nougatti": "#D95F02",
    "Wuschel-Olympia": "#6F3B2E",
}
CALL_TOKEN_COLORS = {
    "A": "#0072B2", "K": "#E69F00", "P": "#009E73",
    "R": "#CC79A7", "S": "#D55E00", "T": "#F0E442",
}
CALL_TOKEN_TEXT_COLORS = {
    "A": "white", "K": "#111111", "P": "white",
    "R": "#111111", "S": "white", "T": "#111111",
}
METRIC_LABELS = {
    "local_alignment": "Local alignment",
    "bigram": "Bigram frequency",
    "repeat_A": "Phee-run distribution",
    "transition_probability": "Transition probability",
}
PANEL_ORDER = (
    ("partner", "before"), ("partner", "after"),
    ("non-partner", "before"), ("non-partner", "after"),
)
PANEL_TITLES = {
    ("partner", "before"): "Partner · Before",
    ("partner", "after"): "Partner · After",
    ("non-partner", "before"): "Non-partner · Before",
    ("non-partner", "after"): "Non-partner · After",
}
CONTEXT_MARKERS = {"partner": "^", "non-partner": "s"}
SEQUENCE_ALPHABET = tuple(CALL_TOKEN_COLORS)
UMAP_SEQUENCE_METRICS = frozenset({"local_alignment", "repeat_A"})
PANEL_B_POSITIONS = {
    ("partner", "before"): (0.518, 0.583, 0.218, 0.292),
    ("partner", "after"): (0.752, 0.583, 0.218, 0.292),
    ("non-partner", "before"): (0.518, 0.205, 0.218, 0.292),
    ("non-partner", "after"): (0.752, 0.205, 0.218, 0.292),
}
PANEL_B_SHARED_YLABEL_Y = 0.540
PANEL_B_SHARED_YLABEL_GAP_POINTS = 0.50
MANUSCRIPT_PLACEMENT_WIDTH_INCHES = 7.0


def transition_matrix(sequences: Sequence[str], alphabet: Sequence[str]) -> np.ndarray:
    matrix = pd.DataFrame(0.0, index=alphabet, columns=alphabet)
    for sequence in sequences:
        for first, second in zip(sequence[:-1], sequence[1:]):
            matrix.at[first, second] += 1
    row_totals = matrix.sum(axis=1).replace(0, np.nan)
    return np.nan_to_num(matrix.div(row_totals, axis=0).to_numpy().ravel())


def ngram_distribution(
    sequences: Sequence[str], n: int, all_ngrams: Sequence[tuple[str, ...]]
) -> np.ndarray:
    counts = {ngram: 0 for ngram in all_ngrams}
    for sequence in sequences:
        if len(sequence) < n:
            continue
        for index in range(len(sequence) - n + 1):
            counts[tuple(sequence[index:index + n])] += 1
    total = sum(counts.values()) + 1e-3
    return np.asarray([counts[ngram] / total for ngram in all_ngrams])


def repeat_distribution(
    sequences: Sequence[str], element: str, max_length: int
) -> np.ndarray:
    counts = {length: 0 for length in range(2, max_length + 1)}
    for sequence in sequences:
        index = 0
        while index < len(sequence):
            run_length = 1
            while index + 1 < len(sequence) and sequence[index] == sequence[index + 1]:
                index += 1
                run_length += 1
            if 1 < run_length <= max_length and sequence[index] == element:
                counts[run_length] += 1
            index += 1
    total = sum(counts.values()) + 1e-3
    return np.asarray([counts[length] / total for length in sorted(counts)])


def _bigram_features(
    repertoires: Sequence[Sequence[str]], alphabet: Sequence[str]
) -> tuple[np.ndarray, list[str]]:
    ordered = tuple(product(alphabet, repeat=2))
    features = np.vstack([
        ngram_distribution(repertoire, 2, ordered) for repertoire in repertoires
    ])
    return features, ["".join(bigram) for bigram in ordered]


def _repeat_features(
    repertoires: Sequence[Sequence[str]], target: str = "A", maximum_length: int = 5
) -> np.ndarray:
    return np.vstack([
        repeat_distribution(repertoire, target, maximum_length)
        for repertoire in repertoires
    ])


def _transition_features(
    repertoires: Sequence[Sequence[str]], alphabet: Sequence[str]
) -> np.ndarray:
    return np.vstack([
        transition_matrix(repertoire, alphabet) for repertoire in repertoires
    ])
def _embedding_method_label(metric: str) -> str:
    """Return the final manuscript projection method for one sequence metric."""
    if metric not in METRIC_LABELS:
        raise ValueError(f"Unknown sequence metric: {metric!r}")
    return "UMAP" if metric in UMAP_SEQUENCE_METRICS else "PaCMAP"

def _add_call_type_legend(figure: Figure, *, text_scale: float) -> None:
    """Show each call token and full name without a colour-key dependency."""
    handles = [
        Line2D(
            [],
            [],
            marker=None,
            linestyle="None",
            color="none",
            label=CALL_TOKEN_DISPLAY_LABELS[token],
        )
        for token in SEQUENCE_ALPHABET
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.245, 0.895),
        ncol=len(handles),
        frameon=False,
        # Six full call-name labels must remain on one line above Panel A.
        # This role is intentionally denser than the general figure text.
        fontsize=8.20 * float(text_scale),
        columnspacing=0.78,
        handlelength=0.0,
        handletextpad=0.0,
    )

@dataclass(frozen=True)
class FinalSequenceFigureConfig:
    """Projection, example-selection, and plotting settings."""

    umap: UmapConfig = field(default_factory=UmapConfig)
    pacmap: PacmapConfig = field(default_factory=PacmapConfig)
    figure_size_inches: tuple[float, float] = (16.0, 9.5)
    dpi: int = 300
    # The landscape canvas is reduced to journal page width. A 1.68 scale is
    # calibrated for a full-width, approximately half-page placement and is
    # the single typography control exposed in script_3_sequence_spaces.ipynb.
    text_scale: float = 1.68
    examples_per_side: int = 5
    heatmap_examples_per_side: int | None = None
    repeat_examples_per_side: int | None = None
    point_size: float = 82.0
    point_alpha: float = 0.97
    centroid_size: float = 235.0
    example_highlight_size: float = 130.0
    distribution_ellipse_radius_scale: float = 1.0
    connector_linewidth: float = 0.90
    connector_alpha: float = 0.30

    @property
    def example_count(self) -> int:
        return 4 * int(self.examples_per_side)

    def examples_per_side_for_metric(self, metric: str) -> int:
        """Return the border density requested for one metric."""
        if metric in {"bigram", "transition_probability"}:
            value = self.heatmap_examples_per_side
            return int(self.examples_per_side if value is None else value)
        if metric == "repeat_A":
            value = self.repeat_examples_per_side
            return int(self.examples_per_side if value is None else value)
        return int(self.examples_per_side)

def _font_size(config: FinalSequenceFigureConfig, base_points: float) -> float:
    """Scale manuscript typography from one notebook-facing control."""
    if config.text_scale <= 0:
        raise ValueError("text_scale must be positive.")
    return float(base_points) * float(config.text_scale)

def _ordered_sessions(sessions: pd.DataFrame) -> pd.DataFrame:
    required = {
        "group_id",
        "focal ID",
        "conspecific_ID",
        "stage",
        "paired_status",
        "pair_id",
        "sex",
        "sequences_json",
    }
    missing = sorted(required - set(sessions.columns))
    if missing:
        raise KeyError(f"Sequence-session manifest is missing columns: {missing}")
    ordered = sessions.sort_values("group_id").reset_index(drop=True).copy()
    expected_ids = np.arange(len(ordered), dtype=int)
    if not np.array_equal(ordered["group_id"].to_numpy(dtype=int), expected_ids):
        raise ValueError("Sequence group_id must be contiguous and in matrix order.")
    if "seq_list" not in ordered:
        ordered["seq_list"] = ordered["sequences_json"].map(json.loads)
    else:
        ordered["seq_list"] = ordered["seq_list"].map(
            lambda value: value if isinstance(value, list) else json.loads(value)
        )
    unknown = sorted(set(ordered["focal ID"]) - set(FINAL_INDIVIDUAL_COLORS))
    if unknown:
        raise ValueError(f"No final-figure colour is defined for: {unknown}")
    observed_tokens = {
        token
        for repertoire in ordered["seq_list"]
        for sequence in repertoire
        for token in sequence
    }
    unknown_tokens = sorted(observed_tokens - set(CALL_TOKEN_COLORS))
    if unknown_tokens:
        raise ValueError(f"No call-token colour is defined for: {unknown_tokens}")
    return ordered

def _embedding_limits(
    embedding: np.ndarray,
    *,
    padding_fraction: float = 0.065,
) -> tuple[tuple[float, float], tuple[float, float]]:
    if padding_fraction < 0.0:
        raise ValueError("padding_fraction must be non-negative.")
    minimum = np.asarray(embedding).min(axis=0)
    maximum = np.asarray(embedding).max(axis=0)
    span = np.maximum(maximum - minimum, 1e-6)
    padding = span * float(padding_fraction)
    return (
        (float(minimum[0] - padding[0]), float(maximum[0] + padding[0])),
        (float(minimum[1] - padding[1]), float(maximum[1] + padding[1])),
    )

def _normalise_embedding(embedding: np.ndarray) -> np.ndarray:
    minimum = np.asarray(embedding).min(axis=0)
    span = np.maximum(np.asarray(embedding).max(axis=0) - minimum, 1e-9)
    return (np.asarray(embedding) - minimum) / span

def _border_slot_definitions(examples_per_side: int) -> list[tuple[str, float, float]]:
    if examples_per_side < 2:
        raise ValueError("examples_per_side must be at least two.")
    spread = np.linspace(0.08, 0.92, examples_per_side)
    return (
        [("top", float(value), 1.10) for value in spread]
        + [("right", 1.10, float(value)) for value in spread[::-1]]
        + [("bottom", float(value), -0.10) for value in spread[::-1]]
        + [("left", -0.10, float(value)) for value in spread]
    )

def _heatmap_corner_slots(examples_per_side: int) -> dict[int, str]:
    """Map the four outer heatmaps to outward-facing axis locations."""
    count = int(examples_per_side)
    return {
        0: "top_left",
        count - 1: "top_right",
        2 * count: "bottom_right",
        3 * count - 1: "bottom_left",
    }

def _assign_border_slots(
    selected: pd.DataFrame,
    embedding: np.ndarray,
    examples_per_side: int,
    metric: str,
) -> pd.DataFrame:
    slots = _border_slot_definitions(examples_per_side)
    if len(selected) != len(slots):
        raise ValueError("Selected repertoire count does not match border slots.")
    points = _normalise_embedding(embedding)[selected["group_id"].to_numpy(dtype=int)]
    slot_coordinates = np.asarray([(x, y) for _, x, y in slots], dtype=float)
    costs = np.square(points[:, None, :] - slot_coordinates[None, :, :]).sum(axis=2)
    # Match card geometry to border geometry: unusually wide cards are mildly
    # discouraged on the horizontal sides, while unusually tall cards are
    # discouraged on the vertical sides. Geometric proximity still dominates.
    layouts = [
        _metric_card_layout(row, metric, examples_per_side)
        for _, row in selected.iterrows()
    ]
    widths = np.asarray([layout.width_inches for layout in layouts], dtype=float)
    heights = np.asarray([layout.height_inches for layout in layouts], dtype=float)

    def normalised(values: np.ndarray) -> np.ndarray:
        span = float(np.ptp(values))
        return np.zeros_like(values) if span == 0 else (values - values.min()) / span

    horizontal_slot = np.asarray(
        [side in {"top", "bottom"} for side, _, _ in slots], dtype=float
    )
    vertical_slot = 1.0 - horizontal_slot
    costs += 0.075 * normalised(widths)[:, None] * horizontal_slot[None, :]
    costs += 0.075 * normalised(heights)[:, None] * vertical_slot[None, :]
    selected_indices, slot_indices = linear_sum_assignment(costs)
    slot_for_selected = np.empty(len(selected), dtype=int)
    slot_for_selected[selected_indices] = slot_indices
    assigned = selected.copy()
    assigned["border_slot"] = slot_for_selected
    assigned["border_order"] = assigned["border_slot"] + 1
    assigned["border_side"] = [slots[index][0] for index in slot_for_selected]
    return assigned.sort_values("border_slot").reset_index(drop=True)

def select_final_repertoire_examples(
    sessions: pd.DataFrame,
    embedding: np.ndarray,
    *,
    examples_per_side: int = 4,
    example_group_ids: Sequence[int] | None = None,
    metric: str = "local_alignment",
) -> pd.DataFrame:
    """Select balanced, space-covering repertoire examples for Panel A.

    The first twelve examples are the centroid-nearest repertoire for every
    individual x stage combination.  Remaining slots are filled by deterministic
    farthest-point sampling, with at most one extra repertoire per individual
    until all individuals have received that extra slot.
    """
    ordered = _ordered_sessions(sessions)
    embedding = np.asarray(embedding, dtype=np.float32)
    if embedding.shape != (len(ordered), 2) or not np.isfinite(embedding).all():
        raise ValueError("Embedding must be finite with shape (number of sessions, 2).")
    target = 4 * int(examples_per_side)
    if target > len(ordered):
        raise ValueError("More border examples were requested than session repertoires.")
    if example_group_ids is not None:
        group_ids = [int(value) for value in example_group_ids]
        if len(group_ids) != target or len(set(group_ids)) != target:
            raise ValueError(f"example_group_ids must contain {target} unique IDs.")
        unknown = sorted(set(group_ids) - set(ordered["group_id"].astype(int)))
        if unknown:
            raise ValueError(f"Unknown sequence group IDs: {unknown}")
        selected = ordered.set_index("group_id", drop=False).loc[group_ids].reset_index(drop=True)
        selected["selection_reason"] = "manual"
        return _assign_border_slots(selected, embedding, examples_per_side, metric)

    if target < 2 * len(INDIVIDUAL_ORDER):
        raise ValueError(
            "At least 12 examples are needed to cover every individual x stage."
        )
    selected_ids: list[int] = []
    selection_reasons: dict[int, str] = {}
    for individual in INDIVIDUAL_ORDER:
        for stage in ("before", "after"):
            group = ordered.loc[
                ordered["focal ID"].eq(individual) & ordered["stage"].eq(stage)
            ]
            if group.empty:
                raise ValueError(f"No repertoire for {individual} x {stage}.")
            group_ids = group["group_id"].to_numpy(dtype=int)
            group_embedding = embedding[group_ids]
            centroid = group_embedding.mean(axis=0)
            chosen = int(group_ids[np.argmin(np.square(group_embedding - centroid).sum(axis=1))])
            selected_ids.append(chosen)
            selection_reasons[chosen] = "individual-stage centroid-nearest"

    normalised = _normalise_embedding(embedding)
    extra_counts = Counter(ordered.loc[selected_ids, "focal ID"].astype(str))
    while len(selected_ids) < target:
        remaining = [
            group_id
            for group_id in ordered["group_id"].astype(int)
            if group_id not in selected_ids
        ]
        minimum_selected_count = min(extra_counts.values())
        balanced = [
            group_id
            for group_id in remaining
            if extra_counts[str(ordered.iloc[group_id]["focal ID"])] == minimum_selected_count
        ]
        candidates = balanced or remaining
        chosen = max(
            candidates,
            key=lambda group_id: (
                float(
                    np.min(
                        np.square(normalised[selected_ids] - normalised[group_id]).sum(axis=1)
                    )
                ),
                -int(ordered.iloc[group_id]["n_unique_sequences"]),
                -group_id,
            ),
        )
        selected_ids.append(int(chosen))
        individual = str(ordered.iloc[chosen]["focal ID"])
        extra_counts[individual] += 1
        selection_reasons[int(chosen)] = "balanced farthest-point coverage"

    selected = ordered.set_index("group_id", drop=False).loc[selected_ids].reset_index(drop=True)
    selected["selection_reason"] = selected["group_id"].map(selection_reasons)
    return _assign_border_slots(selected, embedding, examples_per_side, metric)

def _style_embedding_axis(
    axis,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
    config: FinalSequenceFigureConfig,
) -> None:
    axis.set_xlim(*x_limits)
    axis.set_ylim(*y_limits)
    axis.grid(False)
    axis.set_facecolor("white")
    axis.tick_params(
        labelsize=_font_size(config, 10), length=3.4, width=1.0, pad=2.0
    )
    for spine in axis.spines.values():
        spine.set_color("#555555")
        spine.set_linewidth(1.0)

def _ordered_sequence_counts(seq_list: Sequence[str]) -> list[tuple[str, int]]:
    return sorted(
        Counter(seq_list).items(),
        key=lambda item: (-item[1], len(item[0]), item[0]),
    )

@dataclass(frozen=True)
class _RepertoireCardLayout:
    """Content-sized repertoire-card geometry in physical inches."""

    columns: tuple[tuple[tuple[str, int], ...], ...]
    column_widths_inches: tuple[float, ...]
    width_inches: float
    height_inches: float
    block_width_inches: float = 0.170
    block_height_inches: float = 0.160
    row_pitch_inches: float = 0.190
    column_gap_inches: float = 0.18
    side_margin_inches: float = 0.10
    header_inches: float = 0.35
    token_font_base: float = 7.0
    count_font_base: float = 6.8
    title_font_base: float = 9.2

def _repertoire_card_layout(
    row: pd.Series,
    *,
    compact_side: bool = False,
) -> _RepertoireCardLayout:
    """Size a card from its exact content rather than a fixed thumbnail box."""
    ordered_counts = _ordered_sequence_counts(row["seq_list"])
    row_count = len(ordered_counts)
    if compact_side:
        # Vertical border cards have only ~1.4 physical inches available. Use
        # one narrower column so their complete repertoire stays on-canvas.
        block_width_inches = 0.140
        count_allowance = 0.24 if any(count > 1 for _, count in ordered_counts) else 0.04
        column_width = (
            max(len(sequence) for sequence, _ in ordered_counts)
            * block_width_inches
            + count_allowance
        )
        return _RepertoireCardLayout(
            columns=(tuple(ordered_counts),),
            column_widths_inches=(column_width,),
            width_inches=max(0.68, 0.15 + column_width),
            height_inches=max(0.70, 0.36 + row_count * 0.145),
            block_width_inches=block_width_inches,
            block_height_inches=0.125,
            row_pitch_inches=0.145,
            column_gap_inches=0.0,
            side_margin_inches=0.075,
            header_inches=0.32,
            token_font_base=6.6,
            count_font_base=6.4,
            title_font_base=8.4,
        )
    # Dense repertoires always use two columns. This preserves large blocks and
    # prevents a many-type repertoire from becoming an intrusive tall strip.
    column_count = 2 if row_count >= 8 else 1
    rows_per_column = int(np.ceil(row_count / column_count))
    columns = tuple(
        tuple(
            ordered_counts[
                column_index * rows_per_column : (column_index + 1) * rows_per_column
            ]
        )
        for column_index in range(column_count)
    )
    block_width_inches = 0.170
    column_widths = tuple(
        max(len(sequence) for sequence, _ in column) * block_width_inches
        + (0.34 if any(count > 1 for _, count in column) else 0.05)
        for column in columns
    )
    width_inches = max(
        0.78,
        0.20 + sum(column_widths) + 0.18 * (column_count - 1),
    )
    height_inches = max(0.78, 0.42 + rows_per_column * 0.190)
    return _RepertoireCardLayout(
        columns=columns,
        column_widths_inches=column_widths,
        width_inches=width_inches,
        height_inches=height_inches,
    )

def _metric_card_layout(
    row: pd.Series,
    metric: str,
    examples_per_side: int = 6,
    border_side: str | None = None,
) -> _RepertoireCardLayout:
    """Return content-aware sequence cards or compact metric-summary cards."""
    if metric == "local_alignment":
        layout = _repertoire_card_layout(
            row,
            compact_side=border_side in {"left", "right"},
        )
        # Local-alignment cards contain the full exact repertoire and otherwise
        # dominate Panel A. A uniform 10% geometry reduction preserves all
        # sequence content while giving the embedding the visual priority.
        scale = 0.90
        return replace(
            layout,
            column_widths_inches=tuple(
                width * scale for width in layout.column_widths_inches
            ),
            width_inches=layout.width_inches * scale,
            height_inches=layout.height_inches * scale,
            block_width_inches=layout.block_width_inches * scale,
            block_height_inches=layout.block_height_inches * scale,
            row_pitch_inches=layout.row_pitch_inches * scale,
            column_gap_inches=layout.column_gap_inches * scale,
            side_margin_inches=layout.side_margin_inches * scale,
            header_inches=layout.header_inches * scale,
        )
    if metric in {"bigram", "transition_probability"}:
        width = (
            1.02
            if examples_per_side <= 4
            else (0.85 if examples_per_side == 5 else 0.70)
        )
        height = min(width, 0.92) if border_side in {"left", "right"} else width
        return _RepertoireCardLayout(
            columns=(),
            column_widths_inches=(),
            width_inches=width,
            height_inches=height,
        )
    if metric == "repeat_A":
        if examples_per_side <= 4:
            # Keep the four axis-bearing horizontal cards large enough for
            # readable scales, while narrower side cards preserve trim-safe
            # outer margins and breathing room before Panel B.
            if border_side in {"left", "right"}:
                width, height = 1.18, 0.95
            else:
                width, height = 1.22, 1.05
        else:
            width = 0.92 if examples_per_side == 5 else 0.76
            height = 0.72 if examples_per_side == 5 else 0.62
        return _RepertoireCardLayout(
            columns=(),
            column_widths_inches=(),
            width_inches=width,
            height_inches=height,
        )
    raise ValueError(f"Unknown sequence metric: {metric!r}")

def _packed_card_positions(
    representatives: pd.DataFrame,
    config: FinalSequenceFigureConfig,
    metric: str = "local_alignment",
) -> list[tuple[tuple[float, float, float, float], str, _RepertoireCardLayout]]:
    """Pack variable-sized cards on the four sides without overlap."""
    examples_per_side = config.examples_per_side_for_metric(metric)
    allowed = {3, 4, 5} if metric == "local_alignment" else {3, 4, 5, 6}
    if examples_per_side not in allowed:
        supported = "3, 4, or 5" if metric == "local_alignment" else "3 through 6"
        raise ValueError(
            f"The manuscript layout supports {supported} examples per side for {metric}."
        )
    if len(representatives) != 4 * examples_per_side:
        raise ValueError(
            f"Expected {4 * examples_per_side} linked examples for {metric}, "
            f"received {len(representatives)}."
        )
    figure_width, figure_height = config.figure_size_inches
    records = []
    for _, row in representatives.iterrows():
        side = str(row["border_side"])
        records.append(
            {
                "slot": int(row["border_slot"]),
                "side": side,
                "layout": _metric_card_layout(
                    row,
                    metric,
                    examples_per_side,
                    border_side=side,
                ),
            }
        )
    positioned = {}
    horizontal_min, horizontal_max = 0.015, 0.485
    # Reserve a readable central plot and a bottom legend band. Fewer, larger
    # example cards are aligned to a common outer edge at publication size.
    if metric == "local_alignment":
        top_edge, bottom_y = 0.870, 0.015
        left_edge, right_edge = 0.095, 0.395
    elif metric in {"bigram", "transition_probability"}:
        top_edge, bottom_y = 0.850, 0.065
        left_edge, right_edge = 0.078, 0.402
    else:
        top_edge = 0.850
        bottom_y = 0.100 if metric == "repeat_A" else 0.065
        left_edge, right_edge = 0.085, 0.400

    for side in ("top", "bottom"):
        items = [record for record in records if record["side"] == side]
        widths = [record["layout"].width_inches / figure_width for record in items]
        available = horizontal_max - horizontal_min
        if metric == "local_alignment":
            gap = (available - sum(widths)) / max(len(items) - 1, 1)
            cursor = horizontal_min
        else:
            # Fixed-size metric summaries read as one compact row rather than
            # isolated thumbnails spread across the complete Panel A width.
            gap = 0.020 if metric == "repeat_A" else 0.006
            occupied = sum(widths) + gap * max(len(items) - 1, 0)
            cursor = horizontal_min + (available - occupied) / 2.0
        if gap < 0.005:
            raise ValueError(
                f"The {side} repertoire cards do not fit; reduce examples_per_side."
            )
        ordered_items = items if side == "top" else list(reversed(items))
        for record in ordered_items:
            layout = record["layout"]
            width = layout.width_inches / figure_width
            height = layout.height_inches / figure_height
            y = top_edge - height if side == "top" else bottom_y
            positioned[record["slot"]] = ((cursor, y, width, height), side, layout)
            cursor += width + gap

    # Fit the vertical rows exactly between the tallest top and bottom cards.
    # This keeps the four corner cards separated even when local-alignment
    # repertoires have content-dependent dimensions.
    top_items = [record for record in records if record["side"] == "top"]
    bottom_items = [record for record in records if record["side"] == "bottom"]
    top_bottom = min(
        positioned[record["slot"]][0][1] for record in top_items
    )
    bottom_top = max(
        positioned[record["slot"]][0][1]
        + positioned[record["slot"]][0][3]
        for record in bottom_items
    )
    corner_gap = (
        0.030
        if metric == "repeat_A"
        else (0.009 if metric in {"bigram", "transition_probability"} else 0.005)
    )
    vertical_min = bottom_top + corner_gap
    vertical_max = top_bottom - corner_gap

    for side in ("right", "left"):
        items = [record for record in records if record["side"] == side]
        heights = [record["layout"].height_inches / figure_height for record in items]
        available = vertical_max - vertical_min
        gap = (available - sum(heights)) / max(len(items) - 1, 1)
        if gap < 0.005:
            raise ValueError(
                f"The {side} repertoire cards do not fit; reduce examples_per_side."
            )
        if side == "right":
            cursor = vertical_max
            for record in items:  # right-side slots are ordered top to bottom
                layout = record["layout"]
                width = layout.width_inches / figure_width
                height = layout.height_inches / figure_height
                y = cursor - height
                positioned[record["slot"]] = ((right_edge, y, width, height), side, layout)
                cursor = y - gap
        else:
            cursor = vertical_min
            for record in items:  # left-side slots are ordered bottom to top
                layout = record["layout"]
                width = layout.width_inches / figure_width
                height = layout.height_inches / figure_height
                positioned[record["slot"]] = (
                    (left_edge - width, cursor, width, height),
                    side,
                    layout,
                )
                cursor += height + gap

    return [positioned[int(row["border_slot"])] for _, row in representatives.iterrows()]

def _draw_repertoire_card(
    axis,
    row: pd.Series,
    layout: _RepertoireCardLayout,
    config: FinalSequenceFigureConfig,
) -> None:
    width = layout.width_inches
    height = layout.height_inches
    x_cursor_inches = layout.side_margin_inches
    for column_index, column in enumerate(layout.columns):
        for row_index, (sequence, count) in enumerate(column):
            center_y_inches = (
                height
                - layout.header_inches
                - row_index * layout.row_pitch_inches
                - layout.block_height_inches / 2.0
            )
            center_y = center_y_inches / height
            lower_y = (
                center_y_inches - layout.block_height_inches / 2.0
            ) / height
            start_x = x_cursor_inches / width
            block_width = layout.block_width_inches / width
            block_height = layout.block_height_inches / height
            for token_index, token in enumerate(sequence):
                left_x = start_x + token_index * block_width
                axis.add_patch(
                    Rectangle(
                        (left_x, lower_y),
                        block_width * 0.93,
                        block_height,
                        facecolor=CALL_TOKEN_COLORS[token],
                        edgecolor="#333333",
                        linewidth=0.35,
                    )
                )
                axis.text(
                    left_x + block_width * 0.465,
                    center_y,
                    token,
                    ha="center",
                    va="center",
                    fontsize=_font_size(config, layout.token_font_base),
                    color=CALL_TOKEN_TEXT_COLORS[token],
                    fontweight="bold",
                )
            if count > 1:
                axis.text(
                    start_x + len(sequence) * block_width + 0.035 / width,
                    center_y,
                    f"×{count}",
                    ha="left",
                    va="center",
                    fontsize=_font_size(config, layout.count_font_base),
                    color="#222222",
                )
        x_cursor_inches += layout.column_widths_inches[column_index]
        if column_index < len(layout.columns) - 1:
            x_cursor_inches += layout.column_gap_inches

    axis.text(
        0.5,
        1.0 - 0.045 / height,
        _example_card_title(row, compact=True),
        ha="center",
        va="top",
        fontsize=_font_size(config, layout.title_font_base),
        color="#222222",
    )
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_facecolor("white")
    color = FINAL_INDIVIDUAL_COLORS[str(row["focal ID"])]
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color(color)
        spine.set_linewidth(1.7)

def _metric_card_heatmap_maximum(sessions: pd.DataFrame, metric: str) -> float:
    """Return the shared colour ceiling used by every example heatmap."""
    if metric == "bigram":
        features, _ = _bigram_features(sessions["seq_list"].tolist(), SEQUENCE_ALPHABET)
        return max(float(features.max()), 1e-12)
    if metric == "transition_probability":
        return 1.0
    return 1.0

def _example_card_title(row: pd.Series, *, compact: bool = False) -> str:
    """Describe both experimental stage and repertoire social context."""
    stage = str(row["stage"]).title()
    context = str(row["paired_status"]).replace("_", " ").lower()
    if compact:
        abbreviation = "P" if context == "partner" else "NP"
        return f"{stage} · {abbreviation}"
    return f"{stage}\n{context.title()}"

def _add_metric_card_colorbar(
    figure: Figure,
    metric: str,
    maximum: float,
    config: FinalSequenceFigureConfig,
) -> None:
    """Add the one shared scale used by all linked example heatmaps."""
    if metric not in {"bigram", "transition_probability"}:
        return
    # Use the protected gap between the bottom example row and the PaCMAP.
    # Placing the scale below the cards collides with the two corner axes once
    # the figure is reduced to journal width.
    # An opaque backing also keeps the linked-example connectors from running
    # through the scale or its labels.
    colorbar_backing = Rectangle(
        (0.145, 0.182),
        0.200,
        0.066,
        transform=figure.transFigure,
        facecolor="white",
        edgecolor="none",
        clip_on=False,
        zorder=18,
    )
    figure.add_artist(colorbar_backing)
    colorbar_axis = figure.add_axes((0.160, 0.207, 0.170, 0.014))
    colorbar_axis.set_zorder(19)
    colorbar_axis.set_facecolor("white")
    mappable = ScalarMappable(
        norm=Normalize(vmin=0.0, vmax=float(maximum)),
        cmap="viridis",
    )
    mappable.set_array([])
    colorbar = figure.colorbar(
        mappable,
        cax=colorbar_axis,
        orientation="horizontal",
    )
    colorbar.set_ticks((0.0, float(maximum) / 2.0, float(maximum)))
    colorbar.ax.tick_params(
        labelsize=_font_size(config, 7.0), length=1.4, pad=0.6
    )
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.set_label(
        "Bigram frequency"
        if metric == "bigram"
        else "Transition probability",
        fontsize=_font_size(config, 7.5),
        labelpad=0.8,
    )

def _format_panel_b_axes(
    movement_axes: Mapping[tuple[str, str], object],
    embedding_label: str,
    config: FinalSequenceFigureConfig,
) -> None:
    """Show a compact numeric scale on every Panel B facet.

    X tick numbers remain below all four panels.  Y tick numbers face outwards:
    left of the Before column and right of the After column.  This keeps the
    narrow inter-column gutter clear and prevents the shared dimension-2 label
    from competing with the right-hand example cards in Panel A.
    """
    for panel_key, axis in movement_axes.items():
        paired_status, stage = panel_key
        axis.tick_params(
            axis="x",
            bottom=True,
            top=False,
            labelbottom=True,
            labeltop=False,
            pad=1.0,
        )

        if stage == "before":
            axis.yaxis.set_ticks_position("left")
            axis.tick_params(
                axis="y",
                left=True,
                right=False,
                labelleft=True,
                labelright=False,
                length=1.5,
                pad=0.0,
            )
        else:
            axis.yaxis.set_ticks_position("right")
            axis.tick_params(
                axis="y",
                left=False,
                right=True,
                labelleft=False,
                labelright=True,
                length=1.5,
                pad=0.0,
            )

        if paired_status == "non-partner":
            axis.set_xlabel(
                f"{embedding_label}-1",
                fontsize=_font_size(config, 10.5),
            )

def _add_panel_b_shared_y_label(
    figure: Figure,
    movement_axes: Mapping[tuple[str, str], object],
    embedding_label: str,
    config: FinalSequenceFigureConfig,
):
    """Place the shared dimension-2 label immediately before the y tick text."""
    label = figure.text(
        0.487,
        PANEL_B_SHARED_YLABEL_Y,
        f"{embedding_label}-2",
        ha="right",
        va="center",
        rotation=90,
        fontsize=_font_size(config, 10.5),
    )
    # Tick-label widths differ between embeddings (for example, ``7`` versus
    # ``-2``).  Positioning from the rendered tick boxes gives every metric the
    # same deliberately tight half-point gap at the intended 7-inch placement.
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    tick_bounds = []
    for panel_key in (("partner", "before"), ("non-partner", "before")):
        axis = movement_axes[panel_key]
        lower, upper = sorted(map(float, axis.get_ylim()))
        for tick in axis.yaxis.majorTicks:
            if (
                lower <= tick.get_loc() <= upper
                and tick.label1.get_visible()
                and tick.label1.get_text().strip()
            ):
                tick_bounds.append(tick.label1.get_window_extent(renderer))
    if not tick_bounds:
        raise RuntimeError("Panel B has no visible left-column y tick labels.")
    tick_left = min(bounds.x0 for bounds in tick_bounds)
    target_gap_fraction = PANEL_B_SHARED_YLABEL_GAP_POINTS / (
        MANUSCRIPT_PLACEMENT_WIDTH_INCHES * 72.0
    )
    label.set_x(
        (tick_left - figure.bbox.x0) / figure.bbox.width - target_gap_fraction
    )
    return label

def _draw_metric_card(
    axis,
    row: pd.Series,
    layout: _RepertoireCardLayout,
    metric: str,
    config: FinalSequenceFigureConfig,
    *,
    heatmap_maximum: float = 1.0,
    border_side: str | None = None,
    scale_corner: str | None = None,
) -> None:
    """Draw the analysis feature represented by one linked example point."""
    if metric == "local_alignment":
        _draw_repertoire_card(axis, row, layout, config)
        return

    repertoire = row["seq_list"]
    individual = str(row["focal ID"])
    border_color = FINAL_INDIVIDUAL_COLORS[individual]
    if metric in {"bigram", "transition_probability"}:
        if metric == "bigram":
            values, _ = _bigram_features([repertoire], SEQUENCE_ALPHABET)
            matrix = values[0].reshape(len(SEQUENCE_ALPHABET), -1)
        else:
            matrix = _transition_features(
                [repertoire], SEQUENCE_ALPHABET
            )[0].reshape(len(SEQUENCE_ALPHABET), -1)
        axis.imshow(
            matrix,
            cmap="viridis",
            vmin=0.0,
            vmax=float(heatmap_maximum),
            interpolation="nearest",
            aspect="auto",
        )
        ticks = np.arange(len(SEQUENCE_ALPHABET))
        if scale_corner is not None:
            axis.set_xticks(ticks, labels=SEQUENCE_ALPHABET)
            axis.set_yticks(ticks, labels=SEQUENCE_ALPHABET)
            axis.tick_params(
                axis="both",
                labelsize=_font_size(config, 6.6),
                length=1.35,
                pad=0.3,
            )
            axis.set_xlabel(
                "Next", fontsize=_font_size(config, 7.0), labelpad=0.35
            )
            axis.set_ylabel(
                "Current", fontsize=_font_size(config, 7.0), labelpad=0.35
            )
            if scale_corner.startswith("top"):
                axis.xaxis.set_ticks_position("top")
                axis.xaxis.set_label_position("top")
                axis.tick_params(axis="x", labeltop=True, labelbottom=False)
            if scale_corner.endswith("right"):
                axis.yaxis.set_ticks_position("right")
                axis.yaxis.set_label_position("right")
                axis.tick_params(axis="y", labelright=True, labelleft=False)
        else:
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_xlabel("")
            axis.set_ylabel("")
        # Reserve a small in-card header so tightly stacked side examples do
        # not collide with one another.
        axis.set_ylim(len(SEQUENCE_ALPHABET) - 0.5, -1.20)
        axis.text(
            (len(SEQUENCE_ALPHABET) - 1) / 2.0,
            -0.76,
            _example_card_title(row, compact=True),
            ha="center",
            va="center",
            fontsize=_font_size(config, 6.6),
            color="#222222",
            clip_on=False,
        )
    elif metric == "repeat_A":
        values = _repeat_features([repertoire])[0]
        run_lengths = np.arange(2, 6)
        axis.plot(
            run_lengths,
            values,
            color="#4A4A4A",
            linewidth=1.15,
            marker="o",
            markersize=3.2,
            markerfacecolor="#4A4A4A",
            markeredgecolor="#4A4A4A",
            zorder=2,
        )
        axis.scatter(
            run_lengths,
            values,
            s=8.0,
            facecolor="white",
            edgecolor="#4A4A4A",
            linewidth=0.65,
            zorder=3,
        )
        axis.set_xlim(1.5, 5.5)
        # Reserve a clean header band above the 0-to-1 profile so the
        # stage/context label never sits on the border or covers the data.
        axis.set_ylim(-0.04, 1.50)
        if scale_corner is not None:
            axis.set_xticks(run_lengths, labels=("2", "3", "4", "5"))
            axis.set_yticks((0.0, 0.5, 1.0), labels=("0", "0.5", "1"))
            axis.tick_params(
                axis="both",
                labelsize=_font_size(config, 7.5),
                length=1.5,
                pad=0.5,
            )
            axis.set_xlabel(
                "Phee-run length", fontsize=_font_size(config, 8.2), labelpad=0.7
            )
            axis.set_ylabel(
                "Rel. freq.",
                fontsize=_font_size(config, 8.2),
                labelpad=0.7,
            )
            if scale_corner.startswith("top"):
                axis.xaxis.set_ticks_position("top")
                axis.xaxis.set_label_position("top")
                axis.tick_params(axis="x", labeltop=True, labelbottom=False)
            if scale_corner.endswith("right"):
                axis.yaxis.set_ticks_position("right")
                axis.yaxis.set_label_position("right")
                axis.tick_params(axis="y", labelright=True, labelleft=False)
        else:
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_xlabel("")
            axis.set_ylabel("")
        axis.text(
            0.5,
            0.91,
            _example_card_title(row, compact=True),
            transform=axis.transAxes,
            ha="center",
            va="top",
            fontsize=_font_size(config, 7.4),
            color="#222222",
        )
    else:
        raise ValueError(f"Unknown sequence metric: {metric!r}")

    axis.set_facecolor("white")
    for spine in axis.spines.values():
        if metric in {"bigram", "transition_probability"}:
            spine.set_visible(False)
        else:
            spine.set_visible(True)
            spine.set_color(border_color)
            spine.set_linewidth(1.7)

def _connect_card_to_point(
    figure: Figure,
    card_axis,
    overview_axis,
    point: np.ndarray,
    individual: str,
    border_side: str,
    config: FinalSequenceFigureConfig,
    marker: str = "o",
) -> None:
    anchor = {
        "top": (0.5, 0.0),
        "right": (0.0, 0.5),
        "bottom": (0.5, 1.0),
        "left": (1.0, 0.5),
    }[border_side]
    color = FINAL_INDIVIDUAL_COLORS[individual]
    figure.add_artist(
        ConnectionPatch(
            xyA=anchor,
            coordsA=card_axis.transAxes,
            xyB=(float(point[0]), float(point[1])),
            coordsB=overview_axis.transData,
            arrowstyle="-",
            color=color,
            linewidth=config.connector_linewidth,
            alpha=config.connector_alpha,
            clip_on=False,
            zorder=5,
        )
    )
    overview_axis.scatter(
        [point[0]],
        [point[1]],
        s=config.example_highlight_size * 1.45,
        marker=marker,
        c="white",
        edgecolors="white",
        linewidths=1.8,
        zorder=7,
    )
    overview_axis.scatter(
        [point[0]],
        [point[1]],
        s=config.example_highlight_size,
        marker=marker,
        c=color,
        edgecolors="#111111",
        linewidths=1.1,
        zorder=8,
    )

def _draw_distribution_overview(
    axis,
    sessions: pd.DataFrame,
    embedding: np.ndarray,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
    config: FinalSequenceFigureConfig,
    embedding_label: str,
) -> None:
    """Draw the complete space with context shapes and no overall centroids."""
    for individual in INDIVIDUAL_ORDER:
        color = FINAL_INDIVIDUAL_COLORS[individual]
        for paired_status in ("partner", "non-partner"):
            indices = np.flatnonzero(
                sessions["focal ID"].eq(individual).to_numpy()
                & sessions["paired_status"].eq(paired_status).to_numpy()
            )
            if not len(indices):
                continue
            axis.scatter(
                embedding[indices, 0],
                embedding[indices, 1],
                s=config.point_size * 1.12,
                marker=CONTEXT_MARKERS[paired_status],
                c=color,
                alpha=config.point_alpha,
                linewidths=0,
                zorder=2,
            )
    _style_embedding_axis(axis, x_limits, y_limits, config)
    axis.set_xlabel(f"{embedding_label}-1", fontsize=_font_size(config, 10))
    axis.set_ylabel(f"{embedding_label}-2", fontsize=_font_size(config, 10))

def _distribution_ellipse_parameters(
    points: np.ndarray,
    embedding: np.ndarray,
    *,
    confidence: float = 0.50,
    minimum_spread_fraction: float = 0.006,
    radius_scale: float = 1.0,
) -> tuple[np.ndarray, float, float, float] | None:
    """Return a regularized empirical coverage ellipse for two or more points.

    Small subsets have very unstable sample covariances (two points otherwise
    produce an arbitrarily thin, very long ellipse).  An n-dependent isotropic
    shrinkage and a small variance floor keep those descriptive regions legible.
    These are display regularizers, not additional observations.
    """
    points = np.asarray(points, dtype=float)
    if len(points) < 2:
        return None
    centroid = points.mean(axis=0)
    covariance = np.asarray(np.cov(points, rowvar=False, ddof=0), dtype=float)
    if covariance.shape != (2, 2) or not np.isfinite(covariance).all():
        return None
    isotropic_target = np.eye(2) * float(np.trace(covariance) / 2.0)
    shrinkage = min(0.70, 2.0 / (len(points) + 1.0))
    covariance = (
        (1.0 - shrinkage) * covariance
        + shrinkage * isotropic_target
    )
    embedding_span = np.maximum(np.ptp(np.asarray(embedding), axis=0), 1e-9)
    floor = np.square(embedding_span * float(minimum_spread_fraction))
    covariance = covariance + np.diag(floor)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 1e-12)
    eigenvectors = eigenvectors[:, order]
    centered = points - centroid
    squared_distances = np.einsum(
        "ij,jk,ik->i",
        centered,
        np.linalg.pinv(covariance),
        centered,
    )
    scale = float(radius_scale) * float(
        np.sqrt(np.quantile(np.maximum(squared_distances, 0.0), confidence))
    )
    width = 2.0 * scale * float(np.sqrt(eigenvalues[0]))
    height = 2.0 * scale * float(np.sqrt(eigenvalues[1]))
    angle = float(np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0])))
    return centroid, width, height, angle

def _distribution_limits(
    sessions: pd.DataFrame,
    embedding: np.ndarray,
    config: FinalSequenceFigureConfig,
    *,
    confidence: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Include every displayed ellipse in the common Panel A/B limits."""
    lower = np.asarray(embedding, dtype=float).min(axis=0)
    upper = np.asarray(embedding, dtype=float).max(axis=0)
    for paired_status, stage in PANEL_ORDER:
        panel_mask = (
            sessions["paired_status"].eq(paired_status)
            & sessions["stage"].eq(stage)
        ).to_numpy()
        for individual in INDIVIDUAL_ORDER:
            indices = np.flatnonzero(
                panel_mask & sessions["focal ID"].eq(individual).to_numpy()
            )
            ellipse = _distribution_ellipse_parameters(
                embedding[indices],
                embedding,
                confidence=confidence,
                radius_scale=config.distribution_ellipse_radius_scale,
            )
            if ellipse is None:
                continue
            centroid, width, height, angle = ellipse
            radians = np.deg2rad(angle)
            x_radius = 0.5 * np.hypot(width * np.cos(radians), height * np.sin(radians))
            y_radius = 0.5 * np.hypot(width * np.sin(radians), height * np.cos(radians))
            lower = np.minimum(lower, centroid - (x_radius, y_radius))
            upper = np.maximum(upper, centroid + (x_radius, y_radius))
    span = np.maximum(upper - lower, 1e-6)
    padding = span * 0.035
    return (
        (float(lower[0] - padding[0]), float(upper[0] + padding[0])),
        (float(lower[1] - padding[1]), float(upper[1] + padding[1])),
    )

def _draw_distribution_panel(
    axis,
    sessions: pd.DataFrame,
    embedding: np.ndarray,
    panel_key: tuple[str, str],
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
    config: FinalSequenceFigureConfig,
    *,
    confidence: float = 0.50,
) -> None:
    """Draw one subset with distributions and bonded-pair centroid lines."""
    paired_status, stage = panel_key
    base_mask = (
        sessions["paired_status"].eq(paired_status)
        & sessions["stage"].eq(stage)
    ).to_numpy()
    centroids: dict[str, np.ndarray] = {}
    for individual in INDIVIDUAL_ORDER:
        indices = np.flatnonzero(
            base_mask & sessions["focal ID"].eq(individual).to_numpy()
        )
        if not len(indices):
            continue
        color = FINAL_INDIVIDUAL_COLORS[individual]
        points = embedding[indices]
        ellipse = _distribution_ellipse_parameters(
            points,
            embedding,
            confidence=confidence,
            radius_scale=config.distribution_ellipse_radius_scale,
        )
        if ellipse is not None:
            centroid, width, height, angle = ellipse
            axis.add_patch(
                Ellipse(
                    xy=centroid,
                    width=width,
                    height=height,
                    angle=angle,
                    facecolor=color,
                    edgecolor="none",
                    alpha=0.11,
                    zorder=1,
                )
            )
            axis.add_patch(
                Ellipse(
                    xy=centroid,
                    width=width,
                    height=height,
                    angle=angle,
                    facecolor="none",
                    edgecolor=color,
                    linewidth=1.45,
                    alpha=0.82,
                    # Keep the boundary visible over raw points and pair
                    # links, while the enlarged centroid remains foremost.
                    zorder=3.5,
                )
            )
        axis.scatter(
            points[:, 0],
            points[:, 1],
            s=config.point_size * 1.20,
            marker=CONTEXT_MARKERS[paired_status],
            c=color,
            alpha=config.point_alpha,
            linewidths=0,
            zorder=2,
        )
        centroids[individual] = points.mean(axis=0)

    for pair_id, members in PAIR_MEMBERS.items():
        if all(member in centroids for member in members):
            axis.plot(
                [centroids[member][0] for member in members],
                [centroids[member][1] for member in members],
                color=FINAL_PAIR_LINE_COLORS[pair_id],
                linewidth=2.2,
                alpha=0.92,
                zorder=3,
            )
    for individual, centroid in centroids.items():
        axis.scatter(
            [centroid[0]],
            [centroid[1]],
            s=config.centroid_size,
            marker="o",
            c=FINAL_INDIVIDUAL_COLORS[individual],
            edgecolors="#111111",
            linewidths=1.4,
            zorder=4,
        )
    _style_embedding_axis(axis, x_limits, y_limits, config)
    axis.set_title(
        PANEL_TITLES[panel_key], fontsize=_font_size(config, 11.5), pad=5
    )

def render_sequence_distribution_metric_figure(
    sessions: pd.DataFrame,
    embedding: np.ndarray,
    metric: str,
    *,
    config: FinalSequenceFigureConfig | None = None,
    example_group_ids: Sequence[int] | None = None,
    confidence: float = 0.50,
) -> tuple[Figure, pd.DataFrame]:
    """Render the context-shaped, distribution-ellipse companion figure."""
    config = config or FinalSequenceFigureConfig()
    if metric not in METRIC_LABELS:
        raise ValueError(f"Unknown sequence metric: {metric!r}")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie strictly between zero and one.")
    if config.distribution_ellipse_radius_scale <= 0.0:
        raise ValueError("distribution_ellipse_radius_scale must be positive.")
    ordered = _ordered_sessions(sessions)
    embedding = np.asarray(embedding, dtype=np.float32)
    if embedding.shape != (len(ordered), 2) or not np.isfinite(embedding).all():
        raise ValueError("Embedding must be finite with shape (number of sessions, 2).")
    representatives = select_final_repertoire_examples(
        ordered,
        embedding,
        examples_per_side=config.examples_per_side_for_metric(metric),
        example_group_ids=example_group_ids,
        metric=metric,
    )
    heatmap_maximum = _metric_card_heatmap_maximum(ordered, metric)
    embedding_label = _embedding_method_label(metric)
    panel_x_limits, panel_y_limits = _distribution_limits(
        ordered,
        embedding,
        config,
        confidence=confidence,
    )
    # Panel A shows the observed repertoire cloud, whereas Panel B needs the
    # wider limits required to contain every displayed spread ellipse. Keeping separate
    # limits lets the embedding pattern fill Panel A without clipping Panel B.
    overview_x_limits, overview_y_limits = _embedding_limits(
        embedding,
        padding_fraction=(
            0.035 if metric in {"bigram", "transition_probability"} else 0.065
        ),
    )

    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": _font_size(config, 10.0),
            "axes.labelcolor": "#222222",
            "xtick.color": "#333333",
            "ytick.color": "#333333",
        }
    ):
        figure = plt.figure(figsize=config.figure_size_inches, facecolor="white")
        if metric == "local_alignment":
            overview_position = (0.125, 0.245, 0.260, 0.475)
        elif metric in {"bigram", "transition_probability"}:
            overview_position = (0.082, 0.255, 0.315, 0.475)
        else:
            overview_position = (0.090, 0.270, 0.305, 0.460)
        overview_axis = figure.add_axes(overview_position)
        _draw_distribution_overview(
            overview_axis,
            ordered,
            embedding,
            overview_x_limits,
            overview_y_limits,
            config,
            embedding_label,
        )
        overview_axis.xaxis.set_label_coords(0.5, 0.045)
        overview_axis.yaxis.set_label_coords(0.055, 0.5)
        overview_axis.tick_params(axis="both", labelbottom=False, labelleft=False)
        movement_axes = {}
        for panel_key in PANEL_ORDER:
            axis = figure.add_axes(PANEL_B_POSITIONS[panel_key])
            _draw_distribution_panel(
                axis,
                ordered,
                embedding,
                panel_key,
                panel_x_limits,
                panel_y_limits,
                config,
                confidence=confidence,
            )
            movement_axes[panel_key] = axis
        _format_panel_b_axes(movement_axes, embedding_label, config)
        _add_panel_b_shared_y_label(
            figure, movement_axes, embedding_label, config
        )

        metric_axis_slots = _heatmap_corner_slots(
            config.examples_per_side_for_metric(metric)
        )
        for (_, row), (position, side, card_layout) in zip(
            representatives.iterrows(),
            _packed_card_positions(representatives, config, metric),
        ):
            card_axis = figure.add_axes(position)
            _draw_metric_card(
                card_axis,
                row,
                card_layout,
                metric,
                config,
                heatmap_maximum=heatmap_maximum,
                border_side=side,
                scale_corner=(
                    metric_axis_slots.get(int(row["border_slot"]))
                    if metric in {"bigram", "transition_probability", "repeat_A"}
                    else None
                ),
            )
            group_id = int(row["group_id"])
            _connect_card_to_point(
                figure,
                card_axis,
                overview_axis,
                embedding[group_id],
                str(row["focal ID"]),
                side,
                config,
                marker=CONTEXT_MARKERS[str(row["paired_status"])],
            )

        _add_metric_card_colorbar(figure, metric, heatmap_maximum, config)
        if metric != "repeat_A":
            _add_call_type_legend(figure, text_scale=config.text_scale)

        individual_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                markerfacecolor=FINAL_INDIVIDUAL_COLORS[individual],
                markeredgecolor="#111111",
                markersize=_font_size(config, 8),
                label=INDIVIDUAL_DISPLAY_LABELS[individual],
            )
            for individual in INDIVIDUAL_LEGEND_ORDER
        ]
        context_handles = [
            Line2D(
                [0],
                [0],
                marker=CONTEXT_MARKERS[context],
                linestyle="None",
                markerfacecolor="#888888",
                markeredgecolor="none",
                markersize=_font_size(config, 8),
                label="Partner" if context == "partner" else "Non-partner",
            )
            for context in ("partner", "non-partner")
        ]
        ellipse_label = (
            f"{confidence:.0%} ellipse"
            if np.isclose(config.distribution_ellipse_radius_scale, 1.0)
            else (
                f"{confidence:.0%} ellipse "
                f"({config.distribution_ellipse_radius_scale:.0%} radius)"
            )
        )
        summary_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                markerfacecolor="#888888",
                markeredgecolor="#111111",
                markersize=_font_size(config, 9),
                label="Centroid",
            ),
            Line2D(
                [0],
                [0],
                color="#555555",
                linewidth=2.2,
                label="Pair link",
            ),
            Ellipse(
                (0, 0),
                width=1.2,
                height=0.7,
                facecolor="#999999",
                edgecolor="#555555",
                alpha=0.30,
                label=ellipse_label,
            ),
        ]
        figure.legend(
            handles=individual_handles,
            loc="lower right",
            bbox_to_anchor=(0.982, 0.012),
            ncol=2,
            borderaxespad=0,
            frameon=False,
            fontsize=_font_size(config, 9.2),
            handlelength=0.8,
            columnspacing=1.30,
            handletextpad=0.25,
            labelspacing=0.32,
        )
        figure.legend(
            handles=context_handles + summary_handles,
            loc="lower right",
            bbox_to_anchor=(0.795, 0.012),
            ncol=2,
            borderaxespad=0,
            frameon=False,
            fontsize=_font_size(config, 8.6),
            handlelength=1.0,
            labelspacing=0.24,
            columnspacing=0.40,
            handletextpad=0.25,
        )
        figure.text(
            0.012,
            0.985,
            "A",
            va="top",
            fontsize=_font_size(config, 15),
            fontweight="bold",
        )
        figure.text(
            0.245,
            0.985,
            f"{METRIC_LABELS[metric]} space",
            ha="center",
            va="top",
            fontsize=_font_size(config, 13),
            fontweight="bold",
        )
        figure.text(
            0.460,
            0.985,
            "B",
            va="top",
            fontsize=_font_size(config, 15),
            fontweight="bold",
        )
        figure.text(
            0.7175,
            0.985,
            "Spread by stage and context",
            ha="center",
            va="top",
            fontsize=_font_size(config, 13),
            fontweight="bold",
        )
    return figure, representatives

PUBLIC_METRICS = {
    "transition_probability": "transition_probability",
    "bigram": "bigram",
    "phee_repeat": "repeat_A",
    "local_alignment": "local_alignment",
}


def _config_from_saved_json(saved: Mapping[str, object]) -> FinalSequenceFigureConfig:
    values = dict(saved["figure_config"])
    values["umap"] = UmapConfig(**values["umap"])
    values["pacmap"] = PacmapConfig(**values["pacmap"])
    return FinalSequenceFigureConfig(**values)


def _validate_saved_inputs(
    sessions: pd.DataFrame,
    embedding_table: pd.DataFrame,
    examples: pd.DataFrame,
    subset_counts: pd.DataFrame,
    saved_config: Mapping[str, object],
    internal_metric: str,
) -> np.ndarray:
    """Check that every staged plotting input belongs to the same 107 rows."""
    ordered = _ordered_sessions(sessions)
    embedding_table = embedding_table.sort_values("group_id").reset_index(drop=True)
    if not np.array_equal(
        embedding_table["group_id"].to_numpy(dtype=int),
        ordered["group_id"].to_numpy(dtype=int),
    ):
        raise ValueError(f"{internal_metric}: embedding and session orders differ")
    for column in (
        "focal ID", "conspecific_ID", "stage", "paired_status",
        "session_number", "pair_id", "sex",
    ):
        if not np.array_equal(
            embedding_table[column].astype(str).to_numpy(),
            ordered[column].astype(str).to_numpy(),
        ):
            raise ValueError(f"{internal_metric}: mismatched embedding metadata: {column}")
    if saved_config["metric"] != internal_metric:
        raise ValueError(f"{internal_metric}: saved config names a different metric")
    expected_examples = [int(value) for value in saved_config["example_group_ids"]]
    observed_examples = (
        examples.sort_values("border_order")["group_id"].drop_duplicates().astype(int).tolist()
    )
    if observed_examples != expected_examples:
        raise ValueError(f"{internal_metric}: example selections do not match the config")
    observed_counts = (
        ordered.groupby(["focal ID", "stage", "paired_status"], sort=True)
        .size().rename("n_session_repertoires").reset_index()
    )
    observed_counts["ellipse_drawn"] = observed_counts["n_session_repertoires"].ge(2)
    expected_counts = subset_counts.sort_values(
        ["focal ID", "stage", "paired_status"]
    ).reset_index(drop=True)
    observed_counts = observed_counts.sort_values(
        ["focal ID", "stage", "paired_status"]
    ).reset_index(drop=True)
    columns = ["focal ID", "stage", "paired_status", "n_session_repertoires", "ellipse_drawn"]
    if not observed_counts[columns].equals(expected_counts[columns]):
        raise ValueError(f"{internal_metric}: saved subset counts do not match the sessions")
    coordinates = embedding_table[["embedding_1", "embedding_2"]].to_numpy(dtype=np.float32)
    if coordinates.shape != (107, 2) or not np.isfinite(coordinates).all():
        raise ValueError(f"{internal_metric}: expected 107 finite two-dimensional coordinates")
    return coordinates


def reproduce_sequence_distribution_figures(project_dir: str | Path) -> dict[str, Path]:
    """Recreate the four final 50% distribution figures from bundled inputs."""
    project_dir = Path(project_dir).resolve()
    input_dir = project_dir / "data" / "sequence" / "figure_inputs"
    output_dir = project_dir / "results" / "sequence_spaces"
    output_dir.mkdir(parents=True, exist_ok=True)
    sessions = pd.read_csv(project_dir / "data" / "sequence" / "session_order_107.csv")
    outputs: dict[str, Path] = {}

    for public_metric, internal_metric in PUBLIC_METRICS.items():
        prefix = input_dir / public_metric
        saved_config = json.loads(prefix.with_name(f"{public_metric}_config.json").read_text())
        embedding_table = pd.read_csv(prefix.with_name(f"{public_metric}_embedding.csv"))
        examples = pd.read_csv(prefix.with_name(f"{public_metric}_examples.csv"))
        subset_counts = pd.read_csv(prefix.with_name(f"{public_metric}_subset_counts.csv"))
        coordinates = _validate_saved_inputs(
            sessions, embedding_table, examples, subset_counts,
            saved_config, internal_metric,
        )
        config = _config_from_saved_json(saved_config)
        figure, _ = render_sequence_distribution_metric_figure(
            sessions,
            coordinates,
            internal_metric,
            config=config,
            example_group_ids=saved_config["example_group_ids"],
            confidence=float(saved_config["confidence"]),
        )
        stem = output_dir / f"sequence_space_{public_metric}_distribution50_reproduced"
        path = stem.with_suffix(".png")
        figure.savefig(path, dpi=config.dpi, facecolor="white")
        plt.close(figure)
        outputs[public_metric] = path
    return outputs
