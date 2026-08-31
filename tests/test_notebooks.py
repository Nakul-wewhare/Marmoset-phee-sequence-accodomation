from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_five_numbered_notebooks_are_clean_and_syntax_valid() -> None:
    notebooks = sorted((ROOT / "notebooks").glob("[0-9][0-9]_*.ipynb"))
    assert [path.name[:2] for path in notebooks] == ["01", "02", "03", "04", "05"]

    for path in notebooks:
        payload = json.loads(path.read_text(encoding="utf-8"))
        combined_source = "\n".join(
            "".join(cell.get("source", [])) for cell in payload["cells"]
        )
        assert "EXPENSIVE_STEPS" in combined_source
        assert "recompute_dtw\": True" not in combined_source
        assert "retrain_vae\": True" not in combined_source
        assert "refit_models\": True" not in combined_source
        assert "C:\\Users\\" not in combined_source
        assert "/Users/" not in combined_source
        assert "trusted read" not in combined_source.casefold()

        code_cells = [
            cell for cell in payload["cells"] if cell.get("cell_type") == "code"
        ]
        assert code_cells
        for index, cell in enumerate(code_cells):
            assert cell.get("outputs", []) == []
            assert cell.get("execution_count") is None
            compile(
                "".join(cell.get("source", [])),
                f"{path.name}:code-cell-{index}",
                "exec",
            )


def test_notebooks_link_to_public_specification() -> None:
    for path in sorted((ROOT / "notebooks").glob("[0-9][0-9]_*.ipynb")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        markdown = "\n".join(
            "".join(cell.get("source", []))
            for cell in payload["cells"]
            if cell.get("cell_type") == "markdown"
        )
        assert "manuscript" in markdown.casefold()
        assert "docs/" in markdown or "../docs/" in markdown
