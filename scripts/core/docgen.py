"""Auto-documentation generator — maintains per-comparison Markdown reports.

Each analysis step calls `update_section()` to write/replace its section
in the report.  Sections are delimited by HTML comment markers so that
re-running a step overwrites only its own content.

Typical usage inside a step's run() function:

    from ..docgen import update_section
    update_section(cfg, "3a", "DE Screening", content)

where `content` is a Markdown string (without the heading — it's added
automatically).
"""

import os
import re
from datetime import datetime

from .config_runtime import get_batch_id, get_conditions, get_pair_name, get_path, get_timepoints
from .logger import setup_logger
from .plotting import plot_filename

log = setup_logger("docgen")


def _load_step_order():
    """Derive step order from the registry YAML (public steps only)."""
    try:
        from .runner import _ensure_registry, STEPS
        _ensure_registry()
        return [sid for sid, _, _ in STEPS]
    except Exception:
        # Fallback if runner fails to load
        return [
            "1a",
            "2a", "2b", "2c", "2d",
            "3a", "3b", "3c", "3d",
            "4a", "4b", "4c", "4d",
            "5a", "5b",
        ]


STEP_ORDER = _load_step_order()


def _doc_path(cfg, pair_name):
    """Return absolute path to the comparison Markdown file."""
    pipeline_root = cfg["_pipeline_root"]
    batch_id = get_batch_id(cfg)
    doc_dir = cfg.get("documents_dir", "results/{batch_id}/documents")
    doc_dir = doc_dir.format(batch_id=batch_id)
    if not os.path.isabs(doc_dir):
        doc_dir = os.path.join(pipeline_root, doc_dir)
    os.makedirs(doc_dir, exist_ok=True)
    return os.path.join(doc_dir, f"{pair_name}.md")


def _doc_dir(cfg):
    """Return absolute document directory for this config."""
    return os.path.dirname(_doc_path(cfg, get_pair_name(cfg)))


def _results_relpath(cfg):
    """Return relative path from documents dir to comparison results dir."""
    doc_dir = _doc_dir(cfg)
    results_dir = get_path(cfg, "output_dir")
    return os.path.relpath(results_dir, doc_dir)


def _shared_relpath(cfg):
    """Return relative path from documents dir to shared results dir."""
    doc_dir = _doc_dir(cfg)
    shared_dir = get_path(cfg, "shared_output_dir")
    return os.path.relpath(shared_dir, doc_dir)


def _start_marker(step_id):
    return f"<!-- SECTION:{step_id}:START -->"


def _end_marker(step_id):
    return f"<!-- SECTION:{step_id}:END -->"


def _header(cfg, pair_name):
    """Generate the document header."""
    project = cfg["project"]
    conds = pair_name.split("-")
    timepoints = get_timepoints(cfg)
    return (
        f"# {project['name']} Analysis Report: {' vs '.join(conds)}\n\n"
        f"**Organism:** {project['organism']} (strain {project['strain']})  \n"
        f"**Batch:** {get_batch_id(cfg)}  \n"
        f"**Conditions:** {conds[0]} vs {conds[1]}  \n"
        f"**Timepoints:** {', '.join(str(t) + 'h' for t in timepoints)}  \n"
        f"**Samples per condition:** {len(timepoints)}  \n\n"
        f"---\n\n"
    )


def _insert_section(full_text, step_id, section_block):
    """Insert or replace a section, keeping step order."""
    start = _start_marker(step_id)
    end = _end_marker(step_id)

    # Replace existing section
    pattern = re.compile(
        re.escape(start) + r".*?" + re.escape(end),
        re.DOTALL,
    )
    if pattern.search(full_text):
        return pattern.sub(section_block, full_text)

    # Insert at correct position
    step_idx = STEP_ORDER.index(step_id) if step_id in STEP_ORDER else len(STEP_ORDER)
    for later_id in STEP_ORDER[step_idx + 1:]:
        later_start = _start_marker(later_id)
        pos = full_text.find(later_start)
        if pos != -1:
            return full_text[:pos] + section_block + "\n" + full_text[pos:]

    # Append at end
    return full_text.rstrip("\n") + "\n\n" + section_block + "\n"


def get_comparison_pairs(cfg):
    """Return list of comparison pair names from config."""
    comparisons = cfg.get("comparisons")
    if comparisons:
        return [f"{c[0]}-{c[1]}" for c in comparisons]
    # Default: first two conditions
    conds = get_conditions(cfg)
    if len(conds) >= 2:
        return [f"{conds[0]}-{conds[1]}"]
    return []


def update_section(cfg, step_id, step_name, content):
    """Write or replace a step's section in all comparison documents.

    Parameters
    ----------
    cfg : dict
        Pipeline config.
    step_id : str
        Step identifier, e.g. "3a".
    step_name : str
        Human-readable step name, e.g. "DE Screening".
    content : str
        Markdown body for this section (without heading).
    """
    pairs = get_comparison_pairs(cfg)
    if not pairs:
        return

    start = _start_marker(step_id)
    end = _end_marker(step_id)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    section_block = (
        f"{start}\n"
        f"## {step_id}. {step_name}\n\n"
        f"*Updated: {timestamp}*\n\n"
        f"{content.rstrip()}\n\n"
        f"{end}"
    )

    for pair_name in pairs:
        doc_file = _doc_path(cfg, pair_name)

        if os.path.exists(doc_file):
            with open(doc_file, "r") as f:
                full_text = f.read()
        else:
            full_text = _header(cfg, pair_name)

        full_text = _insert_section(full_text, step_id, section_block)

        with open(doc_file, "w") as f:
            f.write(full_text)

    log.info(f"  Doc updated: {step_id} ({step_name})")


def img(cfg, subdir, filename, caption=""):
    """Image reference to a comparison-specific results file (steps 1a-5b)."""
    rel = _results_relpath(cfg)
    output_name = plot_filename(cfg, filename)
    path = f"{rel}/{subdir}/{output_name}"
    alt = caption or output_name
    return f"![{alt}]({path})"


def shared_img(cfg, subdir, filename, caption=""):
    """Image reference to a shared results file."""
    rel = _shared_relpath(cfg)
    output_name = plot_filename(cfg, filename)
    path = f"{rel}/{subdir}/{output_name}"
    alt = caption or output_name
    return f"![{alt}]({path})"
