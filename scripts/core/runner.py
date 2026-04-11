"""Pipeline runner — orchestrates all steps in sequence.

Step definitions are loaded from configs/step_registry.yaml.
To register a new module, edit that YAML file and create the step .py files.
"""

import importlib
import os
import re
import time
from collections import OrderedDict
import yaml
from .config_runtime import load_config, get_batch_id, get_conditions, get_pair_name, get_timepoints
from .logger import setup_logger
from .pipeline_state import create_state, load_state, save_state, state_matches_config, _utc_now


# ── Registry loader ──────────────────────────────────────────────────────────

def _load_registry(registry_path=None):
    """Load step registry from YAML and build runtime data structures."""
    if registry_path is None:
        from .config_runtime import PLUGIN_ROOT
        registry_path = os.path.join(PLUGIN_ROOT, "configs", "step_registry.yaml")
    try:
        with open(registry_path) as f:
            reg = yaml.safe_load(f)
    except FileNotFoundError:
        raise RuntimeError(
            f"Step registry not found: {registry_path}\n"
            "Create configs/step_registry.yaml or pass a custom path."
        )

    preprocessing_steps = []
    public_steps = []
    public_modules = {}
    all_modules = {}
    module_display_names = {}
    hidden_step_ids = set()
    step_depends = {}
    seen_ids = set()

    # Sort modules by module_number to guarantee ordering
    sorted_modules = sorted(
        reg["modules"].items(),
        key=lambda kv: kv[1].get("module_number", 999),
    )

    for mod_key, mod_def in sorted_modules:
        is_public = mod_def.get("public", True)
        mod_num = str(mod_def["module_number"])
        step_ids = []

        for step in mod_def.get("steps", []):
            sid = step["id"]
            if sid in seen_ids:
                raise ValueError(f"Duplicate step ID '{sid}' in registry")
            seen_ids.add(sid)
            entry = (sid, step["module_path"], step["display_name"])
            step_ids.append(sid)
            step_depends[sid] = step.get("depends_on", [])
            if is_public:
                public_steps.append(entry)
            else:
                preprocessing_steps.append(entry)
                hidden_step_ids.add(sid)

        all_modules[mod_key] = step_ids
        if is_public:
            public_modules[mod_key] = step_ids
            module_display_names[mod_num] = mod_def["display_name"]

    return {
        "PREPROCESSING_STEPS": preprocessing_steps,
        "STEPS": public_steps,
        "ALL_STEPS": preprocessing_steps + public_steps,
        "PUBLIC_MODULES": public_modules,
        "ALL_MODULES": all_modules,
        "MODULE_DISPLAY_NAMES": module_display_names,
        "HIDDEN_STEP_IDS": hidden_step_ids,
        "STEP_DEPENDS": step_depends,
    }


# Lazy-loaded registry — populated on first access via _ensure_registry().
_registry = None
PREPROCESSING_STEPS = None
STEPS = None
ALL_STEPS = None
PUBLIC_MODULES = None
ALL_MODULES = None
MODULE_DISPLAY_NAMES = None
HIDDEN_STEP_IDS = None
STEP_DEPENDS = None


def _ensure_registry():
    """Load step registry on first use (lazy init)."""
    global _registry, PREPROCESSING_STEPS, STEPS, ALL_STEPS, PUBLIC_MODULES
    global ALL_MODULES, MODULE_DISPLAY_NAMES, HIDDEN_STEP_IDS, STEP_DEPENDS
    if _registry is not None:
        return
    _registry = _load_registry()
    PREPROCESSING_STEPS = _registry["PREPROCESSING_STEPS"]
    STEPS = _registry["STEPS"]
    ALL_STEPS = _registry["ALL_STEPS"]
    PUBLIC_MODULES = _registry["PUBLIC_MODULES"]
    ALL_MODULES = _registry["ALL_MODULES"]
    MODULE_DISPLAY_NAMES = _registry["MODULE_DISPLAY_NAMES"]
    HIDDEN_STEP_IDS = _registry["HIDDEN_STEP_IDS"]
    STEP_DEPENDS = _registry["STEP_DEPENDS"]

log = setup_logger("runner")


# ── Step parsing ─────────────────────────────────────────────────────────────

def _import_step(module_path):
    """Dynamically import a step module from scripts.steps."""
    # __name__ is "scripts.core.runner" → go up to "scripts", then into "steps"
    scripts_pkg = __name__.rsplit(".", 2)[0]
    full_path = f"{scripts_pkg}.steps.{module_path}"
    return importlib.import_module(full_path)


def _parse_steps(step_spec, allow_hidden=False):
    """Parse step specification string into list of step IDs.

    Supported formats:
      - "all"           → all steps
      - "1a,2a,3d"      → specific steps
      - "1a-2d"         → range within visible modules
      - "normalization" → all steps in a module
      - "1-3"           → visible modules 1 through 3
    """
    _ensure_registry()
    step_defs = ALL_STEPS if allow_hidden else STEPS
    module_map = ALL_MODULES if allow_hidden else PUBLIC_MODULES

    if step_spec is None or step_spec.lower() == "all":
        return [s[0] for s in step_defs]

    all_ids = [s[0] for s in step_defs]
    id_positions = {sid: idx for idx, sid in enumerate(all_ids)}
    visible_modules = sorted({
        int(m.group(1)) for sid in all_ids
        if (m := re.match(r"(\d+)", sid))
    })
    requested = []

    for part in step_spec.split(","):
        part = part.strip().lower()
        if not part:
            continue

        # Module name: "normalization", "preprocessing"
        if part in module_map:
            requested.extend(module_map[part])
            continue

        # Module number range: "1-3"
        if re.fullmatch(r"\d+-\d+", part):
            start_module, end_module = int(part.split("-")[0]), int(part.split("-")[1])
            if start_module > end_module:
                raise ValueError(f"Invalid module range '{part}': start is after end")
            module_steps = [
                sid for sid in all_ids
                if (m := re.match(r"(\d+)", sid)) and start_module <= int(m.group(1)) <= end_module
            ]
            if not module_steps:
                raise ValueError(
                    f"Invalid module range '{part}'. Available modules: "
                    f"{', '.join(str(m) for m in visible_modules)}"
                )
            requested.extend(module_steps)
            continue

        # Step range: "1a-2d"
        if re.fullmatch(r"\d[a-z]-\d[a-z]", part):
            start, end = part.split("-")
            if start not in id_positions or end not in id_positions:
                if not allow_hidden and (start.startswith("0") or end.startswith("0")):
                    raise ValueError(
                        "Preprocessing is now automatic and hidden. "
                        "Run public analysis steps only, for example '--steps 1a-2d' or '--steps all'."
                    )
                raise ValueError(
                    f"Invalid step range '{part}'. Available steps: {', '.join(all_ids)}"
                )
            if id_positions[start] > id_positions[end]:
                raise ValueError(f"Invalid step range '{part}': start is after end")
            requested.extend(all_ids[id_positions[start]:id_positions[end] + 1])
            continue

        # Single step ID
        if part in id_positions:
            requested.append(part)
            continue
        if not allow_hidden and (part in HIDDEN_STEP_IDS or part == "preprocessing"):
            raise ValueError(
                "Preprocessing is now automatic and hidden. "
                "Run public analysis steps only, for example '--steps 1a-2d' or '--steps all'."
            )
        raise ValueError(
            f"Unrecognized step selector '{part}'. Available steps: {', '.join(all_ids)}"
        )

    # Deduplicate while preserving order
    seen = set()
    result = []
    for sid in requested:
        if sid not in seen:
            seen.add(sid)
            result.append(sid)
    return result


def _resolve_step_plan(step_ids, allow_hidden=False):
    """Expand dependencies and return a topologically ordered execution plan."""
    _ensure_registry()
    available_ids = {
        sid for sid, _, _ in (ALL_STEPS if allow_hidden else STEPS)
    }
    registry_order = [sid for sid, _, _ in ALL_STEPS if sid in available_ids]
    closure = OrderedDict()
    visiting = set()

    def visit(sid):
        if sid not in available_ids:
            if not allow_hidden and sid in HIDDEN_STEP_IDS:
                raise ValueError(
                    "Preprocessing is now automatic and hidden. "
                    "Run public analysis steps only, for example '--steps 1a-2d' or '--steps all'."
                )
            raise ValueError(f"Step '{sid}' is not available in this execution mode")
        if sid in closure:
            return
        if sid in visiting:
            raise ValueError(f"Cyclic dependency detected at step '{sid}'")
        visiting.add(sid)
        for dep in STEP_DEPENDS.get(sid, []):
            visit(dep)
        visiting.remove(sid)
        closure[sid] = None

    for sid in step_ids:
        visit(sid)

    ordered = [sid for sid in registry_order if sid in closure]
    requested_set = set(step_ids)
    auto_added = [sid for sid in ordered if sid not in requested_set]
    return ordered, auto_added


# ── Pipeline execution ───────────────────────────────────────────────────────

def run_pipeline(config_path=None, steps="all", dry_run=False, overrides=None, allow_hidden=False, resume=False):
    """Run the pipeline."""
    _ensure_registry()
    cfg = load_config(config_path, overrides=overrides)
    requested_steps = _parse_steps(steps, allow_hidden=allow_hidden)
    step_ids, auto_added = _resolve_step_plan(requested_steps, allow_hidden=allow_hidden)

    step_lookup = {s[0]: s for s in ALL_STEPS}
    os.makedirs(cfg["paths"]["output_dir"], exist_ok=True)
    prior_state = load_state(cfg)
    resume_state = prior_state if resume and state_matches_config(cfg, prior_state) else None
    if resume and prior_state and resume_state is None:
        log.warning("Existing pipeline_run_manifest.json does not match the current config; starting a fresh run")
    state = create_state(
        cfg=cfg,
        config_path=config_path,
        requested_steps=requested_steps,
        step_ids=step_ids,
        auto_added_steps=auto_added,
        resume_enabled=resume,
    )
    if resume_state:
        state["steps"] = dict(resume_state.get("steps", {}))
        state["started_at"] = resume_state.get("started_at", state["started_at"])
    save_state(cfg, state)

    conds = get_conditions(cfg)
    tps = get_timepoints(cfg)

    log.info("=" * 60)
    log.info(f"RNA-seq Pipeline — Batch {get_batch_id(cfg)}")
    log.info(f"Conditions: {' vs '.join(conds)}  |  Timepoints: {', '.join(str(t)+'h' for t in tps)}")
    log.info(f"Output: {get_pair_name(cfg)}/  |  Steps: {', '.join(step_ids)}")
    if auto_added:
        log.info(f"Auto-added dependencies: {', '.join(auto_added)}")
    if resume_state:
        completed_steps = [
            sid for sid, data in resume_state.get("steps", {}).items()
            if data.get("status") == "ok" and sid in step_ids
        ]
        if completed_steps:
            log.info(f"Resume checkpoint found for completed steps: {', '.join(completed_steps)}")
    log.info("=" * 60)

    if dry_run:
        resume_candidates = []
        for sid in step_ids:
            _, mod_path, name = step_lookup[sid]
            if resume_state and resume_state.get("steps", {}).get(sid, {}).get("status") == "ok":
                resume_candidates.append(sid)
            log.info(f"  [DRY RUN] {sid}: {name} ({mod_path})")
        return {
            "status": "dry_run",
            "requested_steps": requested_steps,
            "steps": step_ids,
            "auto_added_steps": auto_added,
            "resume_candidates": resume_candidates,
        }

    results = {}
    t_start = time.time()

    for sid in step_ids:
        _, mod_path, name = step_lookup[sid]
        log.info(f"\n{'─' * 40}")
        t0 = time.time()

        prior_step = state.get("steps", {}).get(sid, {})
        if resume_state and prior_step.get("status") == "ok":
            results[sid] = {
                "status": "skipped",
                "reason": "resume_checkpoint",
                "time": 0.0,
            }
            state["steps"][sid] = {
                **prior_step,
                "status": "ok",
                "module_path": mod_path,
                "display_name": name,
                "resumed_at": _utc_now(),
            }
            save_state(cfg, state)
            log.info(f"  ↷ {name} skipped (resume checkpoint)")
            continue

        try:
            mod = _import_step(mod_path)
            mod.run(cfg)
            elapsed = time.time() - t0
            results[sid] = {"status": "ok", "time": elapsed}
            state["steps"][sid] = {
                "status": "ok",
                "module_path": mod_path,
                "display_name": name,
                "time": elapsed,
                "completed_at": _utc_now(),
            }
            save_state(cfg, state)
            log.info(f"  ✓ {name} completed ({elapsed:.1f}s)")
        except Exception as e:
            elapsed = time.time() - t0
            results[sid] = {"status": "error", "error": str(e), "time": elapsed}
            state["steps"][sid] = {
                "status": "error",
                "module_path": mod_path,
                "display_name": name,
                "time": elapsed,
                "error": str(e),
                "completed_at": _utc_now(),
            }
            save_state(cfg, state)
            log.error(f"  ✗ {name} FAILED: {e}", exc_info=True)
            if not cfg.get("continue_on_error", False):
                log.error("  Pipeline halted (set continue_on_error: true to keep going)")
                break

    total_time = time.time() - t_start
    n_ok = sum(1 for r in results.values() if r["status"] == "ok")
    n_err = sum(1 for r in results.values() if r["status"] == "error")
    n_skipped = sum(1 for r in results.values() if r["status"] == "skipped")

    state["status"] = "complete"
    state["completed_at"] = _utc_now()
    state["results"] = results
    state["total_time"] = total_time
    save_state(cfg, state)

    log.info(f"\n{'=' * 60}")
    log.info(f"Pipeline complete: {n_ok} succeeded, {n_err} failed, {n_skipped} skipped ({total_time:.1f}s)")
    log.info("=" * 60)

    return {
        "status": "complete",
        "requested_steps": requested_steps,
        "steps": step_ids,
        "auto_added_steps": auto_added,
        "results": results,
        "total_time": total_time,
    }
