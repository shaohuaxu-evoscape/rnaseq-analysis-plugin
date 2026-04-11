"""Persisted runner state and resume/checkpoint helpers."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone

from .config_runtime import get_batch_id, get_pair_name, get_path

STATE_FILENAME = "pipeline_run_manifest.json"
STATE_FORMAT_VERSION = 1


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def _sanitize(value):
    if isinstance(value, dict):
        return {
            str(key): _sanitize(item)
            for key, item in value.items()
            if str(key) != "_data" and not str(key).startswith("__")
        }
    if isinstance(value, list):
        return [_sanitize(item) for item in value]
    if isinstance(value, tuple):
        return [_sanitize(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def config_fingerprint(cfg):
    """Return a stable fingerprint for the resolved config."""
    payload = _sanitize(cfg)
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def get_state_path(cfg):
    """Return the persisted run-manifest path for this output directory."""
    return os.path.join(get_path(cfg, "output_dir"), STATE_FILENAME)


def load_state(cfg):
    """Load persisted run manifest, or return None if unavailable."""
    state_path = get_state_path(cfg)
    if not os.path.isfile(state_path):
        return None
    try:
        with open(state_path, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return None


def create_state(cfg, config_path, requested_steps, step_ids, auto_added_steps, resume_enabled):
    """Create a new persisted run manifest."""
    return {
        "format_version": STATE_FORMAT_VERSION,
        "config_fingerprint": config_fingerprint(cfg),
        "config_path": os.path.abspath(config_path) if config_path else "",
        "batch_id": get_batch_id(cfg),
        "pair_name": get_pair_name(cfg),
        "output_dir": get_path(cfg, "output_dir"),
        "requested_steps": list(requested_steps),
        "resolved_steps": list(step_ids),
        "auto_added_steps": list(auto_added_steps),
        "resume_enabled": bool(resume_enabled),
        "status": "running",
        "started_at": _utc_now(),
        "updated_at": _utc_now(),
        "completed_at": "",
        "steps": {},
    }


def save_state(cfg, state):
    """Persist run manifest to disk."""
    state_path = get_state_path(cfg)
    os.makedirs(os.path.dirname(state_path), exist_ok=True)
    state["updated_at"] = _utc_now()
    with open(state_path, "w") as f:
        json.dump(state, f, indent=2, sort_keys=True)


def state_matches_config(cfg, state):
    """Return True when a persisted state matches the current resolved config."""
    if not state:
        return False
    if state.get("format_version") != STATE_FORMAT_VERSION:
        return False
    return state.get("config_fingerprint") == config_fingerprint(cfg)
