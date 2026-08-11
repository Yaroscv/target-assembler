from __future__ import annotations

import copy
import csv
import json
import re
import time
from pathlib import Path
from typing import Any

import yaml

from .config import (
    DEFAULTS,
    apply_overrides,
    config_from_mapping,
    read_config_mapping,
    set_dotted_value,
)
from .errors import ConfigError
from .pipeline import Pipeline


ROW_ALIASES = {
    "reads1": "input.reads1",
    "reads2": "input.reads2",
    "layout": "input.layout",
    "hmm": "input.hmm",
    "clean_reads1": "input.clean_reads1",
    "clean_reads2": "input.clean_reads2",
    "sample_reads1": "input.sample_reads1",
    "sample_reads2": "input.sample_reads2",
    "variants_fasta": "input.variants_fasta",
    "validation_reads1": "input.validation_reads1",
    "validation_reads2": "input.validation_reads2",
    "mode": "project.mode",
    "outdir": "project.outdir",
    "threads": "runtime.threads",
    "seed": "runtime.seed",
    "validation_profile": "validation.profile",
    "assembly_technology": "assembly.technology",
}

PATH_KEYS = {
    "input.reads1",
    "input.reads2",
    "input.hmm",
    "input.clean_reads1",
    "input.clean_reads2",
    "input.sample_reads1",
    "input.sample_reads2",
    "input.variants_fasta",
    "input.validation_reads1",
    "input.validation_reads2",
    "project.outdir",
    "project.temp_dir",
}


def safe_sample_name(value: str) -> str:
    cleaned = re.sub(r"[^\w.-]+", "_", value, flags=re.UNICODE).strip("._")
    if not cleaned or cleaned in {".", ".."}:
        raise ConfigError(f"Sample name cannot be converted to a safe directory name: {value!r}")
    return cleaned


def _read_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with open(path, encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            try:
                raw_header = next(reader)
            except StopIteration as exc:
                raise ConfigError(f"Batch table is empty: {path}") from exc
            header = [item.strip() for item in raw_header]
            if any(not item for item in header):
                raise ConfigError(f"Batch table contains an empty column name: {path}")
            if len(set(header)) != len(header):
                raise ConfigError(f"Batch table contains duplicate column names: {path}")
            if "sample" not in header:
                raise ConfigError("Batch TSV must contain a 'sample' column")
            rows: list[dict[str, str]] = []
            for row_number, values in enumerate(reader, start=2):
                if not values or all(not item.strip() for item in values):
                    continue
                if len(values) != len(header):
                    raise ConfigError(
                        f"Batch TSV row {row_number} has {len(values)} fields; expected {len(header)}"
                    )
                rows.append({key: value.strip() for key, value in zip(header, values)})
    except OSError as exc:
        raise ConfigError(f"Cannot read batch table {path}: {exc}") from exc
    if not rows:
        raise ConfigError(f"Batch table contains no sample rows: {path}")
    return rows


def _parse_cell(text: str) -> Any:
    try:
        return yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise ConfigError(f"Cannot parse TSV value {text!r} as a YAML scalar") from exc


def _row_key(column: str) -> str:
    dotted = ROW_ALIASES.get(column, column)
    if dotted.startswith("batch.") or dotted == "batch":
        raise ConfigError("Batch settings cannot be overridden inside a sample row")
    top_level = dotted.split(".", 1)[0]
    if "." not in dotted or top_level not in DEFAULTS:
        aliases = ", ".join(sorted(ROW_ALIASES))
        raise ConfigError(
            f"Unknown batch column {column!r}; use a dotted YAML key or one of: sample, {aliases}"
        )
    return dotted


def _absolute_row_path(value: Any, table_dir: Path, column: str) -> Any:
    if value is None:
        return None
    if not isinstance(value, str):
        raise ConfigError(f"Path column {column!r} must contain a path or null")
    path = Path(value).expanduser()
    return str(path if path.is_absolute() else (table_dir / path).resolve())


def load_batch_configs(
    path: str | Path,
    overrides: list[str] | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Create validated per-sample configurations from one YAML file and TSV table."""
    config_path, raw = read_config_mapping(path)
    apply_overrides(raw, overrides)
    batch_raw = copy.deepcopy(DEFAULTS["batch"])
    batch_raw.update(copy.deepcopy(raw.get("batch", {})))
    table_value = batch_raw.get("samples_tsv")
    if not table_value:
        raise ConfigError("batch.samples_tsv is required for the batch command")
    table_path = Path(str(table_value)).expanduser()
    if not table_path.is_absolute():
        table_path = (config_path.parent / table_path).resolve()
    continue_on_error = batch_raw.get("continue_on_error", True)
    if not isinstance(continue_on_error, bool):
        raise ConfigError("batch.continue_on_error must be true or false")

    base_mapping = copy.deepcopy(raw)
    base_mapping.pop("batch", None)
    project_mapping = base_mapping.setdefault("project", {})
    root_value = project_mapping.get("outdir", DEFAULTS["project"]["outdir"])
    root_outdir = Path(str(root_value)).expanduser()
    if not root_outdir.is_absolute():
        root_outdir = (config_path.parent / root_outdir).resolve()

    rows = _read_tsv(table_path)
    configs: list[dict[str, Any]] = []
    seen_samples: set[str] = set()
    seen_outdirs: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        sample = row["sample"].strip()
        if not sample:
            raise ConfigError(f"Batch TSV row {row_number} has an empty sample name")
        if sample in seen_samples:
            raise ConfigError(f"Duplicate sample name in batch TSV: {sample!r}")
        seen_samples.add(sample)

        sample_mapping = copy.deepcopy(base_mapping)
        set_dotted_value(sample_mapping, "project.name", sample)
        set_dotted_value(sample_mapping, "project.outdir", str(root_outdir / safe_sample_name(sample)))
        for column, text in row.items():
            if column == "sample" or text == "":
                continue
            dotted = _row_key(column)
            value = _parse_cell(text)
            if dotted in PATH_KEYS:
                value = _absolute_row_path(value, table_path.parent, column)
            set_dotted_value(sample_mapping, dotted, value)

        try:
            config = config_from_mapping(sample_mapping, config_path)
        except ConfigError as exc:
            raise ConfigError(f"Batch TSV row {row_number} ({sample}): {exc}") from exc
        outdir = str(Path(config["project"]["outdir"]).resolve())
        if outdir in seen_outdirs:
            raise ConfigError(f"Samples resolve to the same output directory: {outdir}")
        seen_outdirs.add(outdir)
        config["_batch_sample"] = sample
        config["_batch_row"] = row_number
        configs.append(config)

    settings = {
        "config": str(config_path),
        "samples_tsv": str(table_path),
        "outdir": str(root_outdir),
        "continue_on_error": continue_on_error,
    }
    return settings, configs


def _event_value(summary: dict[str, Any], name: str) -> Any:
    for event in reversed(summary.get("events", [])):
        if event.get("name") == name:
            return event.get("value")
    return None


def _write_batch_summary(outdir: Path, payload: dict[str, Any]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "status",
        "outdir",
        "elapsed_seconds",
        "target_units",
        "verified_candidates",
        "genome_proportion_percent",
        "error",
    ]
    tsv_path = outdir / "batch_summary.tsv"
    with open(tsv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(payload["samples"])
    (outdir / "batch_summary.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )


def run_batch(
    settings: dict[str, Any],
    configs: list[dict[str, Any]],
    dry_run: bool = False,
) -> dict[str, Any]:
    """Run samples serially and update compact batch summaries after every sample."""
    outdir = Path(settings["outdir"])
    payload: dict[str, Any] = {
        "status": "RUNNING",
        "config": settings["config"],
        "samples_tsv": settings["samples_tsv"],
        "continue_on_error": settings["continue_on_error"],
        "total_samples": len(configs),
        "processed_samples": 0,
        "successful_samples": 0,
        "failed_samples": 0,
        "samples": [],
    }
    _write_batch_summary(outdir, payload)
    for config in configs:
        started = time.monotonic()
        sample = config["_batch_sample"]
        summary: dict[str, Any] = {}
        error = ""
        try:
            summary = Pipeline(config, dry_run=dry_run).run()
            status = summary["status"]
        except Exception as exc:
            status = "FAILED"
            error = str(exc)
        proportion = summary.get("genome_proportion", {})
        row = {
            "sample": sample,
            "status": status,
            "outdir": config["project"]["outdir"],
            "elapsed_seconds": f"{time.monotonic() - started:.3f}",
            "target_units": _event_value(summary, "target_units") or "",
            "verified_candidates": _event_value(summary, "verified_candidates") or "",
            "genome_proportion_percent": proportion.get("mean", ""),
            "error": error or summary.get("error", ""),
        }
        payload["samples"].append(row)
        payload["processed_samples"] += 1
        if status in {"SUCCESS", "DRY_RUN_OK"}:
            payload["successful_samples"] += 1
        else:
            payload["failed_samples"] += 1
        payload["status"] = "RUNNING"
        _write_batch_summary(outdir, payload)
        if status == "FAILED" and not settings["continue_on_error"]:
            break

    if payload["failed_samples"]:
        payload["status"] = "PARTIAL_FAILURE"
    elif dry_run:
        payload["status"] = "DRY_RUN_OK"
    else:
        payload["status"] = "SUCCESS"
    _write_batch_summary(outdir, payload)
    return payload
