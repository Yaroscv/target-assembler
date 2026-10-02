from __future__ import annotations

import copy
import re
from pathlib import Path
from typing import Any

import yaml

from .errors import ConfigError


DEFAULTS: dict[str, Any] = {
    "project": {
        "name": "target_assembly",
        "mode": "all",
        "outdir": "target_assembler_output",
        "temp_dir": None,
    },
    "input": {
        "reads1": None,
        "reads2": None,
        "layout": "auto",
        "hmm": None,
        "clean_reads1": None,
        "clean_reads2": None,
        "sample_reads1": None,
        "sample_reads2": None,
        "variants_fasta": None,
        "validation_reads1": None,
        "validation_reads2": None,
    },
    "runtime": {
        "threads": 8,
        "seed": 42,
        "keep_work": False,
        "temp_space_check": True,
        "temp_space_safety_factor": 1.20,
        "temp_space_margin_gb": 1.0,
    },
    "batch": {
        "samples_tsv": None,
        "continue_on_error": True,
    },
    "preprocessing": {
        "enabled": True,
        "min_length": 50,
        "cut_tail": True,
        "detect_adapters": True,
        "extra_args": [],
    },
    "sampling": {
        "mode": "adaptive",
        "fixed_units": 0,
        "micro_units": 500000,
        "target_hmm_hits": 1000,
        "minimum_units": 100000,
        "maximum_units": None,
        "zero_hit_fallback_units": 1000000,
    },
    "selection": {
        # HMM recruitment stringency. nhmmer's own default is -E 10, which for a
        # conserved rRNA profile recruits across kingdoms: an 18S/5.8S/28S model
        # built from a plant 35S unit scores fungal rDNA highly (measured: 52 of
        # 60 read-length windows of a fungal rDNA contig scored above the default
        # threshold against such a profile), so an unthresholded run pulls those
        # reads into the target pool and assembles them as candidate variants. Both keys are applied to
        # nhmmer AND re-applied after parsing the table, so a threshold is
        # enforced even for a table written by an earlier run.
        # null on both = previous behaviour (nhmmer default).
        "hmm_evalue": None,      # -E  : report/keep hits at or below this E-value
        "hmm_score_min": None,   # -T  : keep hits at or above this bit score
        "flank_bp": 100,
        "minimum_flank_length": 60,
        "blast_identity": 80.0,
        "blast_min_query_coverage": 0.0,
        "blast_evalue": 1e-5,
        "blast_max_targets": 1000000,
    },
    "assembly": {
        "technology": "solexa",
        "rounds": 2,
        "mira_max_units": 30000,
        "contig_min_length": 250,
        "insert_size_min": 200,
        "insert_size_max": 400,
        "segment_placement": "---> <---",
        "job": "denovo,est,accurate",
        "parameters": ["COMMON_SETTINGS", "-NW:cmrnl=no", "-NW:ctp=no", "-NW:cac=no"],
        "fallback_to_round1": True,
    },
    "expansion": {
        "enabled": True,
        "alignment": "local",
        "sensitivity": "very-sensitive",
        "minimum_mapq": 0,
        "max_alignments": 1,
        # Largest soft/hard clip allowed where the reference contig CONTINUES
        # (bp). null = no test, the upstream behaviour. See pipeline.internal_clip.
        "max_internal_clip": None,
    },
    "extraction": {
        "mode": "tandem",
        "flank5": 0,
        "flank3": 0,
        "minimum_model_fraction": 0.90,
    },
    "validation": {
        "profile": "auto",
        "reference_mode": "tandem",
        "tandem_min_copies": 3,
        "compare_profiles": ["lenient", "balanced", "strict"],
        "comparison_max_units": 1000000,
        "kmer_size": 27,
        "window_size": 50,
        "softclip_minimum_bases": 12,
        "auto": {
            "balanced_min_target_units": 100,
            "strict_min_target_units": 1000,
        },
    },
    "proportion": {
        "enabled": True,
        "sample_units": 1000000,
        "replicates": 1,
        "report_metric": "aligned_nt",
        "include_envelope_metric": True,
    },
    "outputs": {"save_intermediate_libraries": False},
}


VALIDATION_COMMON: dict[str, Any] = {
    "source": "expanded",
    "max_units": 250000,
    "alignment": "local",
    "sensitivity": "very-sensitive",
    "max_alignments": 10,
    "primary_only": True,
    "minimum_mapq": 0,
    "require_unique_kmers": False,
    "reference_mode": "tandem",
    "tandem_min_copies": 3,
}


VALIDATION_PROFILES: dict[str, dict[str, Any]] = {
    "lenient": {
        "source": "expanded",
        "max_units": 100000,
        "minimum_average_coverage": 2.0,
        "minimum_covered_percent": 40.0,
        "minimum_mapped_fragments": 3,
        "pair_rule": "disabled",
        "minimum_proper_pairs": 0,
        "minimum_proper_pair_percent": 0.0,
        "maximum_zero_run": 150,
        "maximum_adjacent_window_fold": 100.0,
        "maximum_softclip_percent": 100.0,
        "minimum_unique_kmer_fraction": 0.0,
    },
    "balanced": {
        "source": "expanded",
        "max_units": 250000,
        "minimum_average_coverage": 3.0,
        "minimum_covered_percent": 70.0,
        "minimum_mapped_fragments": 8,
        "pair_rule": "any",
        "minimum_proper_pairs": 5,
        "minimum_proper_pair_percent": 5.0,
        "maximum_zero_run": 75,
        "maximum_adjacent_window_fold": 15.0,
        "maximum_softclip_percent": 60.0,
        "minimum_unique_kmer_fraction": 0.10,
    },
    "strict": {
        "source": "clean",
        "max_units": 1000000,
        "minimum_average_coverage": 5.0,
        "minimum_covered_percent": 90.0,
        "minimum_mapped_fragments": 20,
        "pair_rule": "any",
        "minimum_proper_pairs": 10,
        "minimum_proper_pair_percent": 10.0,
        "maximum_zero_run": 25,
        "maximum_adjacent_window_fold": 5.0,
        "maximum_softclip_percent": 35.0,
        "minimum_unique_kmer_fraction": 0.25,
    },
}


def deep_update(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = copy.deepcopy(value)
    return base


def read_config_mapping(path: str | Path) -> tuple[Path, dict[str, Any]]:
    """Read one YAML file without applying defaults or scientific validation."""
    config_path = Path(path).expanduser().resolve()
    try:
        raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    except OSError as exc:
        raise ConfigError(f"Cannot read configuration: {config_path}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML in {config_path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError("The YAML root must be a mapping.")
    return config_path, raw


def set_dotted_value(mapping: dict[str, Any], dotted: str, value: Any) -> None:
    if not dotted or dotted.startswith("_"):
        raise ConfigError(f"Invalid configuration key: {dotted!r}")
    cursor = mapping
    keys = dotted.split(".")
    for key in keys[:-1]:
        cursor = cursor.setdefault(key, {})
        if not isinstance(cursor, dict):
            raise ConfigError(f"Cannot set nested value beneath {key!r}")
    cursor[keys[-1]] = value


def apply_overrides(mapping: dict[str, Any], overrides: list[str] | None = None) -> None:
    for item in overrides or []:
        if "=" not in item:
            raise ConfigError(f"Override must be KEY=VALUE: {item}")
        dotted, text = item.split("=", 1)
        value = yaml.safe_load(text)
        set_dotted_value(mapping, dotted, value)


def config_from_mapping(raw: dict[str, Any], config_path: str | Path) -> dict[str, Any]:
    """Resolve defaults and paths for an in-memory YAML mapping."""
    config_path = Path(config_path).expanduser().resolve()

    merged = deep_update(copy.deepcopy(DEFAULTS), raw)
    merged["_config_path"] = str(config_path)
    merged["_config_dir"] = str(config_path.parent)
    merged["_validation_overrides"] = copy.deepcopy(raw.get("validation", {}))
    _resolve_paths(merged)
    validate_config(merged)
    return merged


def load_config(path: str | Path, overrides: list[str] | None = None) -> dict[str, Any]:
    config_path, raw = read_config_mapping(path)
    apply_overrides(raw, overrides)
    return config_from_mapping(raw, config_path)


def _resolve_paths(config: dict[str, Any]) -> None:
    base = Path(config["_config_dir"])
    for key in (
        "reads1", "reads2", "hmm", "clean_reads1", "clean_reads2",
        "sample_reads1", "sample_reads2", "variants_fasta",
        "validation_reads1", "validation_reads2",
    ):
        value = config["input"].get(key)
        if value:
            path = Path(value).expanduser()
            config["input"][key] = str(path if path.is_absolute() else (base / path).resolve())
    for section, key in (("project", "outdir"), ("project", "temp_dir")):
        value = config[section].get(key)
        if value:
            path = Path(value).expanduser()
            config[section][key] = str(path if path.is_absolute() else (base / path).resolve())
    samples_tsv = config["batch"].get("samples_tsv")
    if samples_tsv:
        path = Path(samples_tsv).expanduser()
        config["batch"]["samples_tsv"] = str(path if path.is_absolute() else (base / path).resolve())


def infer_layout(config: dict[str, Any]) -> str:
    requested = config["input"]["layout"]
    if requested == "auto":
        has_r2 = bool(
            config["input"].get("reads2")
            or config["input"].get("clean_reads2")
            or config["input"].get("sample_reads2")
            or config["input"].get("validation_reads2")
        )
        return "paired" if has_r2 else "single"
    return requested


def validate_config(config: dict[str, Any]) -> None:
    mode = config["project"]["mode"]
    if mode not in {"all", "variants", "proportion", "validate"}:
        raise ConfigError("project.mode must be all, variants, proportion, or validate")
    if config["input"]["layout"] not in {"auto", "paired", "single"}:
        raise ConfigError("input.layout must be auto, paired, or single")
    layout = infer_layout(config)
    inp = config["input"]
    if mode in {"all", "variants"} and not inp.get("hmm"):
        raise ConfigError("input.hmm is required")
    if layout == "paired":
        pairs = (
            (inp.get("reads1"), inp.get("reads2")),
            (inp.get("clean_reads1"), inp.get("clean_reads2")),
            (inp.get("sample_reads1"), inp.get("sample_reads2")),
            (inp.get("validation_reads1"), inp.get("validation_reads2")),
        )
        for first, second in pairs:
            if bool(first) != bool(second):
                raise ConfigError("Paired input paths must always be supplied as R1/R2 pairs")
    elif inp.get("reads2") or inp.get("clean_reads2") or inp.get("sample_reads2") or inp.get("validation_reads2"):
        raise ConfigError("R2 paths cannot be used with input.layout: single")

    has_raw = bool(inp.get("reads1"))
    has_clean = bool(inp.get("clean_reads1"))
    has_sample = bool(inp.get("sample_reads1"))
    has_validation = bool(inp.get("validation_reads1"))
    # A cleaned library WITHOUT a sample is a supported state: prepare_libraries
    # takes clean_reads1/2 when present and, with no sample_reads1/2, builds the
    # assembly sample from it with make_assembly_sample(). That is the natural way
    # to re-run with a different HMM - preprocessing is profile-independent and
    # expensive, adaptive sampling is profile-DEPENDENT and cheap, so the sample
    # must be rebuilt while the cleaned library is reused. Requiring a sample
    # alongside clean made the check stricter than the pipeline and rejected it.
    if mode in {"all", "variants"} and not (has_raw or has_clean):
        raise ConfigError(
            "Variant assembly needs raw reads, or a cleaned library "
            "(input.clean_reads1/2, with input.sample_reads1/2 optional)"
        )
    if mode == "proportion" and not (has_raw or has_clean):
        raise ConfigError("Proportion mode needs raw or clean reads")
    if mode == "validate" and not (has_validation or has_clean or has_raw):
        raise ConfigError("Validation mode needs validation, clean, or raw reads")

    profile = config["validation"]["profile"]
    if profile not in {"auto", *VALIDATION_PROFILES}:
        raise ConfigError("validation.profile must be auto, lenient, balanced, or strict")
    if "source" in config["validation"] and config["validation"]["source"] not in {"clean", "expanded", "assembly"}:
        raise ConfigError("validation.source must be clean, expanded, or assembly")
    if "pair_rule" in config["validation"] and config["validation"]["pair_rule"] not in {"any", "all", "disabled"}:
        raise ConfigError("validation.pair_rule must be any, all, or disabled")
    if config["validation"]["reference_mode"] not in {"linear", "tandem"}:
        raise ConfigError("validation.reference_mode must be linear or tandem")
    compare_profiles = config["validation"].get("compare_profiles", [])
    if not isinstance(compare_profiles, list) or not compare_profiles:
        raise ConfigError("validation.compare_profiles must be a non-empty list")
    if any(profile not in VALIDATION_PROFILES for profile in compare_profiles):
        raise ConfigError("validation.compare_profiles may contain lenient, balanced, and strict")
    technology = str(config["assembly"]["technology"])
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+", technology):
        raise ConfigError("assembly.technology contains unsafe characters")
    if config["assembly"]["rounds"] not in {1, 2}:
        raise ConfigError("assembly.rounds must be 1 or 2")
    if config["extraction"]["mode"] not in {"tandem", "site"}:
        raise ConfigError("extraction.mode must be tandem or site")
    sel = config["selection"]
    if sel.get("hmm_evalue") is not None and float(sel["hmm_evalue"]) <= 0:
        raise ConfigError("selection.hmm_evalue must be greater than 0, or null")
    if sel.get("hmm_score_min") is not None and float(sel["hmm_score_min"]) < 0:
        raise ConfigError("selection.hmm_score_min must be 0 or greater, or null")
    if config["sampling"]["mode"] not in {"adaptive", "fixed", "all"}:
        raise ConfigError("sampling.mode must be adaptive, fixed, or all")
    if config["proportion"]["report_metric"] not in {"aligned_nt", "envelope_nt"}:
        raise ConfigError("proportion.report_metric must be aligned_nt or envelope_nt")
    if int(config["runtime"]["threads"]) < 1:
        raise ConfigError("runtime.threads must be at least 1")
    if float(config["runtime"]["temp_space_safety_factor"]) < 1.0:
        raise ConfigError("runtime.temp_space_safety_factor must be at least 1.0")
    if float(config["runtime"]["temp_space_margin_gb"]) < 0.0:
        raise ConfigError("runtime.temp_space_margin_gb cannot be negative")
    if int(config["proportion"]["sample_units"]) < 1:
        raise ConfigError("proportion.sample_units must be at least 1")
    if int(config["proportion"]["replicates"]) < 1:
        raise ConfigError("proportion.replicates must be at least 1")
    if int(config["validation"]["kmer_size"]) < 1 or int(config["validation"]["window_size"]) < 1:
        raise ConfigError("validation.kmer_size and validation.window_size must be positive")
    if int(config["validation"]["tandem_min_copies"]) < 3:
        raise ConfigError("validation.tandem_min_copies must be at least 3")
    if int(config["validation"]["comparison_max_units"] or 0) < 0:
        raise ConfigError("validation.comparison_max_units cannot be negative")
    fraction = float(config["extraction"]["minimum_model_fraction"])
    if not 0.0 <= fraction <= 1.0:
        raise ConfigError("extraction.minimum_model_fraction must be between 0 and 1")
    if int(config["assembly"]["mira_max_units"] or 0) < 0:
        raise ConfigError("assembly.mira_max_units cannot be negative")
    clip = config.get("expansion", {}).get("max_internal_clip")
    if clip is not None and (isinstance(clip, bool) or int(clip) < 0):
        raise ConfigError("expansion.max_internal_clip must be a non-negative integer or null")
    if int(config["assembly"]["insert_size_min"]) > int(config["assembly"]["insert_size_max"]):
        raise ConfigError("assembly.insert_size_min cannot exceed insert_size_max")
    if not isinstance(config["outputs"]["save_intermediate_libraries"], bool):
        raise ConfigError("outputs.save_intermediate_libraries must be true or false")
    if not isinstance(config["runtime"]["keep_work"], bool):
        raise ConfigError("runtime.keep_work must be true or false")
    if not isinstance(config["proportion"]["include_envelope_metric"], bool):
        raise ConfigError("proportion.include_envelope_metric must be true or false")


def choose_validation_profile(config: dict[str, Any], layout: str, target_units: int) -> str:
    requested = config["validation"]["profile"]
    if requested != "auto":
        return requested
    if layout == "single":
        return "lenient"
    auto = config["validation"]["auto"]
    if target_units >= int(auto["strict_min_target_units"]):
        return "strict"
    if target_units >= int(auto["balanced_min_target_units"]):
        return "balanced"
    return "lenient"


def resolve_validation(config: dict[str, Any], chosen_profile: str) -> dict[str, Any]:
    result = copy.deepcopy(VALIDATION_COMMON)
    deep_update(result, VALIDATION_PROFILES[chosen_profile])
    explicit = copy.deepcopy(config.get("_validation_overrides", {}))
    explicit.pop("profile", None)
    explicit.pop("auto", None)
    explicit.pop("compare_profiles", None)
    explicit.pop("comparison_max_units", None)
    deep_update(result, explicit)
    result["profile"] = chosen_profile
    for key in ("kmer_size", "window_size", "softclip_minimum_bases"):
        if key not in result:
            result[key] = config["validation"][key]
    if result["source"] not in {"clean", "expanded", "assembly"}:
        raise ConfigError("validation.source must be clean, expanded, or assembly")
    if result["pair_rule"] not in {"any", "all", "disabled"}:
        raise ConfigError("validation.pair_rule must be any, all, or disabled")
    if int(result["max_alignments"]) < 1:
        raise ConfigError("validation.max_alignments must be at least 1")
    if result.get("max_units") is not None and int(result["max_units"]) < 0:
        raise ConfigError("validation.max_units cannot be negative")
    return result


def template_config(layout: str = "paired", profile: str = "auto") -> dict[str, Any]:
    paired = layout == "paired"
    return {
        "project": {
            "name": "my_repeat",
            "mode": "all",
            "outdir": "results",
            "temp_dir": None,
        },
        "input": {
            "reads1": "reads_R1.fastq.gz" if paired else "reads.fastq.gz",
            "reads2": "reads_R2.fastq.gz" if paired else None,
            "layout": layout,
            "hmm": "target.hmm",
        },
        "runtime": {
            "threads": 8,
            "seed": 42,
            "keep_work": False,
            "temp_space_check": True,
            "temp_space_margin_gb": 1.0,
        },
        "sampling": {"mode": "adaptive", "target_hmm_hits": 1000},
        "assembly": {
            "technology": "solexa",
            "mira_max_units": 30000,
            "insert_size_min": 200,
            "insert_size_max": 400,
        },
        "extraction": {"mode": "tandem", "minimum_model_fraction": 0.90},
        "validation": {"profile": profile, "reference_mode": "tandem"},
        "proportion": {
            "enabled": True,
            "sample_units": 1000000,
            "replicates": 1,
            "include_envelope_metric": True,
        },
        "outputs": {"save_intermediate_libraries": False},
    }
