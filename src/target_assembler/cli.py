from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

from . import __version__
from .batch import load_batch_configs, run_batch
from .config import VALIDATION_PROFILES, infer_layout, load_config, template_config
from .errors import TargetAssemblerError
from .pipeline import Pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="target-assembler",
        description="Recover and validate repeat-sequence candidates with a configurable MIRA workflow.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init = subparsers.add_parser("init", help="write an editable YAML configuration")
    init.add_argument("path", nargs="?", default="target-assembler.yaml")
    init.add_argument("--layout", choices=["paired", "single"], default="paired")
    init.add_argument("--profile", choices=["auto", *VALIDATION_PROFILES], default="auto")
    init.add_argument("--force", action="store_true", help="overwrite an existing configuration")

    run = subparsers.add_parser("run", help="run the pipeline")
    run.add_argument("config")
    run.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override one YAML value")
    run.add_argument("--dry-run", action="store_true", help="validate inputs and show the resolved plan")
    run.add_argument(
        "--save-intermediate-libraries",
        action="store_true",
        help="save clean, sampled, target, and expanded FASTQ libraries",
    )
    run.add_argument("--keep-work", action="store_true", help="retain the complete temporary work directory")

    batch = subparsers.add_parser("batch", help="run TSV-defined libraries sequentially")
    batch.add_argument("config")
    batch.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override one YAML value")
    batch.add_argument("--dry-run", action="store_true", help="validate every sample without running external tools")
    batch.add_argument(
        "--save-intermediate-libraries",
        action="store_true",
        help="save clean, sampled, target, and expanded FASTQ libraries for every sample",
    )
    batch.add_argument("--keep-work", action="store_true", help="retain every complete temporary work directory")

    validate = subparsers.add_parser(
        "validate", help="validate one candidate set with multiple profiles without rerunning MIRA"
    )
    validate.add_argument("config")
    validate.add_argument("--candidates", help="candidate FASTA; defaults to input.variants_fasta/output candidates")
    validate.add_argument("--reads1", help="validation R1/SE FASTQ; bypasses preprocessing")
    validate.add_argument("--reads2", help="validation R2 FASTQ")
    validate.add_argument(
        "--profiles", nargs="+", choices=list(VALIDATION_PROFILES),
        default=None, help="profiles to compare (default: lenient balanced strict)",
    )
    validate.add_argument("--max-units", type=int, help="one common library cap for a fair comparison")
    validate.add_argument("--outdir", help="comparison output directory")
    validate.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")

    check = subparsers.add_parser("check", help="validate configuration, inputs, and dependencies")
    check.add_argument("config")
    check.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")

    profiles = subparsers.add_parser("profiles", help="show validation presets")
    profiles.add_argument("--json", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "init":
            destination = Path(args.path)
            if destination.exists() and not args.force:
                raise TargetAssemblerError(f"Configuration already exists: {destination}; use --force to replace it")
            destination.write_text(
                yaml.safe_dump(template_config(args.layout, args.profile), sort_keys=False, allow_unicode=True),
                encoding="utf-8",
            )
            print(f"Created {destination}")
            return 0
        if args.command == "profiles":
            if args.json:
                print(json.dumps(VALIDATION_PROFILES, indent=2))
            else:
                for name, profile in VALIDATION_PROFILES.items():
                    print(
                        f"{name:8s} source={profile['source']:<9s} max_units={profile['max_units']:<8d} "
                        f"coverage>={profile['minimum_average_coverage']:g}x "
                        f"breadth>={profile['minimum_covered_percent']:g}%"
                    )
            return 0

        overrides = list(args.set)
        if getattr(args, "save_intermediate_libraries", False):
            overrides.append("outputs.save_intermediate_libraries=true")
        if getattr(args, "keep_work", False):
            overrides.append("runtime.keep_work=true")
        if args.command == "batch":
            settings, configs = load_batch_configs(args.config, overrides)
            summary = run_batch(settings, configs, dry_run=args.dry_run)
            print(json.dumps(summary, indent=2, ensure_ascii=False))
            return 0 if summary["status"] in {"SUCCESS", "DRY_RUN_OK"} else 2
        if args.command == "validate":
            overrides.append("project.mode=validate")
            if args.candidates:
                overrides.append(f"input.variants_fasta={json.dumps(args.candidates)}")
            if args.reads1:
                overrides.append(f"input.validation_reads1={json.dumps(args.reads1)}")
            if args.reads2:
                overrides.append(f"input.validation_reads2={json.dumps(args.reads2)}")
            if args.profiles:
                overrides.append(f"validation.compare_profiles={json.dumps(args.profiles)}")
            if args.max_units is not None:
                overrides.append(f"validation.comparison_max_units={args.max_units}")
            if args.outdir:
                overrides.append(f"project.outdir={json.dumps(args.outdir)}")
        config = load_config(args.config, overrides)
        pipeline = Pipeline(config, dry_run=getattr(args, "dry_run", False))
        if args.command == "check":
            try:
                pipeline.preflight(require_programs=True)
                print(f"Configuration OK ({infer_layout(config)} library)")
                print("Dependencies: " + ", ".join(pipeline.dependencies()))
            finally:
                pipeline.cleanup()
            return 0
        summary = pipeline.run()
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return 0
    except (TargetAssemblerError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
