from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

from target_assembler.batch import load_batch_configs, run_batch
from target_assembler.errors import ConfigError
from target_assembler.pipeline import Pipeline


def write_fastq(path: Path, mate: int) -> None:
    path.write_text(f"@read1/{mate}\nACGTACGT\n+\nIIIIIIII\n", encoding="utf-8")


def write_batch_fixture(directory: Path) -> Path:
    sheet_dir = directory / "sheet"
    reads_dir = sheet_dir / "reads"
    reads_dir.mkdir(parents=True)
    (directory / "target.hmm").write_text("HMMER3/f\n", encoding="utf-8")
    for sample in ("alpha", "beta"):
        write_fastq(reads_dir / f"{sample}_R1.fastq", 1)
        write_fastq(reads_dir / f"{sample}_R2.fastq", 2)
    (sheet_dir / "samples.tsv").write_text(
        "sample\treads1\treads2\tassembly.insert_size_min\tvalidation.profile\n"
        "alpha\treads/alpha_R1.fastq\treads/alpha_R2.fastq\t150\tbalanced\n"
        "beta\treads/beta_R1.fastq\treads/beta_R2.fastq\t200\tstrict\n",
        encoding="utf-8",
    )
    config = {
        "project": {"mode": "variants", "outdir": "results"},
        "batch": {"samples_tsv": "sheet/samples.tsv", "continue_on_error": True},
        "input": {"layout": "paired", "hmm": "target.hmm"},
        "proportion": {"enabled": False},
    }
    path = directory / "batch.yaml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    return path


class BatchTests(unittest.TestCase):
    def test_batch_rows_inherit_yaml_and_resolve_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            settings, configs = load_batch_configs(
                write_batch_fixture(directory),
                ["outputs.save_intermediate_libraries=true"],
            )
            self.assertEqual(len(configs), 2)
            self.assertEqual(configs[0]["project"]["name"], "alpha")
            self.assertEqual(Path(configs[0]["project"]["outdir"]), directory / "results" / "alpha")
            self.assertEqual(
                Path(configs[0]["input"]["reads1"]),
                directory / "sheet" / "reads" / "alpha_R1.fastq",
            )
            self.assertEqual(configs[0]["assembly"]["insert_size_min"], 150)
            self.assertEqual(configs[1]["validation"]["profile"], "strict")
            self.assertTrue(configs[1]["outputs"]["save_intermediate_libraries"])
            self.assertEqual(Path(settings["outdir"]), directory / "results")

    def test_duplicate_samples_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config_path = write_batch_fixture(directory)
            table = directory / "sheet" / "samples.tsv"
            table.write_text(
                "sample\treads1\treads2\n"
                "same\treads/alpha_R1.fastq\treads/alpha_R2.fastq\n"
                "same\treads/beta_R1.fastq\treads/beta_R2.fastq\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ConfigError, "Duplicate sample"):
                load_batch_configs(config_path)

    def test_batch_summary_is_written(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            settings, configs = load_batch_configs(write_batch_fixture(directory))
            fake_summary = {
                "status": "SUCCESS",
                "events": [
                    {"name": "target_units", "value": 12},
                    {"name": "verified_candidates", "value": 3},
                ],
                "genome_proportion": {"mean": 0.5},
            }
            with patch("target_assembler.batch.Pipeline") as mocked:
                mocked.return_value.run.return_value = fake_summary
                summary = run_batch(settings, configs)
            self.assertEqual(summary["status"], "SUCCESS")
            with open(directory / "results" / "batch_summary.tsv", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["sample"] for row in rows], ["alpha", "beta"])
            self.assertEqual(rows[0]["target_units"], "12")

    def test_mira_work_is_inside_ephemeral_workdir(self):
        with tempfile.TemporaryDirectory() as temporary:
            _, configs = load_batch_configs(write_batch_fixture(Path(temporary)))
            pipeline = Pipeline(configs[0], dry_run=True)
            workdir = pipeline.workdir
            self.assertEqual(pipeline.mira_dir.parent, workdir)
            summary = pipeline.run()
            self.assertEqual(summary["status"], "DRY_RUN_OK")
            self.assertFalse(workdir.exists())

    def test_intermediate_library_retention_is_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            _, configs = load_batch_configs(write_batch_fixture(directory))
            self.assertFalse(configs[0]["outputs"]["save_intermediate_libraries"])
            pipeline = Pipeline(configs[0], dry_run=True)
            pipeline._save_library(
                (configs[0]["input"]["reads1"], configs[0]["input"]["reads2"]),
                "target",
            )
            saved = Path(configs[0]["project"]["outdir"]) / "intermediate_libraries"
            self.assertTrue((saved / "target_R1.fastq").exists())
            self.assertTrue((saved / "target_R2.fastq").exists())
            pipeline.cleanup()


if __name__ == "__main__":
    unittest.main()
