# target-assembler

`target-assembler` is a configurable Python workflow for recovering a **set of candidate repeat variants** from short-read libraries. It retains the two-round MIRA strategy that performed well for rDNA while making read layout, MIRA technology, sampling, extraction, validation, and genome-proportion estimation explicit in YAML.

The intended result is not a single majority consensus. The pipeline preserves distinct assembled candidates, their provenance, and the evidence supporting each candidate.

Current release: **0.1.0+specieshunter.5** - upstream 0.1.0 plus the five commits
of branch `species-hunter` below. The `+specieshunter.5` part is a PEP 440 *local
version label*: it marks a build of 0.1.0 with local changes, sorts after 0.1.0,
and can never collide with a future upstream release number (0.1.1, 0.2.0, ...),
which remains the upstream author's to choose.

## This fork: branch `species-hunter`

This branch is upstream **v0.1.0 plus five commits**, written while running
target-assembler on 35S rDNA of *Fraxinus* (ash) whole-genome libraries in the
[species-hunter](https://github.com/Yaroscv/species-hunter) workflow. Each commit
adds one configuration key or fixes one failure; **every new key defaults to
`null`, which reproduces the upstream behaviour exactly**, so an existing YAML runs
unchanged. Each commit carries its own tests (32 tests in total pass with
`PYTHONPATH=src python -m unittest discover -s tests -v`).

| # | commit | new key / behaviour | default | recorded as | why it was needed |
|---|---|---|---|---|---|
| 1 | HMM recruitment stringency | `selection.hmm_score_min` (`-T`), `selection.hmm_evalue` (`-E`), applied at every `nhmmer` call and re-applied after the table is parsed | `null` (nhmmer's own default, `-E 10`) | event `hmm_hits_below_threshold` | a plant 35S profile scores fungal rDNA: 52 of 60 fungal 150 bp windows passed the default threshold (section *HMM recruitment stringency*) |
| 2 | no self-copy of a resumed library | when the library being saved is the file the run read, it is left in place | - | event `saved_<prefix>_library_already_in_place` | a run resumed from `intermediate_libraries/` failed with `SameFileError` |
| 3 | cleaned library without a sample | `input.clean_reads1/2` accepted without `input.sample_reads1/2`; sampling is rebuilt | - | - | reuse the expensive cleaned library with a different HMM; a sample sized for one profile under-recruits for another |
| 4 | round-2 internal-clip filter | `expansion.max_internal_clip` (bp): reject a local alignment clipped where the contig continues | `null` (no test) | event `expansion_internal_clip_rejected` | reads sharing only the universal 5.8S block (foreign rDNA) were recruited and broke MIRA's contigs at the core (section *Round-2 expansion*) |
| 5 | MIRA's own thread count | `assembly.mira_threads`, passed as `-GE:not` after `parameters` | `null` (= `runtime.threads`) | - | MIRA on several threads returned different contig sets from identical input; one thread is repeatable |

Install the branch:

```bash
git clone -b species-hunter https://github.com/Yaroscv/target-assembler.git
cd target-assembler && python -m pip install -e .
target-assembler --version                                    # target-assembler 0.1.0+specieshunter.5
python -c "from target_assembler.config import DEFAULTS as D; print(D['selection']['hmm_score_min'], D['expansion']['max_internal_clip'], D['assembly']['mira_threads'])"   # None None None
```

The settings used with these keys for a 35S rDNA target, and why, are in
species-hunter's `params/ta_35S.yaml`; the notes below (*MIRA settings for rDNA
targets*, *Tandem extraction completeness*) record what those runs showed about
the upstream options themselves.

## Capabilities

- paired-end and single-end libraries;
- fixed, adaptive, or full-library sampling;
- target selection with HMM and flank-based recruitment;
- one- or two-round MIRA assembly with configurable sequencing technology;
- tandem-repeat or single-site extraction;
- `lenient`, `balanced`, `strict`, and data-aware `auto` validation profiles;
- tandem-aware validation using coverage, PE support, soft clipping, multimapping categories, and variant-specific k-mers;
- genome-proportion estimation with deterministic sampling replicates;
- sequential processing of many libraries from a TSV sample sheet;
- automatic removal of heavy working files after each library, with explicit opt-in retention.

## Scientific interpretation

The main FASTA files deliberately distinguish three levels:

1. `Repeats_extracted_raw.fasta` — all extracted repeat instances.
2. `Repeats_candidates.fasta` — exactly dereplicated assembled candidates.
3. `Repeats_verified.fasta` — candidates passing the selected evidence profile.

Candidates that do not pass are retained in `Repeats_rejected.fasta` and `variants.tsv`; a low-copy or small library therefore does not silently lose all potentially useful sequences. A short-read assembly can still be unphased across distances longer than the physical library insert. The scientifically precise term is therefore **assembled candidate ribotype/repeat variant**, particularly for long ITS or IGS regions.

## Requirements

Python 3.10 or later and PyYAML are required. Depending on the selected mode, the following command-line programs must be on `PATH`:

- `fastp`
- HMMER: `nhmmer`, `hmmstat`, `hmmbuild`
- BLAST+: `makeblastdb`, `blastn`
- `mira`
- `bowtie2`, `bowtie2-build`
- `mafft`

## Installation

From the project directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
target-assembler --help
```

Install the external bioinformatics programs separately using your preferred Conda/Mamba or system setup.

## Quick start

Create a paired-end configuration:

```bash
target-assembler init analysis.yaml --layout paired --profile auto
```

Edit the paths and parameters, then check and run it:

```bash
target-assembler check analysis.yaml
target-assembler run analysis.yaml
```

Override one value without editing YAML:

```bash
target-assembler run analysis.yaml --set assembly.technology=solexa \
  --set validation.profile=balanced
```

Ready-to-edit examples are in `examples/`:

- `paired-balanced.yaml`
- `paired-strict.yaml`
- `single-lenient.yaml`
- `proportion-only.yaml`
- `validate-profiles.yaml`
- `batch.yaml` with `samples.tsv`

### HMM recruitment stringency

`nhmmer` is invoked with its own default reporting threshold unless a cutoff is
configured, and that default is `-E 10`. For a profile built over a ribosomal
repeat this is rarely what you want: 18S, 5.8S and 28S are conserved across
eukaryotes, so a profile built from a plant 35S unit scores **fungal** rDNA
highly. Measured, not assumed: shredding a fungal rDNA contig into 150 bp
windows and scanning it with such a profile put **52 of 60 windows above the
default reporting threshold**, median bit score 66 and p90 108. In a
genome-skimming library that contains any fungal rDNA at all, an unthresholded
recruitment therefore pulls those reads into the target pool, MIRA assembles
them, and they arrive in `Repeats_candidates.fasta` as candidate repeat variants
— and `proportion` measures them.

```yaml
selection:
  hmm_score_min: 98.4   # -T, a bit score; null = nhmmer's default
  hmm_evalue: null      # -E, an E-value; null = nhmmer's default
```

Set one, not both: `-T` and `-E` are mutually exclusive in HMMER, so a score
threshold replaces the E-value model entirely and `hmm_score_min` takes
precedence when both are given. Both are also re-applied after the table is
parsed, so a threshold holds even for a `--tblout` written by an earlier run, and
the number of hits dropped is recorded as the `hmm_hits_below_threshold` event in
the run summary rather than passing silently.

A bit score threshold is only meaningful when it was **measured on the profile in
question**, at the read length recruitment actually sees: nhmmer scores
accumulate over the matching region, so a cutoff calibrated against multi-kb
queries does not transfer to 150 bp reads. Shred a positive and a negative
control into read-length windows, compare the score distributions, and take the
threshold from the gap between them.

### Round-2 expansion: reads clipped inside a contig

Round 2 recruits every read of the clean library that bowtie2 aligns locally to a
round-1 contig. A *local* alignment can stop anywhere, so a read that shares only
a conserved stretch with a contig - foreign rDNA on the universal 5.8S block, for
example - is recruited too, with the rest of it soft-clipped. MIRA then sees
several different continuations out of the same block, which is a repeat signal.

`expansion.max_internal_clip` (bp, default `null` = no test) rejects an alignment
whose clip exceeds the limit **where the contig still continues**. Clipping where
the read runs off a contig end is never counted, so reads that extend a contig -
the point of round 2 - are kept.

```yaml
expansion:
  max_internal_clip: 20    # null keeps the upstream behaviour
```

The number of rejected alignments is recorded as the
`expansion_internal_clip_rejected` event.

### Tandem extraction completeness

In `extraction.mode: tandem`, a unit is cut from the end of one HMM hit to the end of the next hit in repeat orientation. `extraction.minimum_model_fraction` specifies how far the downstream hit must reach toward the 3-prime end of the HMM model (default `0.90`). The first hit defines only the left boundary and may be partial at a contig end; it is therefore not required to span the same fraction of the model. Set the value to `0` only to disable this downstream-completeness check explicitly.

Tandem mode therefore needs **two consecutive hits on one contig**, i.e. a contig
longer than about two repeat periods. Consecutive hits are one repeat period
apart, which is a property of the genome, not of the profile: for 5S rDNA (period
300-500 bp) a 2-3 kb contig holds several copies, which is why the 5S examples
work; for 35S rDNA the period is at least the unit length (~6.4 kb), so tandem
extraction needs a contig of ~13 kb or more. With a profile over the whole 35S unit
a contig of one complete unit plus
part of its spacer has a single hit, and the run reports `tandem_hit_pairs: 0`
with nothing extracted, although the unit itself was assembled (6,778 bp contig
covering model 2-6366 on an ash library). For such targets use
`extraction.mode: site`, which extracts one hit covering at least
`minimum_model_fraction` of the model. A profile over a short part of the unit
(ITS only) does not help tandem mode either: consecutive hits are still one
repeat period apart.

## Sequential batch analysis

The `batch` command combines shared settings from YAML with sample-specific values from a TSV table and processes samples one at a time. This keeps CPU and temporary-space demand close to that of a single run.

Minimal YAML:

```yaml
project:
  mode: all
  outdir: batch_results
  temp_dir: /mnt/d/target_tmp

batch:
  samples_tsv: samples.tsv
  continue_on_error: true

input:
  layout: paired
  hmm: 5S_profile.hmm

runtime:
  threads: 16
```

Minimal TSV:

```text
sample	reads1	reads2	assembly.insert_size_min	assembly.insert_size_max	validation.profile
A_variegatum	reads/A_var_R1.fastq.gz	reads/A_var_R2.fastq.gz	150	450	balanced
A_napellus	reads/A_nap_R1.fastq.gz	reads/A_nap_R2.fastq.gz	180	500	strict
```

`sample` is required and must be unique. Blank cells inherit the YAML value. Common columns have short aliases (`reads1`, `reads2`, `layout`, `hmm`, `mode`, `outdir`, `threads`, `seed`, and `validation_profile`); any setting can be supplied with its dotted YAML name, such as `sampling.maximum_units` or `assembly.mira_max_units`. File paths in the TSV are resolved relative to the TSV file.

Each sample is written to `<project.outdir>/<sample>/`. Compact `batch_summary.tsv` and `batch_summary.json` files are updated after every sample, so completed results remain recorded if a later library fails.

Check every row without invoking external programs, then start the serial run:

```bash
target-assembler batch batch.yaml --dry-run
target-assembler batch batch.yaml
```

## Validation philosophy

| Profile | Default evidence source | Default cap | Intended use |
|---|---:|---:|---|
| `lenient` | expanded target pool | 100,000 units | SE, small libraries, or rare repeats; candidates remain explicitly low-confidence |
| `balanced` | expanded target pool | 250,000 units | routine PE genome-skimming data; useful speed/evidence compromise |
| `strict` | cleaned library | 1,000,000 units | well-supported PE libraries; stricter MAPQ, breadth, uniformity, pair, and k-mer requirements |
| `auto` | chosen from evidence | profile-dependent | SE → lenient; PE with ≥100 selected target units → balanced; PE with ≥1,000 → strict |

A “unit” means one read pair for PE data or one read for SE data. The limits are configurable, and `max_units: null` can be used when a full validation scan is genuinely desired.

The decision for each candidate uses:

- mean depth and breadth of coverage;
- number of mapped fragments;
- longest uncovered segment;
- abrupt depth changes across windows;
- excessive soft clipping;
- proper-pair support for PE only;
- support for k-mers unique among the candidate set.

Alignment support is intentionally split in `variants.tsv`:

- `mapped_fragments` — general support, including alignments shared by similar candidates;
- `primary_mapped_fragments` — fragments contributing the primary coverage track;
- `single_variant_mapped_fragments` — fragments whose reported alignments reach only one candidate;
- `cross_variant_multimapped_fragments` — fragments reported for more than one candidate;
- `unique_kmers_supported` and `unique_kmer_fraction` — sequence evidence specific to one candidate in the current candidate set.

Bowtie2 reports up to `validation.max_alignments` alignments per read (default 10), so the two mapping categories describe the reported alignment set rather than proving global uniqueness.

### Tandem-aware validation

The default `validation.reference_mode: tandem` writes enough artificial copies of each repeat to provide internal mapping context longer than the configured PE insert. Coverage coordinates are then reduced modulo the candidate length and reported for one repeat unit. Proper pairs and soft clipping are therefore assessed away from the artificial boundary of a single-copy linear reference.

For genuinely linear targets, opt out explicitly:

```yaml
validation:
  reference_mode: linear
```

These metrics detect weak support and several common chimera warning patterns, but they do not prove complete molecular phasing beyond the insert span.

### Custom thresholds

Profiles are starting points. Any value can be overridden:

```yaml
validation:
  profile: balanced
  reference_mode: tandem # tandem | linear
  tandem_min_copies: 3
  source: expanded       # clean | expanded | assembly
  max_units: 300000      # null = all available units
  minimum_average_coverage: 3
  minimum_covered_percent: 75
  minimum_mapped_fragments: 8
  pair_rule: any         # any | all | disabled
  minimum_proper_pairs: 5
  minimum_proper_pair_percent: 5
  minimum_mapq: 0
  # Defaults to assembly.insert_size_min/max for PE validation:
  insert_size_min: 150
  insert_size_max: 450
  maximum_zero_run: 75
  maximum_adjacent_window_fold: 15
  maximum_softclip_percent: 60
  kmer_size: 27
  minimum_unique_kmer_fraction: 0.10
  require_unique_kmers: false
  max_alignments: 10
```

With `pair_rule: any`, either the absolute pair threshold **or** the percentage threshold is sufficient. `all` requires both.

## Compare profiles without assembly

The `validate` command measures one candidate set against one common read sample and then applies several threshold profiles. MIRA, target selection, and candidate extraction are not rerun:

```bash
target-assembler validate paired-balanced.yaml \
  --candidates Deg_balanced/Repeats_candidates.fasta \
  --reads1 cleaned_R1.fastq.gz \
  --reads2 cleaned_R2.fastq.gz \
  --profiles lenient balanced strict \
  --max-units 1000000 \
  --outdir Deg_validation_comparison
```

If `--reads1/--reads2` are omitted, the command uses `input.clean_reads1/2` or preprocesses `input.reads1/2`. It writes `validation_comparison.tsv`, profile-specific `variants_<profile>.tsv`, and profile-specific verified/rejected FASTA files.

## Temporary space and WSL

All sampled libraries, indexes, MIRA assemblies, and other heavy working files are created under `project.temp_dir`. At the end of a run they are deleted by default; this also happens after each sample in batch mode. When `temp_dir` is `null`, Python uses the operating system default, normally `/tmp` on Linux/WSL. Free space on a Windows drive is not necessarily the free space visible to `/tmp` inside WSL. Check both the capacity and inode count with:

```bash
df -h /tmp /mnt/d
df -i /tmp /mnt/d
```

Place temporary data on a roomy mounted drive when needed:

```yaml
project:
  outdir: Deg_balanced
  temp_dir: /mnt/d/Spivrobitnyki/Tynkevich/Science/Articles/Val-po-planu/2026/Heracleum_Final/Target_test/target_tmp

runtime:
  temp_space_check: true
  temp_space_safety_factor: 1.20
  temp_space_margin_gb: 1.0
```

The directory is created automatically. `sampling.maximum_units` remains the scientific/speed cap; for example, `2000000` limits a very low-copy library even when the adaptive estimate requests nearly all reads.

To retain only the reusable FASTQ intermediates (`clean`, `sample`, `target`, and `expanded`), opt in with YAML or a command-line flag:

```yaml
outputs:
  save_intermediate_libraries: true
```

```bash
target-assembler run analysis.yaml --save-intermediate-libraries
```

They are written to `intermediate_libraries/` inside the sample result directory.

A later run can feed those files straight back in through `input.clean_reads1/2`
and `input.sample_reads1/2` to skip preprocessing and sampling. **Supplying
`clean_reads1/2` alone is also valid**, and is what you want when re-running with a
different HMM: preprocessing is independent of the profile and expensive, while
adaptive sampling is a function of the profile and cheap, so the cleaned library
should be reused and the sample rebuilt. A sample sized against one profile will
under-recruit for a profile covering a different fraction of the target. The save step then
targets the file it has just read; it detects that the library is already in place
and leaves it alone, recording `saved_<prefix>_library_already_in_place` instead of
failing with `SameFileError`. For complete debugging state, including MIRA work directories and indexes, use `runtime.keep_work: true` or `--keep-work`. Both options can require substantial disk space and are disabled by default.

## Single-end libraries

For SE reads, omit `reads2` and set `layout: single`:

```yaml
input:
  reads1: reads.fastq.gz
  reads2: null
  layout: single
  hmm: target.hmm
```

The MIRA manifest contains one data file and omits template size and segment placement. Proper-pair thresholds never reject an SE candidate. Because SE evidence is less informative about long-range continuity, the default `auto` profile chooses `lenient`, but every other threshold remains editable.

## MIRA technology and assembly controls

`assembly.technology` is an ordinary YAML value rather than a hard-coded CLI choice:

```yaml
assembly:
  technology: solexa
  rounds: 2
  mira_max_units: 30000
  insert_size_min: 150
  insert_size_max: 450
  segment_placement: "---> <---"
  job: denovo,est,accurate
  parameters:
    - COMMON_SETTINGS
    - -NW:cmrnl=no
    - -NW:ctp=no
    - -NW:cac=no
```

The value is passed to MIRA after a conservative character-safety check, allowing other MIRA-supported technologies without changing Python code. For Illumina data, use `solexa`.

`assembly.mira_threads` (default `null` = `runtime.threads`) sets MIRA's own thread count, passed as `-GE:not` after `parameters` - so a `-GE:not` written into `parameters` has no effect. On several threads MIRA can return a different contig set from identical input; `mira_threads: 1` makes a run repeatable while fastp, nhmmer, BLAST and bowtie2 keep `runtime.threads`.

### MIRA settings for rDNA targets

MIRA 5 has a **built-in rRNA read filter** (`-CL:frrna`). It was on with the
upstream `job: denovo,est,accurate` settings (MIRA's own parameter dump: `Filter
rRNA (frrna) : yes`), and for an rDNA target it removes the target itself. On one ash
library (43,458 round-2 reads) MIRA put ~97 % of its input into the debris list,
30,599 reads of them as `CLIP_RRNA` / `CLIP_RRNA_PAIR`; no contig crossed the 5.8S
core. With the filter off and nothing else changed, contigs carried the whole
ITS1-5.8S-ITS2 for the first time. Two practical points:

- the option belongs under `COMMON_SETTINGS`; under `SOLEXA_SETTINGS` MIRA 5 aborts
  ("Aborting process probably due to an error in the input data or parameters");
- MIRA's lossless digital normalisation (`-KS:ldn`) then becomes the largest
  read reducer (90 % of the same input). It was left on in these runs, which then
  returned full-length ITS in every sample; whether it costs minor ribotypes was
  not tested separately.

```yaml
assembly:
  parameters:
    - COMMON_SETTINGS
    - -NW:cmrnl=no
    - -NW:ctp=no
    - -NW:cac=no
    - -CL:frrna=no       # rDNA targets: do not let MIRA filter rRNA reads
  mira_threads: 1        # fork commit 5: repeatable contigs
```

MIRA prints the values it actually used at the top of its stdout log
(`Filter rRNA (frrna) : no`); check them there rather than in the manifest.

`mira_max_units` controls only the target-enriched library entering each MIRA round. It does not change the genome-proportion sample size. A value of `0` disables MIRA downsampling.

## Genome proportion

The estimator uses the following procedure:

1. sample cleaned reads/pairs;
2. align the accepted variants and build a profile HMM;
3. search sampled read ends with `nhmmer`;
4. merge overlapping hit intervals per read;
5. divide repeat-associated nucleotides by all sampled nucleotides.

By default, two estimates are written to `genome_proportion.tsv`:

- `aligned_nt_percent` uses exact alignment coordinates (`ali_from/ali_to`) and is the default reported metric;
- `envelope_nt_percent` uses the HMM envelope coordinates (`env_from/env_to`) as an optional alternative estimate.

Set `proportion.include_envelope_metric: false` to omit the envelope columns when `aligned_nt` is the reported metric.

R1 and R2 receive unique `/1` and `/2` identifiers before interval merging, preventing the two ends from being treated as one sequence. `replicates: 3` or more adds a sampling SD and approximate 95% interval in `genome_proportion_summary.json`.

This is an estimate of the total genomic representation recognized by the profile HMM; it is not an abundance estimate for each individual candidate variant.

For a proportion-only run, provide an accepted FASTA explicitly (or place `Repeats_verified.fasta` in the configured output directory):

```yaml
project:
  mode: proportion
input:
  clean_reads1: cleaned_R1.fastq.gz
  clean_reads2: cleaned_R2.fastq.gz
  layout: paired
  variants_fasta: accepted_variants.fasta
```

`input.hmm` is not required in this mode because a new profile is built from `variants_fasta`.

## Main outputs

| File | Meaning |
|---|---|
| `Repeats_extracted_raw.fasta` | all repeat instances extracted from final MIRA contigs |
| `Repeats_candidates.fasta` | exact dereplication representatives with stable IDs |
| `Repeats_verified.fasta` | candidates passing the active validation profile |
| `Repeats_rejected.fasta` | candidates retained for inspection but failing one or more criteria |
| `variants.tsv` | provenance, support metrics, decision, and explicit failure reasons |
| `validation_comparison.tsv` | same measured evidence with side-by-side decisions from standalone `validate` |
| `genome_proportion.tsv` | per-replicate aligned and envelope estimates |
| `genome_proportion_summary.json` | selected estimate, variability, and interval |
| `effective_config.yaml` | fully resolved settings used for the run |
| `pipeline_summary.json` | machine-readable status and key counts |
| `pipeline.log` | commands and processing events |

## Testing

Pure Python logic is tested without requiring the external bioinformatics programs:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python -m compileall -q src
```

`target-assembler check analysis.yaml` additionally verifies actual inputs and external dependencies before a biological run.

## Reproducibility

- All random sampling uses `runtime.seed` and independent deterministic offsets per stage/replicate.
- The final resolved YAML is saved with the results.
- Every scientific failure raises a non-zero process exit code.
- Each run uses a newly created MIRA work directory, preventing stale contigs from being accepted.