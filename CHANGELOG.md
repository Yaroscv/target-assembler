# Changelog - branch `species-hunter` of Yaroscv/target-assembler

Upstream: [yurij-tynkevich/target-assembler](https://github.com/yurij-tynkevich/target-assembler) v0.1.0 (MIT).
Every key below defaults to `null`, which is the upstream behaviour.

## 0.1.0+specieshunter.5 (2026-10-05)

Version string in `pyproject.toml` and `src/target_assembler/__init__.py`; a PEP 440
local version label on upstream 0.1.0, so it never collides with an upstream release.
`CITATION.cff` is unchanged: it cites the upstream software release.

1. **`selection.hmm_score_min` / `selection.hmm_evalue`** - nhmmer `-T` / `-E` at
   every call site, re-applied after the table is parsed; event
   `hmm_hits_below_threshold`. `-T` wins when both are set.
2. **No self-copy of a resumed library** - `_save_library` leaves a library that is
   already the file being saved in place; event
   `saved_<prefix>_library_already_in_place` (was: `SameFileError`).
3. **`input.clean_reads1/2` without `input.sample_reads1/2`** - accepted; the
   sample is rebuilt for the current profile.
4. **`expansion.max_internal_clip`** - round 2 rejects local alignments clipped by
   more than this many bp where the contig continues; clipping at contig ends is
   never counted; event `expansion_internal_clip_rejected`.
5. **`assembly.mira_threads`** - MIRA's thread count, appended as `-GE:not` after
   `assembly.parameters`; `1` gives repeatable contigs.

Documentation (README): the five keys; MIRA's rRNA filter (`-CL:frrna=no` under
`COMMON_SETTINGS`) for rDNA targets; why tandem extraction needs two repeat
periods on one contig, and `site` mode for a whole-35S profile.