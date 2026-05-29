# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Working guidelines (Karpathy)

Behavioral guidelines (after [Andrej Karpathy](https://x.com/karpathy/status/2015883857489522876)) to reduce common LLM coding mistakes; they apply to all work in this repo. The same content is also installed as the on-demand skill `.claude/skills/karpathy-guidelines/`.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

### 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

### 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

### 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

### 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

## What this is

`minisv` is a lightweight mosaic/somatic structural-variation (SV) caller for long
genomic reads (PacBio HiFi / ONT). Its defining idea: align reads against **multiple**
reference genomes / pangenome graphs and keep an SV on a read only if it appears in
the alignments against *all* references. This filters alignment errors and germline/
population SVs without a matched normal.

This is a Python port of [minisv.js](https://github.com/lh3/minisv) plus extra features.
**Output parity with the JS implementation is a design goal** — many functions carry
`NOTE:` comments documenting JS quirks that are deliberately preserved. The shell tests
in `tests/` diff Python output against `minisv.js` run via `k8`.

## Commands

```sh
make                        # == poetry install (the only build step)
poetry run minisv --help    # CLI entry point is minisv.cli:cli (rich-click)
poetry run pytest tests/    # Python unit tests (only test_constant.py is meaningful;
                            # test_cigar_ds.py is mostly unimplemented stubs)
pre-commit run --all-files  # black (py3.10) + isort(--profile black) + ruff lint/format
```

There is **no lint/test target in the Makefile** and no CI; quality is enforced only by
the pre-commit hooks above.

The `tests/*.sh` scripts are integration checks, **not runnable here**: they hardcode
lab-cluster paths (`/hlilab/...`, `/homes6/...`) and require the `k8`/`minisv.js`
reference binary and large alignment files. Read them to understand expected CLI usage
and the diff-against-JS validation pattern, not to execute.

### Cython (optional, currently disabled)

`build.py` cythonizes `minisv/*.pyx` (only `cyminisv.pyx`) with `-march=native -O3`, but
the `[tool.poetry.build]` hook is commented out in `pyproject.toml`, so **the package
installs as pure Python**. The Cython path (`cy_GafParser`) is experimental and unused by
the active CLI.

## The core pipeline

The standalone caller is three streaming steps (see README for the germline / somatic /
mosaic / trio variants — they differ only in inputs and which calls are kept):

```
extract  →  isec  →  (sort -k1,1 -k2,2n)  →  merge  →  genvcf
```

1. **`extract`** (`read_parser.load_reads` → `breakpoint.get_breakpoint` +
   `indel.get_indel`): streams one read group at a time from PAF/GAF, emits raw SVs in a
   BED-like minisv format (`.rsv`/`.gsv`).
2. **`isec`** (`cli.py` + `type.get_type`): keeps only SVs on a read that appear in
   *every* input file (the multi-reference filter).
3. **`merge`** (`merge.merge_sv`): clusters overlapping records into calls, counts
   supporting reads per sample/strand, applies count/strand/centromere-distance filters.
   **Requires input pre-sorted** with `sort -k1,1 -k2,2n`.
4. **`genvcf`** (`io.write_vcf`): minisv `.msv` → VCF. minisv never infers genotypes.

### Input requirements (easy to get wrong)

Alignments **must** carry the `ds:Z` tag — minigraph ≥0.21 (`-cxlr --ds`) or minimap2
≥2.28 (`-cx map-hifi -s50 --ds` for HiFi, `-cxlr:hq` for ONT). In `load_reads`, records
without a `cg:Z` tag or that are not primary (`tp:A:P`) are silently skipped.

### File formats / conventions

- `.rsv` / `.gsv`: raw per-read SVs. Two record shapes distinguished by column 3:
  **breakpoint** records have an orientation token (`>>`,`<<`,`><`,`<>`) and use
  `col_info=8` (graph) or `6`; **indel** records use `col_info=4`. This `col_info`
  offset logic recurs across `isec`, `merge.parse_sv`, and `type.get_type`.
- `.msv`: merged calls (`merge` output).
- Info field is a `;`-delimited string: `SVTYPE=,SVLEN=,qoff_l=,qoff_r=,source=,count=,reads=,...`.
- `source=` is the sample name (set via `extract -n`); somatic calling = grep the wanted
  sample out of merged output (e.g. `grep TUMOR | grep -v NORMAL`).
- SV type encoding (`type.py`): INS/DUP→flag 1, DEL→2, INV→4, cross-contig BND→8.

## Module map

`cli.py` is the hub: every subcommand lives here and builds one of the option dataclasses
(`opt` for extract, `mergeopt`, `unionopt`, `EvalOpt`, `viewopt`) that thread through the
pipeline.

| Module | Role |
|---|---|
| `read_parser.py` | `load_reads` — group reads, dispatch to indel + breakpoint extractors |
| `breakpoint.py` | split-read breakends; `infer_svtype` classifies DEL/INS/DUP/INV/BND from orientation + gaps; `get_end_coor` maps path coords |
| `indel.py` | cigar/`ds:Z` parsing for long INDELs; TSD/polyA (retrotransposon) detection |
| `merge.py` | `merge_sv`/`same_sv`/`write_sv` — cluster records, count support, filter |
| `type.py` | SV-type flags for `isec` overlap tests |
| `eval.py` | `eval` command (callset comparison); shared VCF/minisv parsing (`gc_parse_sv`), interval tree (`iit_*`), bed reader (`gc_read_bed`) |
| `io.py` | `write_vcf` (genvcf), `gc_cmd_view` (view), `parseNum` ("500k"→int) |
| `union.py` | `union_sv`/`advunion_sv`/`union_sv_with_tr` — ensemble callsets into binary-membership truth sets |
| `ensemble.py` | `insilico_truth` (ensembleunion / collapse), `double_strand_break` |
| `filtercaller.py` | heavyweight filtering orchestration — see below |
| `annot.py` | `annot` command |
| `phase.py` | `extracthp` / `annotatehp` haplotype-tag handling |
| `annotation.py`, `graph_genome_coor.py`, `regex.py` | centromere distance/overlap, GAF-path→contig coords, shared compiled regexes |
| `minisv.py` | legacy `GafParser` behind `getindel`/`getsv` (older interface) |
| `identify_breaks_v5.py`, `merge_break_pts_v3.py`, `paired_merge_v2.py` | superseded versioned implementations — not on the active code path |

### `filtercaller.py` — the realtime cross-reference filter

This drives `sv-cross-ref-filter` (somatic), `sv-trio-filter` (de novo in a trio), and
`test-sv-filter` (cutoff sweep on an existing workdir) via three near-parallel classes:
`MinisvReads`, `MinisvReadsTrio`, `MSVTestCutOff`. They share a method sequence:

```
extract_read_ids (from caller VCF + read-id TSV)
  → extract_reads (samtools/seqtk from BAM/CRAM)
  → align_reads_to_self / _to_graph (mappy + external minimap2/minigraph via --mm2/--mg)
  → parse_raw_sv_* (reuse the extract pipeline)
  → isec_* → othercaller_filterasm → union_filtered_vcf
```

The goal is to re-derive SV evidence for a third-party caller's calls (Severus, SAVANA,
nanomonsv, Sniffles2, etc.) by realigning only the supporting reads to a de novo assembly
and pangenome, then keep only calls still supported. The three classes are largely
copy-pasted with per-mode tweaks — **changes to one usually need mirroring in the others.**
External tool paths come from `--mm2`/`--mg`; intermediate `.gsv.gz`/`.paf.gz`/`.gaf.gz`
and per-caller `*_filtered.{vcf,stat}` files land in the output workdir.

## Packaged data

`data/*.bed` and `minisv/*.bed` are centromere masks (`*.cen-mask.bed`) and confident
regions (`*.reg.bed`) for hs38 and chm13v2 — passed via `-b`/`--maskb` to suppress
spurious calls in centromeric/satellite repeats. `minisv/paired.chm13.*.vcf` is a packaged
breakpoint reference.

## Planned changes (in progress)

Three improvements agreed with the user. Full plan file:
`~/.claude/plans/now-let-s-plan-first-partitioned-orbit.md`. Design decisions were made via
explicit user choices — do not re-litigate them without asking.

**Status (synced 2026-05-27):** 3 of 7 code pieces landed in the working tree (uncommitted).
**Task 1 complete** (function + CLI command + unit tests in `tests/test_gnomad_filter.py`),
**Task 2 design locked, ready to execute next session** (sub-features 2a + 2b below, all
decisions made — do not re-litigate), Task 3 partial. Docs/README not yet updated.

### 1. gnomAD-SV population filter — new standalone `gnomadfilter` command
Drop a caller's calls that fall in common population SV regions. **Decision: BED-based**
(type-agnostic breakpoint overlap, no allele-frequency parsing); **standalone command**
(not wired into the orchestrators), caller-agnostic because VCF parsing is shared.
- **[DONE]** New `gnomad_filter(vcf_file, gnomad_bed, opt, both_ends=False, pad=0, out=None)`
  in `minisv/filtercaller.py:173`: loads BED with `gc_read_bed` (`eval.py:881`), parses the
  caller VCF **permissively** with `gc_parse_sv(vcf_file, 0, 0, opt.ignore_flt, False)`
  (`eval.py:184`, which sets `svid=t[2]` and fills `ctg/pos/ctg2/pos2`), builds a **drop-set**
  of `parse_svid(t.svid)` for SVs whose breakpoint overlaps the BED via `iit_overlap`
  (`eval.py:773`, guarded by `t.ctg in bed`) — default drop if **either** end overlaps;
  `both_ends` requires both — then re-emits the original VCF minus the drop-set.
- **[DONE]** New CLI command `gnomadfilter` in `cli.py` (after `filterasm`; args `gnomad_bed`,
  `vcffile`; options `-F/--ignoreflt`, `--both`, `--pad`; filtered VCF to stdout). Builds an
  `EvalOpt(ignore_flt=...)` and calls `gnomad_filter`.

### 2. Cross-caller consensus export — `MinisvReads.output_consensus`  **[TODO — execute next session]**
Export SVs called by **all of severus/savana/nanomonsv** from the **raw** VCFs, regardless
of asm/pangenome filtering, **and** the subset of that consensus the de-novo-assembly filter
dropped. **Decision: first-3 callers, strict** (`self.som_vcfs[:3]`, present in every one).
Both products are emitted by the **single `--output_consensus` flag**. Runs **inside
`MinisvReads` only** (not trio / cutoff classes), and **after `union_filtered_vcf`** because
it consumes that step's `l+s_union.msv`.

**2a. Consensus set (raw):**
- Add keyword `min_file_count=None` to `union_sv` in `minisv/union.py`; in the per-group
  loop (after the `in_bed`/length/count skips ~`union.py:97-99`) add
  `if min_file_count is not None and bin(x).count("1") < min_file_count: continue`. Default
  `None` keeps all existing callers unchanged.
- New `MinisvReads.output_consensus(self, read_min_len, opt)` (next to `union_filtered_vcf`
  ~`filtercaller.py:947`): `union_sv(self.som_vcfs[:3], ..., min_file_count=3)` →
  `consensus_union.msv`, then `insilico_truth` (`ensemble.py:6`) →
  `consensus_3caller_dedup.msv`.

**2b. Consensus-lost-by-assembly subset (the requested feature):** emit consensus calls the
de-novo assembly filtered out → `consensus_lost_by_l+s.msv`. **Locked decisions:**
  - **Scope = `l+s` only** (reads re-aligned to the self de-novo assembly; *not* `l+g` /
    `l+g+s`). "Lost" is measured against `l+s_union.msv` (union of the 3 callers' l+s-filtered
    survivors, already produced by `union_filtered_vcf` ~`filtercaller.py:963`).
  - **Drop rule = lost everywhere**: a consensus call qualifies **iff it overlaps no call in
    `l+s_union.msv`** — all 3 callers agreed in raw, but after re-aligning supporting reads to
    the assembly no caller's filtered set still carries it.
  - **Output = `.msv`**, folded into `--output_consensus` (no separate flag).
  - Set difference via a **small local helper** that reuses the existing
    `gc_cmp_same_sv1(opt.win_size, opt.min_len_ratio, …)` coordinate test and the same
    sorted-scan/`win_size` window already used in `union_sv` — invent no new matching logic.
    Test consensus reps (`consensus_3caller_dedup.msv`) against `l+s_union.msv`; emit reps with
    no overlap.

**Wiring:** add `--output_consensus` flag to the `sv_cross_ref_filter` command; after the
`union_filtered_vcf(...)` call (~`cli.py:807`) do `if output_consensus: reads.output_consensus(...)`.

**Tests** (style of `tests/test_gnomad_filter.py`, synthetic `.msv` inputs): (i) `min_file_count=3`
keeps only all-three-caller groups, `None` unchanged; (ii) subtraction helper excludes a
consensus call with an overlapping l+s survivor, emits one with none. No JS-parity check (new
functionality, no `minisv.js` counterpart).

**Pre-existing inconsistency to leave alone** (flagged, do not "fix" silently): filtered VCFs are
*written* with `opt.min_count` in the name (`filtercaller.py:932`) but *read* with
`opt.read_min_count` (`filtercaller.py:960`). This design consumes the `l+s_union.msv` product,
so it is insulated; do not touch unless asked.

### 3. CLI hygiene
- **[DONE]** Replaced `__version__ = "0.1.2"` (`cli.py:23`) with
  `importlib.metadata.version("minisv")` + a `PackageNotFoundError` fallback to `"0.1.3"` so
  it tracks `pyproject.toml` and can't drift again.
- **[DONE]** Deleted the duplicate `getsv` command (commit `68dc9a1`). The commented-out
  `#def getsv(` stub and the unrelated `# getsv options` comments in `cli.py` were left.

The Module-map note at the top has been updated (the obsolete "version disagrees / getsv is a
duplicate" sentence is removed). `gnomadfilter` + `--gnomadaf` are now in the README
(commit `2e9e2d8`); adding `--output_consensus` to the README is still TODO.

## Planned: Rust port of the `extract` step (`load_reads` → `.rsv`/`.gsv`)  **[TODO — design locked, not started]**

Port the **entire `extract` pipeline step** — `read_parser.load_reads` + `indel.get_indel` +
`breakpoint.get_breakpoint` — to Rust, exposed as a **PyO3 extension** that the `extract` CLI
command calls in place of the Python `load_reads`. This is the successor to the abandoned Cython
path (`cyminisv.pyx` / `cy_GafParser`, wired in `build.py`). **Decisions locked via AskUserQuestion
(2026-05-27) — do not re-litigate without asking:**

- **Scope = the whole extract step.** Rust reads gzipped PAF/GAF, does tokenize +
  `mapq < min_mapq` / primary-only (`tp:A:P`) / `cg:Z`-required filtering + group-by-qname
  (`read_parser.py:42-108`), runs the indel + breakend SV inference, and emits the same
  `.rsv`/`.gsv` records. SAM input stays out of scope (already stubbed/commented in `load_reads`).
- **Integration = PyO3 extension** (maturin/PyO3 native module, e.g. `minisv_rs`), imported into
  Python; `cli.py`'s `extract` calls `minisv_rs.extract(...)`. Adds a cargo+maturin build step — a
  departure from the current pure-Python install; the commented-out `[tool.poetry.build]` Cython
  hook in `pyproject.toml` is the template for wiring an optional native build.

**HARD CONSTRAINT — output parity.** Rust output must be **byte-identical** to the current Python
`load_reads`→`get_indel`/`get_breakpoint`, which is itself diffed against `minisv.js`. The Rust port
must faithfully reproduce every JS-quirk `NOTE:` behavior in `indel.py` / `breakpoint.py` (incl.
the `col_info` 8-vs-6 graph-vs-linear offset, `infer_svtype` orientation/gap classification,
`get_end_coor` path→contig coord mapping, and TSD/polyA retrotransposon detection). Keep Python as
the **reference oracle**: add a locally-runnable test that runs Rust and Python on the same
PAF/GAF and diffs the `.gsv` (extends the `tests/*.sh` diff-against-JS pattern, but without needing
the `k8`/`minisv.js` binary).

**Suggested sequencing (de-risk parity first):**
1. Rust crate + maturin/PyO3 skeleton; no-op `extract` round-tripping a PAF/GAF line. Verify:
   `import minisv_rs` works after `maturin develop`.
2. Port `load_reads` parse/group/filter; diff read grouping against Python on a sample PAF/GAF.
3. Port `get_indel` (cigar/`ds:Z` parse, TSD/polyA) → diff `.gsv` indel records byte-for-byte.
4. Port `get_breakpoint` (`infer_svtype`, `get_end_coor`) → diff `.gsv` breakpoint records.
5. Wire `cli.py` `extract` to call Rust when `minisv_rs` is importable, else fall back to Python.
   Verify: `extract` output unchanged on test inputs.

**Recommendation:** make the extension **optional with a Python fallback** — preserves pure-Python
installability and keeps Python as a permanent parity oracle. Decide the build backend (maturin vs
poetry + separate maturin build) at step 1.

**Bonus:** `filtercaller.py`'s realtime filter reuses this same extract pipeline via its
`parse_raw_sv_*` methods, so they inherit the speedup for free once `extract` is ported.

<!-- CODEGRAPH_START -->
## CodeGraph

This project has a CodeGraph MCP server (`codegraph_*` tools) configured. CodeGraph is a tree-sitter-parsed knowledge graph of every symbol, edge, and file. Reads are sub-millisecond and return structural information grep cannot.

### When to prefer codegraph over native search

Use codegraph for **structural** questions — what calls what, what would break, where is X defined, what is X's signature. Use native grep/read only for **literal text** queries (string contents, comments, log messages) or after you already have a specific file open.

| Question | Tool |
|---|---|
| "Where is X defined?" / "Find symbol named X" | `codegraph_search` |
| "What calls function Y?" | `codegraph_callers` |
| "What does Y call?" | `codegraph_callees` |
| "How does X reach/become Y? / trace the flow from X to Y" | `codegraph_trace` (one call = the whole path, incl. callback/React/JSX dynamic hops) |
| "What would break if I changed Z?" | `codegraph_impact` |
| "Show me Y's signature / source / docstring" | `codegraph_node` |
| "Give me focused context for a task/area" | `codegraph_context` |
| "See several related symbols' source at once" | `codegraph_explore` |
| "What files exist under path/" | `codegraph_files` |
| "Is the index healthy?" | `codegraph_status` |

### Rules of thumb

- **Answer directly — don't delegate exploration.** For "how does X work" / architecture questions, answer with 2-3 codegraph calls: `codegraph_context` first, then ONE `codegraph_explore` for the source of the symbols it surfaces. For a specific **flow** ("how does X reach Y") start with `codegraph_trace` from→to — one call returns the whole path with dynamic hops bridged — then ONE `codegraph_explore` for the bodies; don't rebuild the path with `codegraph_search` + `codegraph_callers`. Codegraph IS the pre-built index, so spawning a separate file-reading sub-task/agent — or running a grep + read loop — repeats work codegraph already did and costs more for the same answer.
- **Trust codegraph results.** They come from a full AST parse. Do NOT re-verify them with grep — that's slower, less accurate, and wastes context.
- **Don't grep first** when looking up a symbol by name. `codegraph_search` is faster and returns kind + location + signature in one call.
- **Don't chain `codegraph_search` + `codegraph_node`** when you just want context — `codegraph_context` is one call.
- **Don't loop `codegraph_node` over many symbols** — one `codegraph_explore` call returns several symbols' source grouped in a single capped call, while each separate node/Read call re-reads the whole context and costs far more.
- **Index lag — check the staleness banner, don't guess a wait.** When a codegraph response starts with "⚠️ Some files referenced below were edited since the last index sync…", the listed files are pending re-index — Read those specific files for accurate content. Files NOT in that banner are fresh and codegraph is authoritative for them. `codegraph_status` also lists pending files under "Pending sync".

### If `.codegraph/` doesn't exist

The MCP server returns "not initialized." Ask the user: *"I notice this project doesn't have CodeGraph initialized. Want me to run `codegraph init -i` to build the index?"*
<!-- CODEGRAPH_END -->
