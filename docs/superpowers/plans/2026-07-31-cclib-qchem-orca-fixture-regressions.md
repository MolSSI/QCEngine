# CCLib Q-Chem and ORCA Fixture Regressions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add compact, version-grouped Q-Chem and ORCA input/output regressions and one deferred live H2 HF smoke calculation per harness.

**Architecture:** Each dedicated module owns small representative parametrization tables. Input rows generate compact `AtomicInput` requests and assert only the native regular-expression fragments that exercise unique driver/model/basis, resource, coordinate/unit, and program-keyword behavior. Output rows replay fixtures through the existing conversion path, resolve one dotted QCSchema path, and compare a scalar or small array with `numpy.allclose(atol=1e-6)`. They do not retain full native-input snapshots, result JSON, or fixture-derived `AtomicInput` blobs. Live availability checks remain inside live test bodies, so fixture collection does not probe executables.

**Tech Stack:** Python, pytest, QCEngine, QCElemental, cclib fixture corpus.

## Global Constraints

- Fixture rows are grouped and identified by the native version encoded in their cclib data directory: Q-Chem 5.1/5.4/6.0 and ORCA 5.0/6.0/6.1.
- Input rows are `(version, AtomicInput, regex_fragments)` and output rows are `(version, relative_output, dotted_qcschema_path, expected_scalar_or_small_array)`.
- Input/output parametrized tests must never execute Q-Chem or ORCA.
- Output comparisons use `numpy.allclose(atol=1e-6)`; no result dictionaries are retained.
- Output fixtures are loaded from `$CCLIB_SOURCE_ROOT/data`; unavailable fixture roots skip clearly.
- Live HF/STO-3G single-point tests use the existing program availability marks.

---

### Task 1: Build and record the fixture inventory

**Files:**
- Create: `docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md`

**Interfaces:**
- Consumes: `/home/awallace43/gits/cclib/data/{QChem,ORCA}`.
- Produces: explicit per-version inclusion/exclusion inventory used to form test rows.

- [x] **Step 1: Enumerate paired fixture candidates and parser outcomes**

  Ran the prescribed `.out` pairing inventory against the approved version directories. All 91 `.out` candidates are paired; the complete results and `.log` exclusions are recorded below.

- [x] **Step 2: Record collector failures as exclusions**

  Ran the collector through `ccopen`, expected parser identity, parse, `QCSchemaWriter(...).as_dict(validate=False)`, native input generation, and `_parse_and_convert` for every paired `.out`. Exact per-side statuses are recorded below.

- [x] **Step 3: Commit the inventory baseline**

  Committed as `docs: track cclib fixture regressions`.

## Task 1 inventory baseline

**Status: complete.** The inventory is constrained to Q-Chem 5.1/5.4/6.0 and ORCA 5.0/6.0/6.1. The `.out` enumeration found 91 candidates, all paired with their required native input; therefore there are no missing-pair exclusions. The 18 `.log`-only candidates are explicitly excluded below. Collector status was obtained with cclib `1.8.1.post1235+5c3639de` from `/home/awallace43/gits/cclib`, `TaskConfig(ncores=1, nnodes=1, memory=1.0, scratch_directory=None, retries=0, mpiexec_command=None)`, and the current harness definitions.

For each paired `.out`, the collector ran `ccopen`, checked the exact expected parser class, parsed the output, called `QCSchemaWriter(...).as_dict(validate=False)`, built the v2 `AtomicInput`, generated native input, and replayed the output through `_parse_and_convert`. `included-input` and `included-output` are independently recorded so a future row collector can use either side.

### QChem 5.1

Paired `.out` candidates: 21; included input rows: 17; included output rows: 17.

| Source → output | Input status | Output status |
| --- | --- | --- |
| `QChem/basicQChem5.1/C_bigbasis.in` → `QChem/basicQChem5.1/C_bigbasis.out` | included-input | included-output |
| `QChem/basicQChem5.1/MoOCl4_sp.in` → `QChem/basicQChem5.1/MoOCl4_sp.out` | excluded: parser result validation: metadata.success is not true | excluded: parser result validation: metadata.success is not true |
| `QChem/basicQChem5.1/Trp_polar.in` → `QChem/basicQChem5.1/Trp_polar.out` | included-input | included-output |
| `QChem/basicQChem5.1/dvb_bomd.in` → `QChem/basicQChem5.1/dvb_bomd.out` | included-input | included-output |
| `QChem/basicQChem5.1/dvb_dispersion_bp86_d3zero.in` → `QChem/basicQChem5.1/dvb_dispersion_bp86_d3zero.out` | included-input | included-output |
| `QChem/basicQChem5.1/dvb_gopt.in` → `QChem/basicQChem5.1/dvb_gopt.out` | included-input | included-output |
| `QChem/basicQChem5.1/dvb_ir.in` → `QChem/basicQChem5.1/dvb_ir.out` | included-input | included-output |
| `QChem/basicQChem5.1/dvb_raman.in` → `QChem/basicQChem5.1/dvb_raman.out` | included-input | included-output |
| `QChem/basicQChem5.1/dvb_sp.in` → `QChem/basicQChem5.1/dvb_sp.out` | included-input | included-output |
| `QChem/basicQChem5.1/dvb_sp_un.in` → `QChem/basicQChem5.1/dvb_sp_un.out` | included-input | included-output |
| `QChem/basicQChem5.1/dvb_td.in` → `QChem/basicQChem5.1/dvb_td.out` | included-input | included-output |
| `QChem/basicQChem5.1/water_ccd.in` → `QChem/basicQChem5.1/water_ccd.out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCD | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCD |
| `QChem/basicQChem5.1/water_ccsd(t).in` → `QChem/basicQChem5.1/water_ccsd(t).out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCSD(T) | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCSD(T) |
| `QChem/basicQChem5.1/water_ccsd.in` → `QChem/basicQChem5.1/water_ccsd.out` | included-input | included-output |
| `QChem/basicQChem5.1/water_cis.in` → `QChem/basicQChem5.1/water_cis.out` | included-input | included-output |
| `QChem/basicQChem5.1/water_ir.in` → `QChem/basicQChem5.1/water_ir.out` | included-input | included-output |
| `QChem/basicQChem5.1/water_ir_anharm.in` → `QChem/basicQChem5.1/water_ir_anharm.out` | included-input | included-output |
| `QChem/basicQChem5.1/water_mp2.in` → `QChem/basicQChem5.1/water_mp2.out` | included-input | included-output |
| `QChem/basicQChem5.1/water_mp3.in` → `QChem/basicQChem5.1/water_mp3.out` | included-input | included-output |
| `QChem/basicQChem5.1/water_mp4.in` → `QChem/basicQChem5.1/water_mp4.out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method MP4 | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method MP4 |
| `QChem/basicQChem5.1/water_mp4sdq.in` → `QChem/basicQChem5.1/water_mp4sdq.out` | included-input | included-output |

### QChem 5.4

Paired `.out` candidates: 21; included input rows: 18; included output rows: 18.

| Source → output | Input status | Output status |
| --- | --- | --- |
| `QChem/basicQChem5.4/C_bigbasis.in` → `QChem/basicQChem5.4/C_bigbasis.out` | included-input | included-output |
| `QChem/basicQChem5.4/MoOCl4_sp.in` → `QChem/basicQChem5.4/MoOCl4_sp.out` | included-input | included-output |
| `QChem/basicQChem5.4/Trp_polar.in` → `QChem/basicQChem5.4/Trp_polar.out` | included-input | included-output |
| `QChem/basicQChem5.4/dvb_bomd.in` → `QChem/basicQChem5.4/dvb_bomd.out` | included-input | included-output |
| `QChem/basicQChem5.4/dvb_dispersion_bp86_d3zero.in` → `QChem/basicQChem5.4/dvb_dispersion_bp86_d3zero.out` | included-input | included-output |
| `QChem/basicQChem5.4/dvb_gopt.in` → `QChem/basicQChem5.4/dvb_gopt.out` | included-input | included-output |
| `QChem/basicQChem5.4/dvb_ir.in` → `QChem/basicQChem5.4/dvb_ir.out` | included-input | included-output |
| `QChem/basicQChem5.4/dvb_raman.in` → `QChem/basicQChem5.4/dvb_raman.out` | included-input | included-output |
| `QChem/basicQChem5.4/dvb_sp.in` → `QChem/basicQChem5.4/dvb_sp.out` | included-input | included-output |
| `QChem/basicQChem5.4/dvb_sp_un.in` → `QChem/basicQChem5.4/dvb_sp_un.out` | included-input | included-output |
| `QChem/basicQChem5.4/dvb_td.in` → `QChem/basicQChem5.4/dvb_td.out` | included-input | included-output |
| `QChem/basicQChem5.4/water_ccd.in` → `QChem/basicQChem5.4/water_ccd.out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCD | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCD |
| `QChem/basicQChem5.4/water_ccsd(t).in` → `QChem/basicQChem5.4/water_ccsd(t).out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCSD(T) | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCSD(T) |
| `QChem/basicQChem5.4/water_ccsd.in` → `QChem/basicQChem5.4/water_ccsd.out` | included-input | included-output |
| `QChem/basicQChem5.4/water_cis.in` → `QChem/basicQChem5.4/water_cis.out` | included-input | included-output |
| `QChem/basicQChem5.4/water_ir.in` → `QChem/basicQChem5.4/water_ir.out` | included-input | included-output |
| `QChem/basicQChem5.4/water_ir_anharm.in` → `QChem/basicQChem5.4/water_ir_anharm.out` | included-input | included-output |
| `QChem/basicQChem5.4/water_mp2.in` → `QChem/basicQChem5.4/water_mp2.out` | included-input | included-output |
| `QChem/basicQChem5.4/water_mp3.in` → `QChem/basicQChem5.4/water_mp3.out` | included-input | included-output |
| `QChem/basicQChem5.4/water_mp4.in` → `QChem/basicQChem5.4/water_mp4.out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method MP4 | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method MP4 |
| `QChem/basicQChem5.4/water_mp4sdq.in` → `QChem/basicQChem5.4/water_mp4sdq.out` | included-input | included-output |

### QChem 6.0

Paired `.out` candidates: 13; included input rows: 13; included output rows: 13.

| Source → output | Input status | Output status |
| --- | --- | --- |
| `QChem/basicQChem6.0/water_hf_solvent_cosmo.in` → `QChem/basicQChem6.0/water_hf_solvent_cosmo.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_onsager.in` → `QChem/basicQChem6.0/water_hf_solvent_onsager.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_pcm_cosmo.in` → `QChem/basicQChem6.0/water_hf_solvent_pcm_cosmo.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_pcm_cpcm.in` → `QChem/basicQChem6.0/water_hf_solvent_pcm_cpcm.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_pcm_iefpcm.in` → `QChem/basicQChem6.0/water_hf_solvent_pcm_iefpcm.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_pcm_ssvpe.in` → `QChem/basicQChem6.0/water_hf_solvent_pcm_ssvpe.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_sm12_chelpg.in` → `QChem/basicQChem6.0/water_hf_solvent_sm12_chelpg.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_sm12_cm5.in` → `QChem/basicQChem6.0/water_hf_solvent_sm12_cm5.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_sm12_mk.in` → `QChem/basicQChem6.0/water_hf_solvent_sm12_mk.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_sm8.in` → `QChem/basicQChem6.0/water_hf_solvent_sm8.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_smd.in` → `QChem/basicQChem6.0/water_hf_solvent_smd.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_smd_cpcm.in` → `QChem/basicQChem6.0/water_hf_solvent_smd_cpcm.out` | included-input | included-output |
| `QChem/basicQChem6.0/water_hf_solvent_smd_iefpcm.in` → `QChem/basicQChem6.0/water_hf_solvent_smd_iefpcm.out` | included-input | included-output |

### ORCA 5.0

Paired `.out` candidates: 17; included input rows: 4; included output rows: 4.

| Source → output | Input status | Output status |
| --- | --- | --- |
| `ORCA/basicORCA5.0/Trp_polar.inp` → `ORCA/basicORCA5.0/Trp_polar.out` | included-input | included-output |
| `ORCA/basicORCA5.0/dvb_coupling_nmr.inp` → `ORCA/basicORCA5.0/dvb_coupling_nmr.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_dispersion_bp86_d3zero.inp` → `ORCA/basicORCA5.0/dvb_dispersion_bp86_d3zero.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_gopt.inp` → `ORCA/basicORCA5.0/dvb_gopt.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_ir.inp` → `ORCA/basicORCA5.0/dvb_ir.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_nmr.inp` → `ORCA/basicORCA5.0/dvb_nmr.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_raman.inp` → `ORCA/basicORCA5.0/dvb_raman.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_rocis.inp` → `ORCA/basicORCA5.0/dvb_rocis.out` | included-input | included-output |
| `ORCA/basicORCA5.0/dvb_scan_relaxed.inp` → `ORCA/basicORCA5.0/dvb_scan_relaxed.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_scan_unrelaxed.inp` → `ORCA/basicORCA5.0/dvb_scan_unrelaxed.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_sp.inp` → `ORCA/basicORCA5.0/dvb_sp.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_sp_un.inp` → `ORCA/basicORCA5.0/dvb_sp_un.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/dvb_td.inp` → `ORCA/basicORCA5.0/dvb_td.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA5.0/water_ccsd.inp` → `ORCA/basicORCA5.0/water_ccsd.out` | included-input | included-output |
| `ORCA/basicORCA5.0/water_ccsd_t.inp` → `ORCA/basicORCA5.0/water_ccsd_t.out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCSD(T) | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCSD(T) |
| `ORCA/basicORCA5.0/water_mp2.inp` → `ORCA/basicORCA5.0/water_mp2.out` | included-input | included-output |
| `ORCA/basicORCA5.0/water_mp3.inp` → `ORCA/basicORCA5.0/water_mp3.out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method MP3 | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method MP3 |
| `ORCA/basicORCA5.0/dvb_adc2.inp` → `ORCA/basicORCA5.0/dvb_adc2.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA5.0/dvb_eom_ccsd.inp` → `ORCA/basicORCA5.0/dvb_eom_ccsd.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA5.0/dvb_perf.inp` → `ORCA/basicORCA5.0/dvb_perf.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA5.0/dvb_pno_eom_ccsd.inp` → `ORCA/basicORCA5.0/dvb_pno_eom_ccsd.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA5.0/dvb_steom_ccsd.inp` → `ORCA/basicORCA5.0/dvb_steom_ccsd.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA5.0/dvb_steom_dlpno_ccsd.inp` → `ORCA/basicORCA5.0/dvb_steom_dlpno_ccsd.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA5.0/water_hf_solvent_cpcm.inp` → `ORCA/basicORCA5.0/water_hf_solvent_cpcm.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA5.0/water_hf_solvent_cpcm_cosmo.inp` → `ORCA/basicORCA5.0/water_hf_solvent_cpcm_cosmo.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA5.0/water_hf_solvent_smd.inp` → `ORCA/basicORCA5.0/water_hf_solvent_smd.log` | excluded: log-only candidate | excluded: log-only candidate |

### ORCA 6.0

Paired `.out` candidates: 18; included input rows: 6; included output rows: 6.

| Source → output | Input status | Output status |
| --- | --- | --- |
| `ORCA/basicORCA6.0/Trp_polar.inp` → `ORCA/basicORCA6.0/Trp_polar.out` | included-input | included-output |
| `ORCA/basicORCA6.0/dvb_dispersion_bp86_d3zero.inp` → `ORCA/basicORCA6.0/dvb_dispersion_bp86_d3zero.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/dvb_gopt.inp` → `ORCA/basicORCA6.0/dvb_gopt.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/dvb_ir.inp` → `ORCA/basicORCA6.0/dvb_ir.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/dvb_nmr.inp` → `ORCA/basicORCA6.0/dvb_nmr.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/dvb_raman.inp` → `ORCA/basicORCA6.0/dvb_raman.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/dvb_rocis.inp` → `ORCA/basicORCA6.0/dvb_rocis.out` | included-input | included-output |
| `ORCA/basicORCA6.0/dvb_scan_relaxed.inp` → `ORCA/basicORCA6.0/dvb_scan_relaxed.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/dvb_scan_unrelaxed.inp` → `ORCA/basicORCA6.0/dvb_scan_unrelaxed.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/dvb_sp_dft.inp` → `ORCA/basicORCA6.0/dvb_sp_dft.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/dvb_sp_hf.inp` → `ORCA/basicORCA6.0/dvb_sp_hf.out` | included-input | included-output |
| `ORCA/basicORCA6.0/dvb_sp_un_dft.inp` → `ORCA/basicORCA6.0/dvb_sp_un_dft.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/dvb_sp_un_hf.inp` → `ORCA/basicORCA6.0/dvb_sp_un_hf.out` | included-input | included-output |
| `ORCA/basicORCA6.0/dvb_td.inp` → `ORCA/basicORCA6.0/dvb_td.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |
| `ORCA/basicORCA6.0/water_ccsd.inp` → `ORCA/basicORCA6.0/water_ccsd.out` | included-input | included-output |
| `ORCA/basicORCA6.0/water_ccsd_t.inp` → `ORCA/basicORCA6.0/water_ccsd_t.out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCSD(T) | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method CCSD(T) |
| `ORCA/basicORCA6.0/water_mp2.inp` → `ORCA/basicORCA6.0/water_mp2.out` | included-input | included-output |
| `ORCA/basicORCA6.0/water_mp3.inp` → `ORCA/basicORCA6.0/water_mp3.out` | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method MP3 | excluded: QCSchemaWriter.as_dict(validate=False): RuntimeError: Don't know what to do with method MP3 |
| `ORCA/basicORCA6.0/dvb_adc2.inp` → `ORCA/basicORCA6.0/dvb_adc2.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA6.0/dvb_eom_ccsd.inp` → `ORCA/basicORCA6.0/dvb_eom_ccsd.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA6.0/dvb_perf.inp` → `ORCA/basicORCA6.0/dvb_perf.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA6.0/dvb_pno_eom_ccsd.inp` → `ORCA/basicORCA6.0/dvb_pno_eom_ccsd.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA6.0/dvb_steom_ccsd.inp` → `ORCA/basicORCA6.0/dvb_steom_ccsd.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA6.0/dvb_steom_dlpno_ccsd.inp` → `ORCA/basicORCA6.0/dvb_steom_dlpno_ccsd.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA6.0/water_hf_solvent_cpcm.inp` → `ORCA/basicORCA6.0/water_hf_solvent_cpcm.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA6.0/water_hf_solvent_cpcm_cosmo.inp` → `ORCA/basicORCA6.0/water_hf_solvent_cpcm_cosmo.log` | excluded: log-only candidate | excluded: log-only candidate |
| `ORCA/basicORCA6.0/water_hf_solvent_smd.inp` → `ORCA/basicORCA6.0/water_hf_solvent_smd.log` | excluded: log-only candidate | excluded: log-only candidate |

### ORCA 6.1

Paired `.out` candidates: 1; included input rows: 0; included output rows: 0.

| Source → output | Input status | Output status |
| --- | --- | --- |
| `ORCA/basicORCA6.1/dvb_coupling_nmr.inp` → `ORCA/basicORCA6.1/dvb_coupling_nmr.out` | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' | excluded: QCSchemaWriter.as_dict(validate=False): KeyError: 'functional' |

### Inventory totals

| Program | Paired `.out` | Included inputs | Included outputs | Excluded `.out` candidates | Excluded `.log`-only candidates | Missing-pair candidates |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| QChem | 55 | 48 | 48 | 7 | 0 | 0 |
| ORCA | 36 | 10 | 10 | 26 | 18 | 0 |
| **Total** | **91** | **58** | **58** | **33** | **18** | **0** |


## Compact coverage decision

The tables above are **eligibility inventory**, not a mandate to snapshot every eligible file. The approved compact suite selects the following representative fixture classes. Everything eligibility-marked but absent from this list is intentionally excluded as redundant corpus coverage; unsupported drivers/features and writer failures retain the concrete reasons recorded above.

| Program/version | Included compact output properties | Compact input focus | Explicit compact exclusions |
| --- | --- | --- | --- |
| Q-Chem 5.1 | HF `C_bigbasis`, B3LYP `dvb_sp`, MP2 `water_mp2`, CCSD `water_ccsd`: `properties.return_energy` | energy + scalar keywords/resources | other eligible files are redundant; non-energy fixture classes are represented by generated gradient/hessian inputs |
| Q-Chem 5.4 | HF `C_bigbasis`, B3LYP `dvb_sp`, MP2 `water_mp2`, CCSD `water_ccsd`: `properties.return_energy` | gradient, method/basis, Bohr coordinates | other eligible files are redundant |
| Q-Chem 6.0 | HF solvent `water_hf_solvent_onsager`: `properties.return_energy` | hessian, method/basis/resources | only eligible fixture class is HF solvent; all other solvent variants are redundant |
| ORCA 5.0 | HF `Trp_polar`, MP2 `water_mp2`, CCSD `water_ccsd`: `properties.return_energy` | simple keyword and `%scf` block/resources | eligible duplicate fixtures are redundant; DFT fixture conversion fails with `KeyError: 'functional'` |
| ORCA 6.0 | HF `dvb_sp_hf`, MP2 `water_mp2`, CCSD `water_ccsd`: `properties.return_energy` | simple keyword, `%cpcm` block, Bohr-to-Angstrom coordinates | DFT fixture conversion fails with `KeyError: 'functional'`; other eligible fixtures are redundant |
| ORCA 6.1 | none | MP2 method/basis and `%scf` block/resources | sole output (`dvb_coupling_nmr`) fails conversion with `KeyError: 'functional'` |

### Task 2: Implement compact fixture suites

**Files:**
- Modify: `qcengine/programs/cclib_programs/tests/test_qchem.py`
- Modify: `qcengine/programs/cclib_programs/tests/test_orca.py`

- [x] Replace full native-input snapshots with three small, version-labeled `AtomicInput` rows per suite and regex fragment assertions.
- [x] Cover Q-Chem driver mapping, method/basis, resources, Bohr coordinates, and scalar `$rem` keywords; cover ORCA simple keywords, user blocks, resources, and Bohr-to-Angstrom coordinates.
- [x] Replace full parsed-result dictionaries with the compact output rows documented above. Resolve dotted QCSchema paths and compare with `numpy.allclose(atol=1e-6)`.
- [x] Keep fixture-root skips and defer executable availability probes to live-test bodies.
- [x] Change the ORCA live HF/STO-3G smoke molecule from He to H2.

### Task 3: Verification and tracking

**Files:**
- Modify: `docs/superpowers/specs/2026-07-31-cclib-qchem-orca-fixture-regressions-design.md`
- Modify: this plan
- Create: `.superpowers/sdd/2026-07-31-cclib-qchem-orca-fixture-regressions/compact-rewrite-report.md`

- [ ] Run the dedicated Q-Chem and ORCA suites in the activated `cclib_qcng` environment with `CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib`.
- [ ] Run the cclib test suite with the same fixture root.
- [ ] Record line counts, diff summary, commands, and results in the compact rewrite report.
- [ ] Commit the compact rewrite after validation. Do not modify the cclib checkout.
