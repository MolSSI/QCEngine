# CCLib Q-Chem and ORCA Fixture Regressions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add version-grouped, inline Q-Chem and ORCA input/output regression tables and one live HF smoke calculation per harness.

**Architecture:** Each dedicated module owns an explicit fixture inventory as parametrization rows. Input rows compare the generator's complete native text without native execution. Output rows replay cclib fixture output through the existing conversion path and compare saved dictionaries excluding `extras`; live availability-marked rows alone invoke native software.

**Tech Stack:** Python, pytest, QCEngine, QCElemental, cclib fixture corpus.

## Global Constraints

- Fixture rows are grouped and identified by the native version encoded in their cclib data directory: Q-Chem 5.1/5.4/6.0 and ORCA 5.0/6.0/6.1.
- Input and output expected values are inline `pytest.mark.parametrize` arguments.
- Input/output parametrized tests must never execute Q-Chem or ORCA.
- Output comparisons remove only `extras`.
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


### Task 2: Collect inline expectations

**Files:**
- Create: `qcengine/programs/cclib_programs/tests/test_qchem.py`
- Create: `qcengine/programs/cclib_programs/tests/test_orca.py`

**Interfaces:**
- Consumes: fixture path, `AtomicInput`, `TaskConfig`, and program `build_input`/`_parse_and_convert` APIs.
- Produces: `(version, input_model, expected_input)` and `(version, relative_output, expected_output)` rows.

- [ ] **Step 1: Write a failing exact-input row**

Start each module with one minimal HF/STO-3G He row:
```python
@pytest.mark.parametrize("version,input_model,expected_input", [("minimum", HE_INPUT, "")])
def test_generated_inputs(version, input_model, expected_input):
    assert build_input(input_model, TASK_CONFIG, "/resolved/program").input_text == expected_input
```

- [ ] **Step 2: Verify the test fails because expected text is empty**

Run:
```bash
pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_generated_inputs -q
pytest qcengine/programs/cclib_programs/tests/test_orca.py::test_generated_inputs -q
```
Expected: each fails with generated native text differing from `""`.

- [ ] **Step 3: Collect exact native input and parsed dictionaries**

Use a local one-shot Python collector to create `AtomicInput` from the cclib writer result, call the program generator with fixed `TaskConfig(ncores=1, memory=1.0)`, replay the complete `.out` through `ExecutionResult` and `_parse_and_convert`, then serialize `result.dict(exclude={"extras"})` using `json.dumps(..., sort_keys=True)`. Preserve the input text and dictionary as Python literals in their corresponding version group.

- [ ] **Step 4: Verify the exact input tests pass without executables**

Run:
```bash
pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_generated_inputs qcengine/programs/cclib_programs/tests/test_orca.py::test_generated_inputs -q
```
Expected: PASS without `qchem` or `orca` invocation.

### Task 3: Implement Q-Chem regressions

**Files:**
- Modify: `qcengine/programs/cclib_programs/tests/test_qchem.py`

**Interfaces:**
- Consumes: Q-Chem 5.1, 5.4, and 6.0 fixture paths and saved expectation rows.
- Produces: `test_generated_inputs`, `test_parsed_outputs`, and `test_live_hf_single_point`.

- [ ] **Step 1: Add a failing output assertion for `basicQChem5.1/water_mp2.out`**

```python
@pytest.mark.parametrize("version,relative_output,expected_output", [
    ("5.1", "QChem/basicQChem5.1/water_mp2.out", {}),
])
def test_parsed_outputs(version, relative_output, expected_output):
    assert normalize(convert_fixture(relative_output)) == expected_output
```

- [ ] **Step 2: Run it with fixture data and verify it fails on `{}`**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_qchem.py::test_parsed_outputs -q
```
Expected: FAIL because the converted result contains fields absent from `{}`.

- [ ] **Step 3: Add all collected Q-Chem version groups and live smoke test**

Populate `5.1`, `5.4`, and `6.0` rows using collector output. Add an availability-guarded `cclib-qchem` HF/STO-3G He (or H2) energy calculation asserting `success` and scalar `return_result`.

- [ ] **Step 4: Run Q-Chem regression tests**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_qchem.py -q
```
Expected: fixture rows PASS; live row PASS or SKIP based on program availability.

- [ ] **Step 5: Commit Q-Chem coverage**

```bash
git add qcengine/programs/cclib_programs/tests/test_qchem.py
git commit -m "test: add qchem cclib fixture regressions"
```

### Task 4: Implement ORCA regressions

**Files:**
- Modify: `qcengine/programs/cclib_programs/tests/test_orca.py`

**Interfaces:**
- Consumes: ORCA 5.0, 6.0, and 6.1 fixture paths and saved expectation rows.
- Produces: `test_generated_inputs`, `test_parsed_outputs`, and `test_live_hf_single_point`.

- [ ] **Step 1: Add a failing output assertion for `basicORCA6.0/dvb_sp_hf.out`**

```python
@pytest.mark.parametrize("version,relative_output,expected_output", [
    ("6.0", "ORCA/basicORCA6.0/dvb_sp_hf.out", {}),
])
def test_parsed_outputs(version, relative_output, expected_output):
    assert normalize(convert_fixture(relative_output)) == expected_output
```

- [ ] **Step 2: Verify it fails on the empty expected dictionary**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_orca.py::test_parsed_outputs -q
```
Expected: FAIL because the converted result contains fields absent from `{}`.

- [ ] **Step 3: Add all collected ORCA version groups and live smoke test**

Populate `5.0`, `6.0`, and `6.1` rows. Include generator-supported energy requests in input rows; record other-driver input exclusions in the inventory. Add an availability-guarded `cclib-orca` HF/STO-3G He energy calculation asserting `success` and scalar `return_result`.

- [ ] **Step 4: Run ORCA regression tests**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_orca.py -q
```
Expected: fixture rows PASS; live row PASS or SKIP based on program availability.

- [ ] **Step 5: Commit ORCA coverage**

```bash
git add qcengine/programs/cclib_programs/tests/test_orca.py
git commit -m "test: add orca cclib fixture regressions"
```

### Task 5: Verify suite boundaries and complete tracking

**Files:**
- Modify: `docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md`
- Test: `qcengine/programs/cclib_programs/tests/test_qchem.py`
- Test: `qcengine/programs/cclib_programs/tests/test_orca.py`
- Test: `qcengine/programs/cclib_programs/tests/test_cclib.py`

- [ ] **Step 1: Mark each completed or excluded inventory row**

For every candidate, retain its version, source/output path, final status, and concrete exclusion reason if omitted.

- [ ] **Step 2: Run the targeted suite with fixture data**

Run:
```bash
CCLIB_SOURCE_ROOT=/home/awallace43/gits/cclib pytest qcengine/programs/cclib_programs/tests/test_qchem.py qcengine/programs/cclib_programs/tests/test_orca.py qcengine/programs/cclib_programs/tests/test_cclib.py -q
```
Expected: PASS, with only availability-guarded live tests skipped when native software is unavailable.

- [ ] **Step 3: Commit tracking and verification updates**

```bash
git add -f docs/superpowers/plans/2026-07-31-cclib-qchem-orca-fixture-regressions.md
git commit -m "docs: complete cclib fixture regression plan"
```
