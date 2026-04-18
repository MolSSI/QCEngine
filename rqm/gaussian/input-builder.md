# Feature: Gaussian Harness — Input Builder <!-- rq-9cdf3f7d -->

This document specifies the input-file construction layer for the Gaussian harness.
It covers two modules:

- `qcengine/programs/gaussian/germinate.py` — maps QCSchema method/driver/multiplicity
  to Gaussian route-section components.
- `qcengine/programs/gaussian/keywords.py` — assembles the Gaussian `.com` file from
  Link 0 directives, the route section, and the molecule block.

The runner calls these modules from `build_input()` (see `runner.md`).

## Gaussian Input File Format <!-- rq-9ec1917a -->

A Gaussian `.com` file has five sections, separated by blank lines:

```
%NProcShared=<N>
%Mem=<M>GB
#P <method>/<basis> <job_type> [<extra_keywords>]

<title>

<charge> <multiplicity>
<symbol> <x> <y> <z>
...

```

- **Link 0** (`%` lines): processor and memory directives.
- **Route section** (`#` line): method, basis, job type, and additional keywords.
- **Title card**: a single non-blank line (may be `"QCEngine Gaussian calculation"`).
- **Charge/multiplicity**: on one line, space-separated integers.
- **Atom block**: element symbol followed by three Cartesian coordinates in Ångström,
  one atom per line.
- A trailing blank line terminates the molecule section.

## Feature API <!-- rq-2ad053f2 -->

### `muster_modelchem(method: str, driver: str, multiplicity: int) -> Dict[str, Any]` <!-- rq-ef62979b -->

Located in `qcengine/programs/gaussian/germinate.py`.

**Purpose**: converts the QCSchema `method`, `driver`, and molecular `multiplicity` into
the Gaussian route-section keywords needed to describe the calculation.

**Parameters**:
- `method` — method string from `AtomicInput.model.method` (case-insensitive). <!-- rq-a88c706f -->
- `driver` — one of `"energy"`, `"gradient"`, `"hessian"`, `"properties"`. <!-- rq-b49bb2a2 -->
- `multiplicity` — from `AtomicInput.molecule.molecular_multiplicity`. <!-- rq-4d763c6d -->

**Returns** a dict with keys:

| Key | Type | Description |
|-----|------|-------------|
| `"method_string"` | `str` | The Gaussian method token, e.g. `"B3LYP"`, `"UMP2"`. |
| `"job_type"` | `str` | Gaussian job-type keyword appended to the route section. |
| `"extra_keywords"` | `Dict[str, str]` | Additional route-section tokens (e.g. `{"Pop": "Full"}`). |

**Method mapping** (case-insensitive input):

| Input method | Gaussian method token | Notes |
|--------------|-----------------------|-------|
| `"hf"` or `"scf"` | `"HF"` | |
| `"mp2"` | `"MP2"` | |
| `"ccsd"` | `"CCSD"` | |
| `"ccsd(t)"` | `"CCSD(T)"` | |
| anything else | passed through verbatim (uppercased) | Treated as a DFT functional or Gaussian-native method; not validated. |

**Driver → job_type mapping**:

| Driver | `job_type` |
|--------|-----------|
| `"energy"` | `""` (no extra keyword; Gaussian defaults to single-point) |
| `"gradient"` | `"Force=NoStep"` |
| `"hessian"` | `"Freq"` |
| `"properties"` | `""` (Pop=Full added via `extra_keywords`) |

**Open-shell handling** (mirrors the Psi4 pattern):

If `multiplicity > 1` and the method string (after uppercasing) does not already begin
with `"U"` or `"RO"`, the prefix `"U"` is automatically prepended (e.g. `"HF"` →
`"UHF"`, `"B3LYP"` → `"UB3LYP"`). This matches Gaussian's unrestricted-reference
naming convention.

If `multiplicity == 1`, the method string is used as-is (closed-shell reference assumed
unless the user explicitly provides `"U..."` or `"RO..."` in the method name).

**Properties driver extras**:

When `driver == "properties"`, `extra_keywords` includes `{"Pop": "Full"}` so that
Mulliken charges and dipole moment are printed in the log.

**Raises**:

- `InputError` — if `driver` is not one of the four supported values. <!-- rq-86c27356 -->

---

### `build_com_file(link0: Dict[str, Any], route_line: str, title: str, charge: int, multiplicity: int, atom_block: str) -> str` <!-- rq-0052a2c7 -->

Located in `qcengine/programs/gaussian/keywords.py`.

**Purpose**: assembles all sections of a Gaussian `.com` file into a single string.

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `link0` | `dict` | Key/value pairs written as `%Key=Value` lines. Keys are `"NProcShared"` and `"Mem"`. |
| `route_line` | `str` | Full route section string (already formatted, starting with `#P`). |
| `title` | `str` | Title card (single line). |
| `charge` | `int` | Molecular charge. |
| `multiplicity` | `int` | Spin multiplicity. |
| `atom_block` | `str` | Newline-separated lines of `"<symbol> <x> <y> <z>"` in Ångström. |

**Returns** the complete `.com` file as a string with a trailing newline.

**Link 0 directives**:

| Keyword | Source | Format |
|---------|--------|--------|
| `%NProcShared` | `config.ncores` | `%NProcShared=<N>` |
| `%Mem` | `max(1, int(config.memory))` GB | `%Mem=<N>GB` |

---

### `build_route_line(method_string: str, basis: str, job_type: str, extra_keywords: Dict[str, str], user_keywords: Dict[str, Any]) -> str` <!-- rq-3c343d2a -->

Located in `qcengine/programs/gaussian/keywords.py`.

**Purpose**: formats the `#P` route section line.

**Behaviour**:
- Starts with `#P <method_string>/<basis>`.
- Appends `job_type` if non-empty.
- Appends each `(key, value)` pair from `extra_keywords` as `Key=Value` (or bare `Key`
  if value is empty).
- Appends pass-through keywords from `user_keywords`: any key not in the reserved set
  `{"memory", "nprocs"}` is appended verbatim as `Key=Value` (or bare `Key` if value is
  falsy).

**Example output**:

```
#P B3LYP/6-31G* Force=NoStep
#P UMP2/cc-pVDZ Freq
#P CCSD(T)/cc-pVTZ
#P HF/STO-3G Pop=Full
```

## Molecule Block Construction <!-- rq-74e0110c -->

The atom block is obtained from `input_model.molecule.to_string(dtype="xyz", units="angstrom")`,
keeping only the coordinate lines (i.e. stripping the atom-count and comment lines
present in the XYZ format), or equivalently by iterating over `molecule.symbols` and
`molecule.geometry` (converting from Bohr to Ångström).

Coordinates must be written with sufficient decimal places (at least 10) to preserve
numerical precision.

Ghost atoms (`real=False` in the molecule) are not supported; `InputError` is raised if
any atom is a ghost atom.

## Gherkin Scenarios <!-- rq-cc549245 -->

```gherkin
Feature: Gaussian input builder

  Background:
    Given a TaskConfig with ncores=4, memory=8.0, scratch_directory=None, scratch_messy=False

  # --- muster_modelchem: method mapping ---

  @rq-0199bb02
  Scenario: HF method maps to "HF" method string
    When muster_modelchem("hf", "energy", 1) is called
    Then method_string is "HF"
    And job_type is ""
    And extra_keywords is {}

  @rq-7049aaa5
  Scenario: SCF is an alias for HF
    When muster_modelchem("scf", "energy", 1) is called
    Then method_string is "HF"

  @rq-546faafb
  Scenario: MP2 method maps to "MP2"
    When muster_modelchem("mp2", "energy", 1) is called
    Then method_string is "MP2"

  @rq-2052a2fa
  Scenario: CCSD method maps to "CCSD"
    When muster_modelchem("ccsd", "energy", 1) is called
    Then method_string is "CCSD"

  @rq-92492893
  Scenario: CCSD(T) method maps to "CCSD(T)"
    When muster_modelchem("ccsd(t)", "energy", 1) is called
    Then method_string is "CCSD(T)"

  @rq-851d782c
  Scenario: Unknown method is passed through verbatim (uppercased)
    When muster_modelchem("B3LYP", "energy", 1) is called
    Then method_string is "B3LYP"

  @rq-079b1e22
  Scenario: Method matching is case-insensitive
    When muster_modelchem("CCSD(T)", "energy", 1) is called
    Then method_string is "CCSD(T)"

  # --- muster_modelchem: driver mapping ---

  @rq-c266b9bb
  Scenario: Energy driver produces no job-type keyword
    When muster_modelchem("hf", "energy", 1) is called
    Then job_type is ""

  @rq-f1c116db
  Scenario: Gradient driver produces "Force=NoStep"
    When muster_modelchem("hf", "gradient", 1) is called
    Then job_type is "Force=NoStep"

  @rq-82942b90
  Scenario: Hessian driver produces "Freq"
    When muster_modelchem("hf", "hessian", 1) is called
    Then job_type is "Freq"

  @rq-f512b60b
  Scenario: Properties driver produces Pop=Full in extra_keywords
    When muster_modelchem("hf", "properties", 1) is called
    Then job_type is ""
    And extra_keywords contains {"Pop": "Full"}

  @rq-e008c6fe
  Scenario: Unsupported driver raises InputError
    When muster_modelchem("hf", "optimization", 1) is called
    Then InputError is raised

  # --- muster_modelchem: open-shell handling ---

  @rq-972b6470
  Scenario: Singlet uses method name as-is (no U prefix)
    When muster_modelchem("hf", "energy", 1) is called
    Then method_string is "HF"

  @rq-6d73bb9e
  Scenario: Doublet automatically prepends U to HF
    When muster_modelchem("hf", "energy", 2) is called
    Then method_string is "UHF"

  @rq-c946c285
  Scenario: Triplet automatically prepends U to B3LYP
    When muster_modelchem("B3LYP", "energy", 3) is called
    Then method_string is "UB3LYP"

  @rq-25c4ad49
  Scenario: Method already starting with U is not double-prefixed
    When muster_modelchem("UHF", "energy", 2) is called
    Then method_string is "UHF"

  @rq-fe13245b
  Scenario: Method starting with RO is not prefixed with U
    When muster_modelchem("ROHF", "energy", 2) is called
    Then method_string is "ROHF"

  # --- build_route_line ---

  @rq-af032242
  Scenario: Route line for closed-shell HF energy
    When build_route_line("HF", "STO-3G", "", {}, {}) is called
    Then the result is "#P HF/STO-3G"

  @rq-9acd1ce5
  Scenario: Route line for gradient calculation
    When build_route_line("HF", "STO-3G", "Force=NoStep", {}, {}) is called
    Then the result starts with "#P HF/STO-3G"
    And the result contains "Force=NoStep"

  @rq-cd0d6637
  Scenario: Route line for properties calculation includes Pop=Full
    When build_route_line("HF", "STO-3G", "", {"Pop": "Full"}, {}) is called
    Then the result contains "Pop=Full"

  @rq-7152455e
  Scenario: User keyword is appended to route line
    When build_route_line("HF", "STO-3G", "", {}, {"SCF": "Tight"}) is called
    Then the result contains "SCF=Tight"

  @rq-e7897f09
  Scenario: Reserved user keywords (memory, nprocs) are not appended to route line
    When build_route_line("HF", "STO-3G", "", {}, {"memory": "8GB", "nprocs": "4"}) is called
    Then the result does not contain "memory"
    And the result does not contain "nprocs"

  # --- build_com_file ---

  @rq-d7e8161d
  Scenario: Link 0 section contains NProcShared and Mem directives
    When build_com_file({"NProcShared": 4, "Mem": "8GB"}, "#P HF/STO-3G", "title", 0, 1, "H 0.0 0.0 0.0\nH 0.0 0.0 0.74") is called
    Then the result starts with "%NProcShared=4\n%Mem=8GB"

  @rq-342a2ca3
  Scenario: Charge and multiplicity are on one line in the molecule section
    When build_com_file is called with charge=1 and multiplicity=2
    Then the result contains a line "1 2"

  @rq-b819d8e0
  Scenario: Atom coordinates appear in Angstrom with at least 10 decimal places
    When build_com_file is called with atom_block for water
    Then each coordinate line contains at least 10 decimal digits

  @rq-23d8b2c6
  Scenario: The .com file ends with a trailing blank line
    When build_com_file is called
    Then the returned string ends with "\n\n"

  # --- Ghost atoms ---

  @rq-5053204b
  Scenario: Ghost atoms in the molecule raise InputError
    Given an AtomicInput whose molecule contains a ghost atom (real=False)
    When build_input() is called
    Then InputError is raised

  # --- BasisSet object ---

  @rq-0aa99076
  Scenario: BasisSet object for model.basis raises InputError
    Given an AtomicInput whose model.basis is a BasisSet object
    When build_input() is called
    Then InputError is raised

  # --- Memory and cores ---

  @rq-db92093a
  Scenario: Memory directive rounds down to nearest integer GB (minimum 1)
    Given config.memory is 0.5 (GiB)
    When build_com_file is called
    Then the Link 0 section contains "%Mem=1GB"

  @rq-4ad16d95
  Scenario: Memory directive uses integer GB for whole-number memory
    Given config.memory is 8.0 (GiB)
    When build_com_file is called
    Then the Link 0 section contains "%Mem=8GB"
```
