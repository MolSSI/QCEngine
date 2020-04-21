# QCEngine Pre-built Conda Environments

base.yaml
nwchem.yaml
openmm.yaml
psi-nightly.yaml
psi.yaml
rdkit.yaml
README.md
torchani.yaml

The QCEngine program has few requirements on its own `meta.yaml` file, however,
you may want to emulate the server side of things on your own. To help make that
possible, we have provided the various YAML files here which can be used
to quickly and mostly automatically build a working environment for to emulate
the server.

These use the `conda env create` commands (examples below) instead of the
more common `conda create` (the commands are slightly different as of writing,
circa Conda 4.3), so note the difference in commands.

* `base.yaml` is environment specification for general use.
* `psi-nightly.yaml` additionally installs the latest development Psi4 and anything it needs.
* `nwchem.yaml` doesn't install NWChem (no conda package), but it installs QCEngine and
  NWChem dependencies.

Similarly, other environments are tailored for other CMS harnesses. These
are in active use in CI and so should be reliable.

## Requirements to use Environments

1. `git`
2. `conda`
3. `conda` installed `pip` (pretty much always available unless you are in
   some custom Python-less Conda environment such as an `R`-based env.)
4. Network access

## Setup/Install

Run the following command to configure a new environment with the replacements:

* `{name}`: Replace with whatever you want to call the new env
* `{file}`: Replace with target file

```bash
conda env create -n {name} -f {file}
```

To access the new environment:
```bash
conda activate {name}
```

