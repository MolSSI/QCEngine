QCEngine
========

![Build Status](https://img.shields.io/github/actions/workflow/status/MolSSI/QCEngine/CI.yml?branch=master&logo=github&link=https%3A%2F%2Fgithub.com%2FMolSSI%2FQCEngine%2Factions%3Fquery%3Dworkflow%253ACI)
[![codecov](https://img.shields.io/codecov/c/github/MolSSI/QCEngine.svg?logo=Codecov&logoColor=white)](https://codecov.io/gh/MolSSI/QCEngine)
[![Documentation
Status](https://img.shields.io/github/actions/workflow/status/MolSSI/QCEngine/CI.yaml?label=docs&logo=readthedocs&logoColor=white)](https://molssi.github.io/QCEngine/)
[![Chat on Slack](https://img.shields.io/badge/chat-on_slack-green.svg?longCache=true&style=flat&logo=slack)](https://join.slack.com/t/qcarchive/shared_invite/zt-3calopudd-2rtUC~XN1tj1Zn9MHkV6GQ)
![python](https://img.shields.io/badge/python-3.10+-blue.svg)

**Documentation:** [GitHub Pages](https://molssi.github.io/QCEngine/)

<!--[![Azure Build Status](https://dev.azure.com/MolSSI/QCArchive/_apis/build/status/MolSSI.QCEngine?branchName=master)](https://dev.azure.com/MolSSI/QCArchive/_build/latest?definitionId=5&branchName=master)-->
<!--Quantum chemistry program executor and IO standardizer ([QCSchema](https://github.com/MolSSI/QCSchema)) for quantum chemistry.-->
Quantum chemistry program executor and IO standardizer (QCSchema) for quantum chemistry.

# Example

A simple example of QCEngine's capabilities is as follows:

```python
>>> import qcengine as qcng
>>> import qcelemental as qcel

>>> mol = qcel.models.Molecule.from_data("""
O  0.0  0.000  -0.129
H  0.0 -1.494  1.027
H  0.0  1.494  1.027
""")

>>> inp = qcel.models.AtomicInput(  # QCSchema v1
    molecule=mol,
    driver="energy",
    model={"method": "SCF", "basis": "sto-3g"},
    keywords={"scf_type": "df"}
    )

>>> inp = qcel.models.AtomicInput(  # QCSchema v2
    molecule=mol,
    specification={
        "driver": "energy",
        "model": {"method": "SCF", "basis": "sto-3g"},
        "keywords": {"scf_type": "df"}
    })
```

These input specifications can be executed with the ``compute`` function along with a program specifier:

```python
>>> ret = qcng.compute(inp, "psi4")
```

The results contain a complete record of the computation:


```python
>>> ret.return_result
-74.45994963230625

>>> ret.properties.scf_dipole_moment
[0.0, 0.0, 0.6635967188869244]

>>> ret.provenance.cpu
Intel(R) Core(TM) i7-7820HQ CPU @ 2.90GHz
```

See the [documentation](https://molssi.github.io/QCEngine/) for more information.

# License

BSD-3C. See the [License File](LICENSE) for more information.
