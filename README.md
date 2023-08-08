QCEngine
========
[![Build Status](https://github.com/MolSSI/QCEngine/workflows/CI/badge.svg?branch=master)](https://github.com/MolSSI/QCEngine/actions?query=workflow%3ACI)
[![codecov](https://img.shields.io/codecov/c/github/MolSSI/QCEngine.svg?logo=Codecov&logoColor=white)](https://codecov.io/gh/MolSSI/QCEngine)
[![Documentation Status](https://img.shields.io/github/actions/workflow/status/MolSSI/QCEngine/CI.yml?label=docs&logo=readthedocs&logoColor=white)](https://molssi.github.io/QCEngine/)
[![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/qcengine?color=blue&logo=anaconda&logoColor=white)](https://anaconda.org/conda-forge/qcengine)
[![Chat on Slack](https://img.shields.io/badge/chat-on_slack-808493.svg?longCache=true&style=flat&logo=slack)](https://join.slack.com/t/qcarchive/shared_invite/enQtNDIzNTQ2OTExODk0LTE3MWI0YzBjNzVhNzczNDM0ZTA5MmQ1ODcxYTc0YTA1ZDQ2MTk1NDhlMjhjMmQ0YWYwOGMzYzJkZTM2NDlmOGM)
![python](https://img.shields.io/badge/python-3.7+-blue.svg)

<!--[![Azure Build Status](https://dev.azure.com/MolSSI/QCArchive/_apis/build/status/MolSSI.QCEngine?branchName=master)](https://dev.azure.com/MolSSI/QCArchive/_build/latest?definitionId=5&branchName=master)-->
Quantum chemistry program executor and IO standardizer ([QCSchema](https://github.com/MolSSI/QCSchema)) for quantum chemistry.

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

>>> inp = qcel.models.AtomicInput(
    molecule=mol,
    driver="energy",
    model={"method": "SCF", "basis": "sto-3g"},
    keywords={"scf_type": "df"}
    )
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
