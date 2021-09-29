import sys
import setuptools
import versioneer

short_description = (
    "QCEngine provides a wrapper to ingest and produce QCSchema for a variety of quantum chemistry programs."
)

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except FileNotFoundError:
    long_description = short_description

if __name__ == "__main__":
    setuptools.setup(
        name="qcengine",
        description="A compute wrapper for Quantum Chemistry.",
        author="The QCArchive Development Team",
        author_email="qcarchive@molssi.org",
        url="https://github.com/MolSSI/QCEngine",
        license="BSD-3C",
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        packages=setuptools.find_packages(),
        setup_requires=[] + pytest_runner,
        install_requires=["pyyaml", "py-cpuinfo", "psutil", "qcelemental>=0.23.0", "pydantic>=1.8.2"],
        entry_points={"console_scripts": ["qcengine=qcengine.cli:main"]},
        extras_require={
            "docs": [
                "sphinx==1.2.3",  # autodoc was broken in 1.3.1
                "sphinxcontrib-napoleon",
                "sphinx_rtd_theme",
                "numpydoc",
            ],
            "tests": ["pytest", "pytest-cov"],
            "lint": ["black"],
        },
        tests_require=["pytest", "pytest-cov"],
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3 :: Only",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
        ],
        zip_safe=False,
        long_description=long_description,
        long_description_content_type="text/markdown",
    )
