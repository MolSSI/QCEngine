import setuptools
import versioneer

if __name__ == "__main__":
    setuptools.setup(
        name='dqm_compute',
        description='Compute module for the DQM project',
        author='Daniel G. A. Smith',
        author_email='dgasmith@vt.edu',
        url="https://github.com/psi4/dqm_compute",
        license='BSD-3C',
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass(),
        packages=setuptools.find_packages(),
        install_requires=[
            'pyyaml',
            'py-cpuinfo',
        ],
        extras_require={
            'docs': [
                'sphinx==1.2.3',  # autodoc was broken in 1.3.1
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
            ],
        },

        tests_require=[
            'pytest',
            'pytest-cov',
        ],

        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=True,
        long_description="""
"""
    )
