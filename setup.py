from setuptools import setup, find_packages

setup(
    name='dcc_sra',
    version='0.0.1',
    description=("Transfer sequence files and metadata from the iHMP DCC "
                 "to NCBI's SRA"),
    packages=["dcc_sra"],
    zip_safe=False,
    install_requires=[
        'anadama',
        'osdf-python',
        'cutlass',
        'paramiko',
        'python-dateutil>=2.2'
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha"
    ],
    entry_points= {
        'anadama.pipeline': [
            ".dcc_sra = dcc_sra.pipeline:DCCSRAPipeline"
        ]
    }
)
