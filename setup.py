from setuptools import setup, find_packages

setup(
    name=go_enrichment,
    version=0.1.0,
    packages=find_packages(),
    install_requires=[
        biopython,
        pandas,
        scipy,
        statsmodels,
        requests,
        goatools
    ],
    entry_points={
        console_scripts: [
            go-enrichment=go_enrichment.enrichment:main,
        ],
    },
)

