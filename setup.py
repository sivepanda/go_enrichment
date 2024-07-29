from setuptools import setup, find_packages

setup(
    name="go_enrichment",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "biopython",
        "pandas",
        "scipy",
        "statsmodels",
        "requests",
        "goatools"
    ],
    entry_points={
        "console_scripts": [
            "run-enrichment=go_enrichment.go_enrichment:main",
        ],
    },
    python_requires='>=3.6',
)

