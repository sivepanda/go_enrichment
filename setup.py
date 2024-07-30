from setuptools import setup, find_packages

setup(
    name='go_enrichment',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'scipy',
        # Add other dependencies as needed
    ],
    entry_points={
        'console_scripts': [
            'go_enrichment=go_enrichment.enrichment:main',
        ],
    },
)

