from setuptools import setup, find_packages
print(find_packages())

setup(
    name='go_enrichment',
    version='0.1.0',
    packages=find_packages(include=['go_enrichment']),
    install_requires=[
        'numpy<2',
        'tqdm',
        'pandas',
        'mygene',
        'scipy',
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'go_enrichment=go_enrichment.enrichment:main',
        ],
    },
)

