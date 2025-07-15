from setuptools import setup, find_packages

setup(
    name='pdbe_sifts',
    version='0.1.0',
    description='Python package to run SIFTS from EMBL-EBI locally.',
    author='EMBL-EBI, Adam Bellaiche, Preeti Choudhary, Sreenath Sasidharan Nair, Jennifer Fleming, Sameer Velankar',
    author_email='adamb@ebi.ac.uk',
    url='https://gitlab.ebi.ac.uk/',
    packages=find_packages(include=['pdbe_sifts', 'pdbe_sifts.*']),
    install_requires=[
        'tqdm',
        'biopython',
        'funcy',
        'lxml',
        'requests',
        'pandas',
        'pymmseqs',
        'coloredlogs',
        'pyyaml',
        'scikit-learn',
        'ete4',
    ],
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
