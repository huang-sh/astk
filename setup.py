from setuptools import setup
from setuptools import find_packages


with open("README.md", "r") as fh:
    description = fh.read()


setup(
    name="astk",
    version="0.1.4",
    url='https://github.com/huang-sh/astk/',
    author='huangsh',
    author_email='hsh-me@outlook.com',
    packages=find_packages(), 
    description="A comprehensive toolkit for alternative splicing analysis",
    long_description=description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    package_data = {
    '': ['data/motif/*/*.meme', "data/motif/ELM/*", "data/maxent/*", "cli/aliases.ini"],
    },
    install_requires=[
        "click>=8.0.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "scikit-learn>=0.24.0",
        "statsmodels>=0.12.0",
        "gseapy==0.10.8",  ## Temporarily specified version, the latest version has bugs when install 
        "statannotations",
        "seaborn",
        "tqdm",
        "pysam",
        "pyfaidx",
        "nease",
        "biopython",
        "pyecharts",
        "click_option_group"
    ],
    entry_points={
        "console_scripts": [
            "astk=astk.cli:cli_fun"
        ]
    },
)
