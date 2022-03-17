from setuptools import setup
from setuptools import find_packages

setup(
    name="astk",
    version="0.1.2",
    url='https://github.com/huang-sh/astk/',
    author='huangsh',
    author_email='hsh-me@outlook.com',
    packages=find_packages(),    
    include_package_data=True,
    package_data = {
    '': ['data/motif/*/*.meme', "data/motif/ELM/*"],
    },
    install_requires=[
        "click>=8.0.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "scikit-learn>=0.24.0",
        "statsmodels>= 0.12.0",
        "pysam",
        "pybedtools",
        "nease"
    ],
    entry_points={
        "console_scripts": [
            "astk=astk.cli:cli_fun"
        ]
    },
)
