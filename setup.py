from setuptools import setup
from setuptools import find_packages

setup(
    name="astk",
    version="0.1.1",
    url='https://github.com/huang-sh/astk/',
    author='huangsh',
    author_email='hsh-me@outlook.com',
    packages =  find_packages(),
    package_dir = {"astk": "astk"},
    include_package_data=True,
    install_requires=[
        "click>=8.0.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "scikit-learn>=0.24.0",
        "statsmodels>= 0.12.0"
    ],
    entry_points="""
        [console_scripts]
        astk=astk.astk:cli
    """,
)