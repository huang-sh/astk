from setuptools import setup

setup(
    name="astk",
    version="0.1",
    packages = ["astk", "astk.suppa"],
    
    include_package_data=True,
    install_requires=["click", "pandas", "matplotlib", "sklearn", "statsmodels"],
    entry_points="""
        [console_scripts]
        astk=astk.astk:cli
    """,
)