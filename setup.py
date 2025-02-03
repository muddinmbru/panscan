  GNU nano 4.8                                                                                                      setup.py                                                                                                                 
from setuptools import setup, find_packages

setup(
    name="panscan",
    version="0.3.2",
    packages=find_packages(),
    package_data={
        "panscan": [
            "scripts/*",              # Include all files in the scripts directory
            "scripts/perlModules/*"    # Include all files in the perlModules directory
        ]
    },
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "panscan=panscan.__main__:main",  # Assuming you have a `__main__.py` in `panscan/`
        ],
    },
    install_requires=[
        "pandas",
    ],
)
