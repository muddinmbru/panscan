from setuptools import setup, find_packages

setup(
    name='panscan',
    version='0.2',
    packages=find_packages(),
    scripts=['scripts/findNovelSeq.pl', 'scripts/preprocessVCF.pl'],
    entry_points={
        'console_scripts': [
            'panscan=panscan.__main__:main',
        ],
    },
    install_requires=[
        "pandas"
    ],
)