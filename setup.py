from setuptools import setup, find_packages

setup(
    name='panscan',
    version='0.3.2',
    packages=find_packages(),
    scripts=['scripts/preprocessVCF.pl', 'scripts/findUniqVariants.pl', 'scripts/findNovelSeq.pl'],
    entry_points={
        'console_scripts': [
            'scanpy=scanpy.__main__:main',
        ],
    },
    install_requires=[
        "pandas"
    ],
)
