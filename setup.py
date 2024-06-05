from setuptools import setup, find_packages

setup(
    name='panscan',
    version='0.1',
    packages=find_packages(),
    scripts=['scripts/script1.pl', 'scripts/script2.pl'],
    entry_points={
        'console_scripts': [
            'scanpy=scanpy.__main__:main',
        ],
    },
    install_requires=[
        "pandas"
    ],
)