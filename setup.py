import os
import sys
import subprocess
import shutil

try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup, find_packages


with open("requirements.txt") as f:
    required = f.read().splitlines()

# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)

setup(
    name='virchap',
    version='1.0',
    description='Viral haplotype reconstruction using long reads',
    url='https://github.com/gaoyun-25/virCHap',
    author='Yun Gao',
    packages=find_packages(),
    install_requires=required,
    entry_points={
        'console_scripts': [
            'virchap = src.phase_pipeline:main',
            'virchap_pipeline = src.virchap_pipeline:main',
        ],
    },
    python_requires=">=3.8",
)

