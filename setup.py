"""
Setup configuration for FAST TCP package
"""

from setuptools import setup, find_packages
import os


# Read the README file
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()


# Read requirements
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        return [
            line.strip() for line in fh if line.strip() and not line.startswith("#")
        ]


setup(
    name="fast-tcp",
    version="1.0.0",
    author="Breno Miranda, Emilio Cruciani, Roberto Verdecchia, Antonia Bertolino",
    author_email="",
    description="FAST Approaches to Scalable Similarity-based Test Case Prioritization",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/icse18-FAST/FAST",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Software Development :: Testing",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=read_requirements(),
    entry_points={
        "console_scripts": [
            "fast-tcp=fast_tcp.cli:main",
        ],
    },
    keywords="test case prioritization, software testing, similarity-based, LSH, FAST",
    project_urls={
        "Bug Reports": "https://github.com/icse18-FAST/FAST/issues",
        "Source": "https://github.com/icse18-FAST/FAST",
        "Documentation": "https://github.com/icse18-FAST/FAST/blob/master/README.md",
    },
    include_package_data=True,
    zip_safe=False,
)
