# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

setup(
    name="nhanes_downloader",
    version="2.06",
    description="NHANES Data Downloader & Validation Engine — A Standardized Preprocessing Pipeline",
    author="李鑫 (Li Xin)",
    author_email="lxddzyx@126.com",
    license="MIT",
    python_requires=">=3.10",
    install_requires=[
        "pandas>=1.5.0",
        "numpy>=1.24.0",
        "openpyxl>=3.0.0",
        "python-docx>=0.8.11",
        "scikit-learn>=1.2.0",
        "statsmodels>=0.14.0",
        "scipy>=1.10.0",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Intended Audience :: Healthcare Researchers",
    ],
)
