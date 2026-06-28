from setuptools import setup, find_packages

setup(
    name="nhanes-data-automator",
    version="2.05",
    description="Automated NHANES data extraction with independent R validation",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/lxthyy/nhanes-data-automator",
    license="MIT",
    python_requires=">=3.10",
    install_requires=[
        "pandas>=1.5",
        "numpy>=1.24",
        "scipy>=1.10",
        "statsmodels>=0.14",
        "matplotlib>=3.7",
        "openpyxl>=3.1",
        "python-docx>=0.8",
    ],
    classifiers=[
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
)
