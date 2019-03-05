import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cvr-sv-annotator",
    version="1.0.0",
    author="Gowtham Jayakumaran, Anoop Balakrishnan Rema",
    author_email="jayakumg@mskcc.org, balakra1@mskcc.org",
    description="sv-annotator is a python package for generating human-readable interpretations of structural variants called by iCallSV (https://github.com/rhshah/iCallSV). It uses the rules and logic that are currently in place in clinbx manual review and annotation.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/andurill/sv_annotator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7.3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
