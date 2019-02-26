# sv-annotator

sv-annotator is a python package for generating human-readable interpretations of structural variants called by iCallSV (https://github.com/rhshah/iCallSV). It uses the rules and logic that are currently in place in clinbx manual review and annotation. 

## Getting Started
Disclaimer: The annotated data that is provided with the current release is limited to IMPACT cv6 panel and refseq canonical transcripts. The scope of this dataset can be expanded or restricted as needed.

To get a copy of the project up and running on your local machine:

git clone https://github.com/andurill/sv_annotator.git

| Tool | Version |
| --- | --- |
| python | >= 2.7.3 |
| pandas | ==0.19.0 |
| numpy | ==1.8.1 |
| requests | >=2.5.0 |
| re | >=2.2.1 |

### Additional dependencies
API:
Ensembl REST API Endpoints - http://grch37.rest.ensembl.org/?content-type=text/html
OncoKB APIs - https://oncokb.org/api/v1/swagger-ui.html

Datasets:
ftp://ftp.ncbi.nih.gov/genomes/MapView/Homo_sapiens/objects/BUILD.37.3/initial_release/ideogram_9606_GCF_000001305.13_850_V1

https://github.com/genome-nexus/genome-nexus-importer/blob/master/data/ensembl_biomart_pfam_grch37.p13.txt


Canonical transcript list (provided in this repo for all refseq coding genes) - refFlat.canonical_all_coding_exons_aa.interval_list

## Tests

```bash
python svannotate.py "DELETION,16:72798880,16:72829520,ZFHX3 / ZFHX3,IGR: 18Kb before ZFHX3(-),Exon 9 of ZFHX3(-),-"

```
Check the test directory for test data and expected results.
