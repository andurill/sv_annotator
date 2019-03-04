# sv-annotator

sv-annotator is a python package for annotating structural variants called by iCallSV (https://github.com/rhshah/iCallSV) for clinical/research reports. It uses the rules and logic that are currently in place in clinbx manual review and annotation. 

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
python svannotate.py "DELETION,21:39863298,21:42869384,ERG / TMPRSS2,Intron of ERG(-):7Kb after exon 1,Intron of TMPRSS2(-):661bp after exon 2,Protein Fusion: in frame {TMPRSS2:ERG}"

# Example output following a successful execution
{'Note': 'The TMPRSS2 - ERG fusion involves TMPRSS2 exons 1 - 2 and ERG exons 2 - 10. The fusion is predicted to be in frame.', 'Position': 'TMPRSS2 exon 2 to ERG exon 2', 'message': None, 'status': 'SUCCESS', 'Annotation': 'TMPRSS2 (NM_001135099) - ERG (NM_182918) fusion: c.126+662:TMPRSS2_c.18+6989:ERGdel'}

python svannotate.py "DELETION,5:1295168,5:10249820,TERT / FAM173B,Promoter of TERT(-):42Kb from tx start,Intron of FAM173B(-):149bp after exon 1,Transcript Fusion {FAM173B:TERT},3to5"

# Example output following a failure
{'Note': None, 'Position': None, 'message': BothBreakpointsNoncoding('Atleast one of the breakpoints need to be in target panel and in coding region.',), 'status': 'FAILED', 'Annotation': None}
```
