import os
import requests, sys
import pandas as pd
import numpy as np

refflat=sys.argv[1]
annotate=sys.argv[2]

def main():
  refFlat_df = gen_data_frame(refflat)
  with open(annotate, 'r') as rf, open("annotated.txt", 'w') as wf:
    for line in rf:
      values = line.strip().split()
      positions = annotate_gDNA_positions(values[3], values[10], values[11])
      exons = get_exons_for_segment(values[6], values[10], values[11], refFlat_df)
      new_values = [values[0], positions[0], positions[1]] +\
      [values[i] for i in [1,2,3,4,5,6,7]] +\
      [exons[0], exons[1]] + [values[i] for i in [8,9,10,11]]
      wf.write("\t".join(new_values) + "\n" )
  return


def get_exons_for_segment(gene, cDNA1, cDNA2, refFlat_df):
  exon1 = get_exon(gene, cDNA1, refFlat_df)
  exon2 = get_exon(gene, cDNA2, refFlat_df)
  return (exon1,exon2)


def get_exon(gene, cDNA, refFlat_df):
  exon_index = 6
  exon = refFlat_df[(refFlat_df['Gene'].values == gene) & \
  (refFlat_df['cDNA_start'].values <= int(cDNA)) & \
  (refFlat_df['cDNA_end'].values >= int(cDNA))].values.item(exon_index)
  return exon


def gen_data_frame(file):
  refFlat_df = pd.read_csv(file, delimiter="\t", header=None, \
        names=["Chr", "Start", "End", "Strand", "Gene", "Transcript",\
        "Exon", "cDNA", "AA"], dtype='U')
  refFlat_df[['cDNA_start','cDNA_end']] = refFlat_df['cDNA'].str.split('-',expand=True).astype(int)
  refFlat_df[['AA_start', 'AA_end']] = refFlat_df['AA'].str.split('-', expand=True).astype(int)
  return refFlat_df


def annotate_gDNA_positions(transcript, cDNA1, cDNA2):
  start = get_position_from_vep(transcript, cDNA1)
  stop = get_position_from_vep(transcript, cDNA2)
  return (start, stop)


def get_position_from_vep(transcript, cDNA):
  server = "http://grch37.rest.ensembl.org"
  ext = "/vep/human/hgvs/" + transcript + ":c." + cDNA + "del?"
  #ENST00000003084:c.1431_1433delTTC?"
  r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
  
  if not r.ok:
    r.raise_for_status()
    sys.exit()
 
  decoded = r.json()
  result = decoded[0]
  return str(result["start"])


if __name__ == "__main__":
  main()