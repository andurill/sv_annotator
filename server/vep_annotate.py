import os, requests, sys, re


def get_cdna_pos(bkp, transcripts):
  dummy_ref = "C"
  query = make_query(bkp, dummy_ref)

  # get request max twice to VEP for annotation
  request = make_get_request(query)
  if not request.ok:
    rx = re.compile("\([A|C|G|T]\)")
    actual_ref = rx.findall(request.text)
    if not actual_ref:
      request.raise_for_status()
    dummy_ref = str(actual_ref[0]).replace("(", "").replace(")", "")
    query = make_query(bkp, dummy_ref)
    request = make_get_request(query)
    if not request.ok:
      request.raise_for_status()
  
  decoded = request.json()
  result = dict(str(s).split(':', 1) for s in decoded[0])
  for tx in transcripts:
    cdna = result[tx]
    if cdna:
      cdna = cdna[:-3]
      return tx, cdna
  raise m.cdnaNotFound()


def make_get_request(query):
  server = "http://grch37.rest.ensembl.org"
  ext = "/variant_recoder/human/" + query
  try:
    request = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
  except ConnectionError
  return request


def make_query(bkp, dummy_ref):
  chrom, pos = bkp.split(":")
  revcomp = {
      "A": "T",
      "T": "A",
      "C": "G",
      "G": "C",
      "N": "N"
  }
  dummy_alt = revcomp[dummy_ref]
  query = chrom + ":g." + pos + dummy_ref + ">" + dummy_alt + "?"
  return query