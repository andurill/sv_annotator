import requests, sys
 
server = "http://grch37.rest.ensembl.org"
ext = "/vep/human/hgvs/ENST00000003084:c.1431_1433delTTC?"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print repr(decoded)
