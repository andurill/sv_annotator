import requests, sys

if len(sys.argv) < 2:
  print 'Please provide a variant'
  sys.exit(1)

server = "http://grch37.rest.ensembl.org"
ext = "/variant_recoder/human/AGT:c.803T>C?"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print repr(decoded)
