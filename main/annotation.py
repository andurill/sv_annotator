
def annotate_variant(sv):
    sv.__bkpOrder = {}
    sv.__translocationCoordOrder = {}
    if self.__isKnownFusion:
        call = "fusion"
    else:
        call = "rearrangement"
    annotation = get_annotation(sv)
    coordinates = get_coordinate(sv)


def get_coordinates(sv):
    if sv.__svtype == "TRANSLOCATION":
        cband1 = get_cytoband(sv.__bkp1)
        cband2 = get_cytoband(sv.__bkp2)
        if sv.__chr1 == "X":
            coordinate = "t(%s,%s)(%s;%s)(chr%s:g.%s::chr%s:g.%s)".\
            format(sv.__chr1, sv.__chr2, cband1, cband2, sv.__chr1, sv.__pos1, sv.__chr2, sv.__pos2)
        elif sv.__chr2 == "X":
            coordinate = "t(%s,%s)(%s;%s)(chr%s:g.%s::chr%s:g.%s)".\
            format(sv.__chr2, sv.__chr1, cband2, cband1, sv.__chr2, sv.__pos2, sv.__chr1, sv.__pos1)  
        elif int(chr1) < int(chr2):
            coordinate = "t(%s,%s)(%s;%s)(chr%s:g.%s::chr%s:g.%s)".\
            format(sv.__chr1, sv.__chr2, cband1, cband2, sv.__chr1, sv.__pos1, sv.__chr2, sv.__pos2)
        else:
            coordinate = "t(%s,%s)(%s;%s)(chr%s:g.%s::chr%s:g.%s)".\
            format(sv.__chr2, sv.__chr1, cband2, cband1, sv.__chr2, sv.__pos2, sv.__chr1, sv.__pos1)
    else:
        if sv.__isFusion or (sv.__isGene1InPanel and sv.__isGene2InPanel):
            cdna1 = get_cdna_pos(self.__bkp1, self.__transcript1)
            cdna2 = get_cdna_pos(self.__bkp2, self.__transcript2)
        elif sv.__isGene1InPanel:
            cdna1 = get_cdna_pos(self.__bkp1, self.__transcript1)
            cdna2 = "chr" + self.__bkp2.replace(":", ":g.")
            coordinate = cdna1 + ":" + self.__gene1 + "_" + cdna2 + refomat(sv.__svtype)
        elif sv.__isGene2InPanel:
        if sv.__isFusion:
            if self.__fusionPartner1 == self.__gene1:




def get_cytoband(bkp):
    chrom, coord = bkp.split(":")
    which_cytoband = cb_df[cb_df['Chr'].values == chrom &
                           cb_df['Bp_start'].values <= coord &
                           cd_df['Bp_stop'].values >= coord]
    if len(which_cytoband) == 0:
        raise MissingCytoBand(bkp)
    elif len(which_cytoband) > 1:
        raise MultipleCytoband(bkp)
    return which_cytoband


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
  raise cdnaNotFound()


def make_get_request(query):
  server = "http://grch37.rest.ensembl.org"
  ext = "/variant_recoder/human/" + query
  request = requests.get(
      server+ext, headers={"Content-Type": "application/json"})
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
