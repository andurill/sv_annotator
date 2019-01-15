
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
        t_format = None
        if sv.__chr1 == "X":
            t_format = sv.__chr1, sv.__chr2, cband1, cband2, sv.__chr1, sv.__pos1, sv.__chr2, sv.__pos2
        elif sv.__chr2 == "X":
            t_format = sv.__chr2, sv.__chr1, cband2, cband1, sv.__chr2, sv.__pos2, sv.__chr1, sv.__pos1  
        elif int(chr1) < int(chr2):
            t_format = sv.__chr1, sv.__chr2, cband1, cband2, sv.__chr1, sv.__pos1, sv.__chr2, sv.__pos2
        else:
            t_format = sv.__chr2, sv.__chr1, cband2, cband1, sv.__chr2, sv.__pos2, sv.__chr1, sv.__pos1
        coordinate = "t(%s,%s)(%s;%s)(chr%s:g.%s::chr%s:g.%s)".format(t_format)
    else:
        if sv.__isFusion or ():
            cdna1 = get_cdna_pos(self.__bkp1, self.__transcript1)
            cdna2 = get_cdna_pos(self.__bkp2, self.__transcript2)
        elif sv.__isGene1InPanel and sv.__isGene2InPanel:
            cdna1 = get_cdna_pos(self.__bkp1, self.__transcript1)
            cdna2 = get_cdna_pos(self.__bkp2, self.__transcript2)

        elif sv.__isGene1InPanel:
            cdna1 = get_cdna_pos(self.__bkp1, self.__transcript1)
            if not cdna1.startswith("chr"):
                cdna1 = cdna1 + ":" + self.__gene1
            cdna2 = "chr" + self.__bkp2.replace(":", ":g.")
        else:
            cdna2 = get_cdna_pos(self.__bkp2, self.__transcript1)
            cdna1 = "chr" + self.__bkp1.replace(":", ":g.")
        if not cdna1.startswith("chr"):
            cdna1 = cdna1 + ":" + self.__gene1

        if sv.__isFusion:
            if self.__fusionPartner1 == self.__gene1




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
    #rx = re.compile("\([A|C|G|T]\)")
    #actual_ref = rx.findall(request.text)
    s = request.text
    actual_ref = s[s.find("(")+1:s.find(")")] 
    if actual_ref not in ["A", "C", "G", "T"] or actual_ref == dummy_ref:
      request.raise_for_status()
    query = make_query(bkp, actual_ref)
    request = make_get_request(query)
    if not request.ok:
      request.raise_for_status()

  decoded = request.json()
  result = dict(str(s).split(':', 1) for s in decoded[0])
  for tx in transcripts:
    cdna = result[tx]
    if cdna:
      cdna = cdna[:-3]
      if cdna.startswith("-"):
          cdna = "chr" + bkp.replace(":", ":g.")
      return tx, cdna
  raise cdnaNotFound()


def make_get_request(query):
  server = "http://grch37.rest.ensembl.org"
  ext = "/variant_recoder/human/" + query
  try:
    request = requests.get(
        server+ext, headers={"Content-Type": "application/json"})
    except requests.exceptions.RequestException as e:
        raise e
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


def reformat(svtype):
    if svtype == "DUPLICATION":
        return "dup"
    elif svtype == "DELETION":
        return "del"
    else:
        return "inv"