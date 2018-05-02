import os
import pandas
import sys
from models import SV


def transcription_annotation(chrom, coord):
    which_cytoband = cb_df[cb_df['Chr'].values == chrom & \
    cb_df['Bp_start'].values <= coord & \
    cd_df['Bp_stop'].values >= coord]
    if len(which_cytoband) == 0:
        raise "Cannot ide"
    which_cytoband = refFlat_df[(refFlat_df['Gene'].values == gene) & \
    (refFlat_df['cDNA_start'].values <= int(cDNA)) & \
    (refFlat_df['cDNA_end'].values >= int(cDNA))].values.item(exon_index)
    return exon




if __name__ == "__main__":
    main()
