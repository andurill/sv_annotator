import os
import pandas
import sys
import definitions as de
import models as m

def main():
    pass

def translation_annotation(chrom, coord):
    error = None
    which_cytoband = cb_df[cb_df['Chr'].values == chrom & \
    cb_df['Bp_start'].values <= coord & \
    cd_df['Bp_stop'].values >= coord]
    if len(which_cytoband) == 0:
        raise m.MissingCytoBandError as error
    elif len(which_cytoband) > 1:
        raise m.MultipleCytobandError
    which_cytoband = refFlat_df[(refFlat_df['Gene'].values == gene) & \
    (refFlat_df['cDNA_start'].values <= int(cDNA)) & \
    (refFlat_df['cDNA_end'].values >= int(cDNA))].values.item(exon_index)
    return exon



def nontrans_annotation(value):
    error = "Classified"
    cell_type = pd.start(values)
    translation_annotation(cell_type, coord)


if __name__ == "__main__":
    main()
