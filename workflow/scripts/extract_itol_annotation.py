#!/usr/bin/python

# This script is used to extract info from iTOL annotation files which is generated by TB-Profiler
## Including 1. Lineage; 2. Drug resistance category; 3. Individual drug resistance
## Extracted info will be used in R::ggtree to generate phylotree visualization result

import pandas as pd
import sys

# Data path
lineage_file = sys.argv[1]
dr_cat_file = sys.argv[2]
dr_indiv_file = sys.argv[3]

lineage_outfile = sys.argv[4]
dr_cat_outfile = sys.argv[5]
dr_indiv_outfile = sys.argv[6]

# 1. Lineage
with open(lineage_file) as f_lineage:
    buf = f_lineage.readline()
    while buf:
        buf = buf.strip()
        # example: LEGEND_LABELS	Lineage2	Lineage4
        if buf.startswith("LEGEND_LABELS"):
            lineage_list = buf.split('\t')[1:]
        # example: LEGEND_COLORS	#ab2323	#f68e51
        if buf.startswith("LEGEND_COLORS"):
            color_list = buf.split('\t')[1:]
            color_lineage_dict = dict(zip(color_list, lineage_list))
        # example:
        # DATA
        # NBCDC45	#f68e51
        if buf.startswith("DATA") and not buf.startswith("DATASET"):
            lineage_matrix = []
            buf = f_lineage.readline()
            while buf:
                buf = buf.strip()
                line = buf.split('\t')
                # add lineage based on color
                line.append(color_lineage_dict[line[1]])
                lineage_matrix.append(line)
                buf = f_lineage.readline()
        # Read next line
        buf = f_lineage.readline()

lineage_df = pd.DataFrame(lineage_matrix, columns = ["sample", "color", "lineage"])
lineage_df = lineage_df.set_index("sample")
lineage_df.to_csv(lineage_outfile, header=True, index=True, sep='\t')

# 2. Drug resistance category
with open(dr_cat_file) as f_dr_cat:
    buf = f_dr_cat.readline()
    while buf:
        buf = buf.strip()
        # example: LEGEND_LABELS	Sensitive	Pre-MDR	MDR	Pre-XDR	XDR	Other
        if buf.startswith("LEGEND_LABELS"):
            dr_cat_list = buf.split('\t')[1:]
        # example: LEGEND_COLORS	#28a745	#007bff	#ffc107	#dc3545	#343a40	#f8f9fa
        if buf.startswith("LEGEND_COLORS"):
            color_list = buf.split('\t')[1:]
            color_dr_cat_dict = dict(zip(color_list, dr_cat_list))
        # example:
        # DATA
        # NBCDC45	#f68e51
        if buf.startswith("DATA") and not buf.startswith("DATASET"):
            dr_cat_matrix = []
            buf = f_dr_cat.readline()
            while buf:
                buf = buf.strip()
                line = buf.split('\t')
                # add dr_cat based on color
                line.append(color_dr_cat_dict[line[1]])
                dr_cat_matrix.append(line)
                buf = f_dr_cat.readline()
        # Read next line
        buf = f_dr_cat.readline()

dr_cat_df = pd.DataFrame(dr_cat_matrix, columns = ["sample", "color", "dr_cat"])
dr_cat_df = dr_cat_df.set_index("sample")
dr_cat_df.to_csv(dr_cat_outfile, header=True, index=True, sep='\t')

# 3. Individual drug resistance
with open(dr_indiv_file) as f_dr_indiv:
    buf = f_dr_indiv.readline()
    while buf:
        buf = buf.strip()
        # example: FIELD_LABELS	rifampicin	isoniazid	ethambutol	pyrazinamide	streptomycin
        # fluoroquinolones	aminoglycosides	kanamycin	amikacin	capreomycin	ethionamide	para-aminosalicylic_acid
        # clofazimine	linezolid	bedaquiline	delamanid
        if buf.startswith("FIELD_LABELS"):
            dr_list = buf.split('\t')[1:]
            # replace drug name "para-aminosalicylic_acid" with "PAS" since it's too long
            dr_list = ["PAS" if x == "para-aminosalicylic_acid" else x for x in dr_list]
        # example:
        # DATA
        # NBCDC45	#f68e51
        if buf.startswith("DATA") and not buf.startswith("DATASET"):
            dr_indiv_matrix = []
            buf = f_dr_indiv.readline()
            while buf:
                buf = buf.strip()
                line = buf.split('\t')
                dr_indiv_matrix.append(line)
                buf = f_dr_indiv.readline()
        # Read next line
        buf = f_dr_indiv.readline()

dr_indiv_df = pd.DataFrame(dr_indiv_matrix, columns=["sample"]+dr_list)
dr_indiv_df = dr_indiv_df.set_index("sample")
dr_indiv_df.to_csv(dr_indiv_outfile, header=True, index=True, sep='\t')