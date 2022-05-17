import pandas as pd


def find_expression_frequency(dge_file):
    dge_df = pd.read_csv(dge_file, delimiter="\t", index_col=0) # the index is the names of the genes
    num_of_genes = len(dge_df.columns)
    # calculate frequency of cells that express each gene:
    dge_df[">0 freq"] = (dge_df.gt(0).sum(axis=1))/num_of_genes*100
    dge_df["num of cells"] = dge_df[dge_df.columns[:-1]].gt(0).sum(axis=1)
    dge_df["num of reads"] = dge_df[dge_df.columns[:-2]].sum(axis=1)
    dge_df = dge_df.sort_values(">0 freq") # sort by frequency
    print("All the genes with frequency of cells that express them:")
    print(dge_df.loc[:, dge_df.columns[-3:]])
    high_exp_df = dge_df.loc[dge_df[">0 freq"] >= 20] # genes that expressed in at least 20% of the cells
    print("Genes that at least 20% of the cells express them:")
    print(high_exp_df.loc[:, dge_df.columns[-3:]])


def DGE_sanity_check(dge_file, type1_markers, type2_markers):
    dge_df = pd.read_csv(dge_file, delimiter="\t", index_col=0)
    num_of_cells = len(dge_df.columns)
    type1_co = check_coexpression(dge_df, type1_markers)
    print("Number of cells to co-express the genes {} is: {}".format(type1_markers, len(type1_co.columns)))
    print("which is {}% of the cells".format(round(len(type1_co.columns)/num_of_cells*100, 2)))
    type1_and_type2_co = check_coexpression(dge_df, type1_markers+type2_markers)
    print("Number of cells to co-express the genes {} is: {}".format(type1_markers+type2_markers, len(type1_and_type2_co.columns)))
    print("which is {}% of the cells".format(round(len(type1_and_type2_co.columns) / num_of_cells * 100), 2))
    print("and {}% of the cells that co-express the genes {}".format(round(len(type1_and_type2_co.columns)/len(type1_co.columns)*100 , 2), type1_markers))


def check_coexpression(dge_df, markers):
    df_markers = dge_df.loc[markers, :]  # only the rows of the markers
    filter = (df_markers > 0).all()  # returns a series of bools, that represent if each col contains only values > 0 or not.
    df_co = df_markers.loc[:, filter] # only the columns of cells that co-express the markers
    return df_co


if __name__ == '__main__':
    dge_file = "bowtie_gene_exon_tagged.dge.txt"
    find_expression_frequency(dge_file)
    epidermal = ['dd_Smed_v6_69_0_1', 'dd_Smed_v6_332_0_1']
    muscle = ['dd_Smed_v6_323_0_1']
    DGE_sanity_check(dge_file, epidermal, muscle)



