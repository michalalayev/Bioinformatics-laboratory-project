import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
pd.options.mode.chained_assignment = None  # default='warn'

output_file = "all_genes_mm_per_pos_info.csv"
genes_file = sys.argv[1]  # "genes_to_use.txt"


def create_mm_info_for_all_genes():
    print("Creating mismatches info file for all genes combined.")
    df = pd.read_csv(genes_file, names=['gene_names'])
    print("The genes to use:")
    print(df)
    mm_p = "good quality mismatch %"
    cols = ['gene_name', 'position', 'coverage', 'reference', mm_p]
    nucs = ['A', 'T', 'C', 'G']
    nucs_p = [nuc + " /total reads (%)" for nuc in nucs]
    cols = cols + nucs + nucs_p
    result = pd.DataFrame(columns=cols)

    for gene_name in df['gene_names']:
        gene_num = gene_name.split('_')[3]
        info_file_name = gene_name + "/" + gene_num + "_sample_good_quality_mismatch_info.csv"
        gene_info = pd.read_csv(info_file_name, header=0)
        tmp = pd.DataFrame(columns=cols)
        for col in cols[1:]:
            tmp[col] = gene_info[col]
        tmp['gene_name'] = gene_name
        result = pd.concat([result, tmp])

    print("Head of the created file:")
    print(result.head())
    print("total length (number of positions):", len(result))
    result.to_csv(output_file, index=False)
    print("Saving file as '{}'".format(output_file))


def plot_distribution(cutoff=0):
    print("Creating distribution of good quality mismatch % per position in all of the genes")
    df = pd.read_csv(output_file, header=0)
    print("The content of the file 'all_genes_mm_per_pos_info.csv':")
    print(df)
    mm_col = "good quality mismatch %"
    if cutoff: # include only values >= the cutoff
        df = df.loc[df[mm_col] >= cutoff]
    mm_precentage = df[mm_col]

    # plot and save as png:
    sns.set_theme()
    sns.set(font_scale=1.7)
    ax = sns.displot(mm_precentage, kind="kde")
    plt.xlabel("good quality mismatch % per position")
    plt.title("Distribution of good quality mismatch percentage per position (number of values: {})".format(len(df)))
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    plot_file_name = 'all genes- distribution of mismatch percentage per position - cutoff {}.png'.format(cutoff)
    plt.savefig(plot_file_name, dpi=100)
    print("The plot was saved as '{}'.".format(plot_file_name))
    plt.show()


def get_positions_by_percentage_range(min, max):
    print("extracting positions with mismatch percentage of {} to {}".format(min, max))
    df = pd.read_csv(output_file, header=0)
    range_df = df.loc[(df["good quality mismatch %"] >= min) & (df["good quality mismatch %"] <= max)]
    range_df = range_df.sort_values(by=["good quality mismatch %"])
    file_name = "mm_percentage_range_{}-{}.csv".format(min, max)
    range_df.to_csv(file_name)
    print("results were saved to '{}'".format(file_name))
    print("The results:")
    print(range_df)
    print("Number of positions found: {}".format(len(range_df)))


# this function was not used. It meant to find positions where there are at least 2 other nucleotides besides the
# reference base, each appears with frequency of at least 'thresh'.
def find_positions_with_3_nucs(thresh):
    df = pd.read_csv(output_file, header=0)
    nucs = ['A', 'T', 'C', 'G']
    nucs_p = [nuc + " /total reads (%)" for nuc in nucs]
    nucs_thresh = [nuc + ' >= ' + str(thresh) for nuc in nucs]
    new_df = df.loc[df["good quality mismatch %"] < (100-thresh)] # so at least thresh % is the reference nucleotide
    new_df.reset_index(inplace=True, drop=True)
    for i in range(4):
        new_df[nucs_thresh[i]] = new_df[nucs_p[i]].apply(lambda x: x >= thresh)
    col = 'how many >= ' + str(thresh)
    new_df[col] = new_df.iloc[:, -4:].sum(axis=1)
    relev_df = new_df.loc[new_df[col] >= 2] # in the nuc col that matches the reference will be 0, so having two columns
                                            # that have mm % >= thresh means we have at lest 3 different nucs
    relev_df.drop(columns=nucs_thresh, inplace=True)
    print("The found positions:")
    print(relev_df)
    print("Amount: ", len(relev_df))
    file_name = "positions_with_at_least_3_nucs_thresh_{}.csv".format(thresh)
    relev_df.to_csv(file_name)
    print("results were saved to '{}'".format(file_name))



if __name__ == '__main__':
    create_mm_info_for_all_genes()
    plot_distribution()  # initially all values
    plot_distribution(2)  # values >= 2
    # get_positions_by_percentage_range(93, 100)
    # get_positions_by_percentage_range(40, 60)
    # get_positions_by_percentage_range(15, 35)
    # get_positions_by_percentage_range(75, 82)
    # find_positions_with_3_nucs(5)



