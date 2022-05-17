import pickle
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def see(pos):
    file = open("pos_good_quality_mm_dict.pkl", "rb")
    "sample_pos_good_quality_mm_dict.pkl"
    good_dict = pickle.load(file)
    print(f"looking at position {pos} in the gene:")
    print(good_dict[pos][:10])
    print(len(good_dict[pos]))
    file = open("sample_pos_good_quality_mm_dict.pkl", "rb")
    good_dict = pickle.load(file)
    print(f"looking at position {pos} in the gene:")
    print(good_dict[pos][:10])
    print(len(good_dict[pos]))


def explore_mm_by_cells_in_position(pos_mm_dict_file, pos):
    file = open(pos_mm_dict_file, "rb")
    good_dict = pickle.load(file)
    print("looking at position {} in the gene:".format(pos))
    print("first 10 items in the dictionary {}:\n".format(pos_mm_dict_file), good_dict[pos][:10])
    print("number of mismatches found at this position: ", len(good_dict[pos]))

    tags_mm_dict = {}  # keys are distinct cell tags, and value is list of tuples (char, qual) of mismatches at the specified position in each cell
    for record in good_dict[pos]:
        cell_tag = record[2]
        if cell_tag in tags_mm_dict.keys():
            tags_mm_dict[cell_tag].append((record[0], record[1]))
        else:
            tags_mm_dict[cell_tag] = [(record[0], record[1])]
    print("dictionary of (mismatch, qual) per each cell tag, first 5 items:\n", list(tags_mm_dict.items())[:5] )

    cell_mm_freq = {}  # keys are distinct cell tags, and value is dict that counts the number of each nucleotide
                       # appearing as mismatch at the specific position in the cell
    for cell in tags_mm_dict.keys():
        mm_list = tags_mm_dict[cell]
        nuc_count = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
        for mm in mm_list:
            nuc_count[mm[0]] += 1
        # nuc_count.update((x, y/len(mm_list)*100) for x, y in nuc_count.items()) # convert from count to frequency
        nuc_count['total'] = len(mm_list)
        cell_mm_freq[cell] = nuc_count
    print("mismatch count of each nucleotide in each cell, first 10 items:")
    print(list(cell_mm_freq.items())[:10])
    print("number of cells:", len(cell_mm_freq))
    # other_nuc_than_c = [item for item in cell_mm_freq.items() if item[1]['C'] != item[1]['total']]
    # print("number of cells that has not only C:", len(other_nuc_than_c))

    # save the results to a file:
    nucs = ['A', 'T', 'C', 'G']
    df = pd.DataFrame(index=cell_mm_freq.keys(), columns=nucs)
    for nuc in nucs:
        df[nuc] = [dict_tup[1][nuc] for dict_tup in list(cell_mm_freq.items())]
    res_file = "position_{}-mm_count_per_cell.csv".format(pos)
    df.to_csv(res_file)
    print("Created a file named '{}' containing the count of each nucleotide that appeared as mismatch at the position, "
          "for each cell separately".format(res_file))


def plot_heatmap(pos):
    file = "position_{}-mm_count_per_cell.csv".format(pos)
    print("Using the previously created file {} for generating a heatmap".format(file))
    df = pd.read_csv(file, index_col=0)

    # comput the matrix height in points and inches
    matrix_height_pt = 0.5 * df.shape[0]
    dpi = 72.27
    matrix_height_in = matrix_height_pt / dpi

    # compute the required figure height
    top_margin = 0.04  # in percentage of the figure height
    bottom_margin = 0.04  # in percentage of the figure height
    figure_height = matrix_height_in / (1 - top_margin - bottom_margin)

    # build the figure instance with the desired height
    fig, ax = plt.subplots(figsize=(15, figure_height), gridspec_kw=dict(top=1 - top_margin, bottom=bottom_margin))

    # create the heatmap using seaborn
    ax = sns.heatmap(df, cmap='YlOrBr', ax=ax)
    v_fig_file = 'position {} heatmap (vertical).png'.format(pos)
    plt.savefig(v_fig_file)

    # my attempt for ploting the heatmap:
    plt.clf()
    df = df.transpose()
    fig, ax = plt.subplots(figsize=(25, 5))
    sns.heatmap(df, cmap='YlOrBr', ax=ax)
    plt.xticks(fontsize='xx-small')
    plt.tight_layout()
    h_fig_file = 'position {} heatmap (horizontal).png'.format(pos)
    plt.savefig(h_fig_file)

    print("Two heatmaps were saved: '{}' and '{}'".format(v_fig_file, h_fig_file))


if __name__ == '__main__':
    position = 808
    pos_mm_dict = "sample_pos_good_quality_mm_dict.pkl"
    explore_mm_by_cells_in_position(pos_mm_dict, position)
    plot_heatmap(position)