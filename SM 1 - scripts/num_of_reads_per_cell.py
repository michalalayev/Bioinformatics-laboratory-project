import pickle
import bamnostic as bs
import matplotlib.pyplot as plt


def count_reads_per_cell(gene_name, bam_file):
    bam = bs.AlignmentFile(bam_file, 'rb')
    # create dictionary where keys are cell tags and values are the number of reads found in each cell
    reads_per_cell_count = dict()
    i = 0
    for read in bam:
        cell_tag = read.tags['XC'][1]
        if cell_tag in reads_per_cell_count.keys():
            reads_per_cell_count[cell_tag] += 1
        else:
            reads_per_cell_count[cell_tag] = 1
        i += 1
        if i % 100000 == 0:
            print("done {} reads".format(i))

    reads_per_cell_count = list(reads_per_cell_count.items()) # convert dict to list of tuples to sort it
    reads_per_cell_count_sorted = sorted(reads_per_cell_count, key=lambda x: x[1])  # sort by number of reads per cell
    print("number of cells:", len(reads_per_cell_count_sorted))

    # save list reads_per_cell_count_sorted in a file:
    output_file_name = gene_name+"_reads_per_cell_count_sorted.pkl"
    output_file = open(output_file_name, "wb")
    pickle.dump(reads_per_cell_count_sorted, output_file)
    output_file.close()
    print("The count of the number of reads per cell was saved as a list in '{}'".format(output_file_name))
    print("First 20 items in the list:")
    print(reads_per_cell_count_sorted[:20])


def plot_reads_per_cell(gene_name):
    count_file = open(gene_name+"/reads_per_cell_count_sorted.pkl", "rb")
    reads_per_cell_count_sorted = pickle.load(count_file)
    # plot read per cell:
    fig = plt.gcf()
    fig.set_size_inches(18.5, 10.5, forward=True)
    x = range(len(reads_per_cell_count_sorted))
    y = [(tup[1]) for tup in reads_per_cell_count_sorted]
    plt.plot(x, y)
    plt.yscale("log")
    plt.grid()
    size = 16
    plt.xlabel("cells (sorted by number of reads)", fontsize=size)
    plt.ylabel("number of reads (log10 scale)", fontsize=size)
    plt.title("Gene {} - Number of reads in each cell".format(gene_name), fontsize=size)
    plot_file_name = "gene {} number of reads per cell.png".format(gene_name)
    plt.savefig(plot_file_name, dpi=100)
    print("The plot was saved as '{}'.".format(plot_file_name))
    plt.show()


    # find the percentage of cells that contain up to 10 reads:
    # g = [tup for tup in reads_per_cell_count_sorted if tup[1] <= 10]
    # print(g, len(g))
    # total = len(reads_per_cell_count_sorted)
    # print(total)
    # print(len(g)/total*100)


if __name__ == '__main__':
    gene_name = 'dd_Smed_v6_332_0_1'
    bam_file = '332_no_indels_mm_sorted.bam'
    count_reads_per_cell(gene_name, bam_file)
    plot_reads_per_cell(gene_name)
