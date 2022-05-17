import bamnostic as bs
import pandas as pd
from Bio import SeqIO
import pickle
from ast import literal_eval
import sys
pd.options.mode.chained_assignment = None  # default='warn'


gene_name = sys.argv[1]
bam_file = sys.argv[2]  # X_no_indels_mm_sorted.bam
fasta_file = sys.argv[3]  # "dd_Smed_v6.fasta"
file_prefix = ''


def set_file_prefix(gene_name):
    num = gene_name.split('_')[3]
    global file_prefix
    file_prefix = gene_name + "/" + num + "_"


def parse_bam_into_df(bam_file):
    print("Parsing BAM file '{}' into a dataframe".format(bam_file))
    bam = bs.AlignmentFile(bam_file, 'rb')
    df = pd.DataFrame(columns=['name', 'pos', 'seq', 'qual', 'tag', 'NM'])
    name_list, pos_list, seq_list, qual_list, tag_list, NM_list = [], [], [], [], [], []
    i = 0
    for read in bam:
        name_list.append(read.read_name)
        pos_list.append(read.pos + 1)  # as written in the read
        seq_list.append(read.query_sequence)  # sequence of the read that was aligned
        qual_list.append(list(read.query_qualities))  # an array of sequencing quality per base
        tag_list.append(read.tags['XC'][1])  # cell tag
        NM_list.append(read.tags['NM'][1])  # number of mismatches
        i += 1
        if i % 100000 == 0:
            print("{} reads were processed".format(i))
    print("Done processing the reads, the dataframe was created")
    df['name'] = name_list
    df['pos'] = pos_list
    df['seq'] = seq_list
    df['qual'] = qual_list
    df['tag'] = tag_list
    df['NM'] = NM_list
    print("Head of the created dataframe:")
    print(df.head())
    output_file_name = file_prefix + "bam_as_df.csv.gz"
    df.to_csv(output_file_name, index=False, compression="gzip")
    print("The dataframe was saved to '{}'".format(output_file_name))


def sample_reads(gene_name, bam_file, amount):
    print("Starting sampling reads")
    df = pd.read_csv(file_prefix + "bam_as_df.csv.gz", compression='gzip', header=0)
    cnt = df['tag'].value_counts() # count number of distinct cell tags
    cnt = cnt.reset_index()  # converts from series to dataframe
    cnt.rename(columns={'index': 'tag', 'tag': 'count'}, inplace=True)

    # divide the dataframe to cells that have up to 'amount' reads, and cells that have more reads:
    cells_no_sampling = cnt.loc[cnt['count'] <= amount]['tag']
    cells_to_sample = cnt.loc[cnt['count'] > amount]['tag']
    reduced_df = df.loc[df['tag'].isin(cells_no_sampling)]  # all reads from cells that have reads up to 'amount'

    # sample 'amount' reads from cells that have more reads this that:
    df_for_sampling = df.loc[df['tag'].isin(cells_to_sample)]
    result = df_for_sampling.groupby('tag').sample(n=amount, random_state=1)

    reduced_df = pd.concat([result, reduced_df]) # combine all selected reads
    output_csv = "bam_as_df_sample.csv.gz"
    reduced_df.to_csv(file_prefix + output_csv, index=False, compression="gzip")
    print("Done sampling")
    print("The sampled reads were saved to '{}'".format(file_prefix + output_csv))

    # save statistics of the gene in a file:
    txt_output_file_name = file_prefix + "BAM_file_stats.txt"
    f = open(txt_output_file_name, "w")
    f.write("Gene '{}' STATS, in BAM file '{}':\n".format(gene_name, bam_file))
    f.write("number of reads (w/o indels, only matches or mismatches) in BAM file: {}\n".format(len(df)))
    f.write("number of distinct cells in BAM file: {}\n".format(len(cnt)))
    f.write("number of cells that have up to {} reads: {}\n".format(amount, len(cells_no_sampling)))
    num = sum(cnt.loc[cnt['count'] <= amount]['count'])
    f.write("number of reads in these cells: {}\n".format(num))
    f.write("number of cells that have more then {} reads and need to be sampled from: {}\n".format(amount, len(cells_to_sample)))
    f.write("number of reads in these cells: {}\n".format(len(df_for_sampling)))
    f.write("number of reads sampled from these cells: {}\n".format(len(result)))
    f.write("total number of reads in '{}': {}\n".format(output_csv, len(reduced_df)))
    f.close()
    print("Created stat file for the gene named '{}', that reads:".format(txt_output_file_name))
    print(open(txt_output_file_name, "r").read())


def count_mm_per_pos(gene_name, fasta_file):
    print("Starting to count mismatches per position, including only good quality mismatches")
    df = pd.read_csv(file_prefix + "bam_as_df_sample.csv.gz", compression='gzip', header=0)
    mm_df = df.loc[df['NM'] > 0]
    # the 'seq_qual' column supposed to be list but read as str, so we convert it:
    mm_df['qual'] = mm_df['qual'].apply(literal_eval)

    # add some more information to the gene stat file:
    txt_output_file_name = file_prefix + "BAM_file_stats.txt"
    f = open(txt_output_file_name, 'a')
    f.write("In the sample file:\n")
    f.write("number of reads with mismatches: {}\n".format(len(mm_df)))
    f.write("contained within {} distinct cells.\n".format(mm_df.tag.nunique()))
    f.close()
    print("Updated stat file for the gene named {}, that now reads:".format(txt_output_file_name))
    print(open(txt_output_file_name, "r").read())

    # find the length of the gene sequence:
    record_dict = SeqIO.index(fasta_file, "fasta")
    gene_record = record_dict[gene_name]
    gene_seq_len = len(gene_record.seq)

    # create a dictionary where keys are positions in the gene reference sequence, values are data of good quality mismatches at each position
    d = {key: [] for key in range(1, gene_seq_len+1)}
    print("Number of reads to process:", len(mm_df))
    for index, read in mm_df.iterrows():
        seq = read['seq']
        qual = read['qual']
        for char_pos, char in enumerate(seq):
            if char != '=' and qual[char_pos] >= 20:
                d[char_pos + read['pos']].append((char, qual[char_pos], read['tag']))
    print("Done processing the reads")

    # save dict d to a pkl file:
    dict_file_name = file_prefix + "sample_pos_good_quality_mm_dict.pkl"
    d_file = open(dict_file_name, "wb")
    pickle.dump(d, d_file)
    d_file.close()
    print("The results were saved to '{}'".format(dict_file_name))

    print("Checking results:")
    file = open(dict_file_name, "rb")
    dict = pickle.load(file)
    print("first 50 positions: ", list(dict.items())[:50])
    file.close()


def save_sample_read_names_to_file():
    df = pd.read_csv(file_prefix + "bam_as_df_sample.csv.gz", compression='gzip', header=0)
    names = df['name']
    output_file_name = file_prefix + "read_names_to_use.csv"
    names.to_csv(output_file_name, index=False, header=False)
    print("Names of reads that were sampled were saved to '{}'".format(output_file_name))


def do_work_for_gene(gene_name):
    if len(sys.argv) != 4:
        print("There should be exactly 3 arguments, but {} were passed. Exiting.".format(len(sys.argv)))
        exit()
    set_file_prefix(gene_name)
    parse_bam_into_df(bam_file)
    sample_reads(gene_name, bam_file, 10)
    count_mm_per_pos(gene_name, fasta_file)
    save_sample_read_names_to_file()


if __name__ == '__main__':
    do_work_for_gene(gene_name)
