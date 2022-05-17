from Bio import SeqIO

def create_gtf(fasta_file, gtf_file_name):
    print("Creating GTF file for the fasta file: ", fasta_file)
    f = open(gtf_file_name, "w")
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        id = seq_record.id
        record = '''{}\tsmed\texon\t1\t{}\t.\t+\t.\tgene_id "{}"; transcript_id "{}"; exon_number "1"; gene_name "{}"; gene_biotype "{}"; transcript_name "{}"; exon_id "{}";\n'''.format(
                id, len(seq_record), id, id, id, id, id, id)
        f.write(record)
    f.close()
    print("The GTF file '{}' was created".format(gtf_file_name))

if __name__ == '__main__':
    fasta_file = "dd_Smed_v6_added.fasta"
    gtf_file_name = "dd_Smed_v6_added.gtf"
    create_gtf(fasta_file, gtf_file_name)

# example for gtf record:
# dd_Smed_v6_1_0_1    smed    exon    1   416 .   +   .   gene_id "dd_Smed_v6_1_0_1"; transcript_id "dd_Smed_v6_1_0_1"; exon_number "1"; gene_name "dd_Smed_v6_1_0_1"; gene_biotype "dd_Smed_v6_1_0_1"; transcript_name "dd_Smed_v6_1_0_1"; exon_id "dd_Smed_v6_1_0_1";

