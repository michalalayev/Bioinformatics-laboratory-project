# Plot the cumulative distribution of reads in each cell, and select the “knee” of the distribution as the estimated 
# number of cells that were sequenced.
# The input file should contain the number of transcripts per each cell tag, ordered from max to min, 
# named 'out_cell_readcounts.txt.gz' as created in the Drop-seq processing.
# The limit of the X axis in the plot can be also be controllod by the argument x_axis_lim, which is set to 30000 by default.

plot_cum_dist <- function(out_cell_readcounts_file, x_axis_lim=30000) {
  data<-read.table(out_cell_readcounts_file, header=F, stringsAsFactors=F) 
  x<-cumsum(data$V1) # claculate comulative count of transcripts
  x<-x/max(x) # the growing frequencies, as y axis
  plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
       ylab="cumulative fraction of reads", main="Cumulative distribution of reads per cell", xlim=c(1,x_axis_lim))
  grid()
}
