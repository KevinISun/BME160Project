import brewtools
import pysam 
import matplotlib.pyplot as plt
import fastAreader

class Main(self):
	self.ref_genome = refGenome
	self.bamfile = bamfile 

	def GenCoverage(self):
		posbp_read = {}
		for i in range(1, len(ref_genome)+1):
			readCount = 0
			for read in bamfile.fetch('ENA|MN908947|MN908947.3', i, i+1):
				readCount += 1
				posbp_read[i] = readCount
			bamfile.close()
			
		coverage = readCount/(len(ref_genome)+1)

	def GenHistogram(self):
		# graphing a histogram of the depth of the genome 
		# x-axis is the base pair position
		# y-axis is number of reads in that position
		# Plotting a basic histogram
        		plt.hist(i, readCount, bins=100, color='skyblue', edgecolor='black')
 
        		# Adding labels and title
        		plt.xlabel('Number of Reads')
       		plt.ylabel('Nucleotide Position')
        		plt.title('Read Depth of Genome')
 
        		# Display the plot
        		plt.show()
	def main(self):
		bamfile = pysam.AlignmentFile("alignment_1.bam", "rb")
		refGenome = fastAreader(GCA...fasta)
		


