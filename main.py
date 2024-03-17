
import pysam 
import matplotlib.pyplot as plt
from sequenceAnalysis import FastAreader 

class Main(self):
	self.ref_genome = refGenome
	self.bamfile = bamfile 

	def GenCoverage(self):
		''' Generates the coverage of the genome. '''
		# Initializes dictionary for the storage of the nucleotide position from the reference and the number of reads aligned in that position
		posbp_read = {}
		# For loop iterates through each nucleotide position in the reference genome 
		for i in range(1, len(ref_genome)+1):
			# Initializes the readCount of the reference genome
			readCount = 0
			# For loop iterates through the .bam file to find the number of reads 
			for read in bamfile.fetch('ENA|MN908947|MN908947.3', i, i+1):
				# As each read is found in the .bam file it is added to the readCount
				readCount += 1
				# Creates a dictionary that stores the nucleotide position from the reference genome as the key 
				# and the number of reads aligned in that position as the value
				posbp_read[i] = readCount
			bamfile.close()
			
		#calculates the coverage of the genome
		coverage = readCount/(len(ref_genome))

	def GenHistogram(self):
		''' Generates a histogram graph of the depth of the genome. ''' 
	
		# Plotting a basic histogram
		plt.hist(i, readCount, bins=100, color='skyblue', edgecolor='black')
 
		# Labels the x-axis as the nucleotide position
		plt.xlabel('Nucleotide Position')
		# Labels the y-axis as the number of reads in that position
		plt.ylabel('Number of Reads')
		# Labels the histogram with the title
		plt.title('Read Depth of Genome')
 
		# Display the plot
		plt.show()

	def main(self):
		bamfile = pysam.AlignmentFile("alignment_1.bam", "rb")
		refGenome = FastAreader(GCA_009858895.3.fasta)
	


