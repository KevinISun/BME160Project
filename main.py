
import pysam 
import matplotlib.pyplot as plt
from sequenceAnalysis import FastAreader 

class depthGraphGenerator():


	def GenCoverage(self, ref_genome, bamfile):
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
		
		return coverage

	def GenHistogram(self, posbp_read, N, patches):
		''' Generates a histogram graph of the depth of the genome. ''' 
		# Defines x as the dictionary keys 
		x = posbp_read.keys()
		# Defines y as the dictionary values 
		y = posbp_read.values()
		# Defines the legend 
		legend_key = ['distribution of the reads in the genome']
		# Plotting a basic histogram
		plt.hist(x, y, bins=100, color='skyblue', edgecolor='black')
 		# Setting color
		fracs = ((N**(1 / 5)) / N.max())

		# Normalize the data to 0-1
		norm = colors.Normalize(fracs.min(), fracs.max())
 
		for thisfrac, thispatch in zip(fracs, patches):
    			color = plt.cm.viridis(norm(thisfrac))
    			thispatch.set_facecolor(color)
 
		# Labels the x-axis as the nucleotide position
		plt.xlabel('Nucleotide Position')
		# Labels the y-axis as the number of reads in that position
		plt.ylabel('Number of Reads')
		# Labels the histogram with the title
		plt.title('Read Depth of Genome')
		# Adds the legend to the histogram
		plt.legend(legend_key)
 
		# Display the plot
		plt.show()
	
	

if __name__ == '__main__':
	bamfile = pysam.AlignmentFile('alignment_1.bam', "rb")
	refGenome = FastAreader('GCA_009858895.3.fasta')
	depthGraph = depthGraphGenerator()
	coverage = depthGraphGenerator.GenCoverage(refGenome, bamfile)
	print(coverage)
	posbp_read = depthGraphGenerator.posbp_read
	depthGraphGenerator.GenHistogram(posbp_read, N, patches)

	
