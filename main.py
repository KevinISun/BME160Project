
import pysam 
import matplotlib.pyplot as plt
from sequenceAnalysis import FastAreader 

class depthGraphGenerator():

	# Constructor
	def __init__(self):
		self.posbp_read = {}

	# Method to get posbp_read
	def getPosbpRead(self):
		return self.posbp_read



	def GenCoverage(self, bamfile, refgenome_header, refgenome_sequence):
		''' Generates the coverage of the genome. '''
		# Initializes dictionary for the storage of the nucleotide position from the reference and the number of reads aligned in that position
		self.posbp_read = {}
		# For loop iterates through each nucleotide position in the reference genome 
		for i in range(1, len(refgenome_sequence)+1):
			# Initializes the readCount of the reference genome
			readCount = 0
			# For loop iterates through the .bam file to find the number of reads 
			for read in bamfile.fetch('ENA|MN908947|MN908947.3', i, i+1):
				# As each read is found in the .bam file it is added to the readCount
				readCount += 1
				# Creates a dictionary that stores the nucleotide position from the reference genome as the key 
				# and the number of reads aligned in that position as the value
				self.posbp_read[i] = readCount			
		#calculates the coverage of the genome
		coverage = readCount/(len(refgenome_sequence))
		
		return coverage

	def GenHistogram(self, posbp_read, N, patches):
		''' Generates a histogram graph of the depth of the genome. ''' 
		'''
		Parameters:
		posbp_read (dict): A dictionary that stores the nucleotide position from the reference genome as the key
		and the number of reads aligned in that position as the value
		N (int): The number of reads aligned in that position
		patches (int): The number of patches
		'''
		# Defines x as the dictionary keys 
		x = posbp_read.keys()
		# Defines y as the dictionary values 
		y = [posbp_read[key] for key in x]
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
	refGenomeFile = FastAreader('GCA_009858895.3.fasta')
	depthGraph = depthGraphGenerator()
	refGenome_header, refGenome_sequence = refGenomeFile.readFasta()
	coverage = depthGraphGenerator.GenCoverage(bamfile, refGenome_header, refGenome_sequence)
	bamfile.close()
	print(coverage)

	posbp_read = depthGraphGenerator.getPosbpRead()
	depthGraphGenerator.GenHistogram(posbp_read, N, patches)

	
