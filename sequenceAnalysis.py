#!/usr/bin/env python3
import sys

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence



class NucParams:
    '''
    NuclParams: a class that takes in a string of nucleotides and returns the composition of the 
    nucleotides, the composition of the codons, and the composition of the amino acids
    '''

    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        '''
        Constructor: initializes the dictionaries for the composition of the nucleotides, codons, and amino acids
        inString: a string, the nucleotide sequence, default is an empty string
        '''
        
        self.aaComp = {}
        self.nucComp = {}
        self.codonComp = {}

        #initialize the dictionaries
        for aa in NucParams.rnaCodonTable.values():
            self.aaComp[aa] = 0
        for nucleotide in ['A','C','G','T','U', 'N']:
            self.nucComp[nucleotide] = 0
        for RNAcodon in NucParams.rnaCodonTable.keys():
            self.codonComp[RNAcodon] = 0

        #add the sequence to the dictionary
        if len(inString) > 0:
            self.addSequence(inString)

        
    def addSequence (self, inSeq):
        '''
        inSeq: a string, the nucleotide sequence
        '''
        #add the sequence to the nucleotide dictionary
        for nucleotide in inSeq:
            #if the nucleotide is in the dictionary, add 1 to the count
            if nucleotide in self.nucComp:
                self.nucComp[nucleotide] += 1

        #add the sequence to the codon dictionary
        for i in range(0, len(inSeq), 3):
            #get the frame of the sequence
            codon = inSeq[i:i+3]
            #replace T with U, read it as RNA
            codon = codon.replace('T', 'U')
            #if the codon is in the dictionary, add 1 to the count
            if codon in self.codonComp:
                self.codonComp[codon] += 1
        

        #add the sequence to the amino acid dictionary
        for i in range(0, len(inSeq), 3):
            #get the frame of the sequence
            codon = inSeq[i:i+3]
            #replace T with U, read it as RNA
            codon = codon.replace('T', 'U')
            #if the codon is in the dictionary, add 1 to the count
            if codon in NucParams.rnaCodonTable:
                aa = NucParams.rnaCodonTable[codon]
                self.aaComp[aa] += 1

    def aaComposition(self):
        '''
        return: a dictionary of the composition of the amino acids
        '''
        return self.aaComp
    
    def nucComposition(self):
        '''return: a dictionary of the composition of the nucleotides'''
        return self.nucComp
    
    def codonComposition(self):
        '''return: a dictionary of the composition of the codons'''
        return self.codonComp
    
    def nucCount(self):
        '''return: the number of valid nucleotides in the sequence'''
        return sum(self.nucComp.values())
    




 

class ProteinParam :
    '''
    protein_string: a string, the amino acid sequence
    '''
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34


    #self defined variables
    protein_string = ''
    protein_dict = {}
    raw_string = ''

    def __init__ (self, protein):
        '''
        protein: a string, the amino acid sequence
        '''
        
        # raw_string is the original string
        self.raw_string = protein

        # add all the amino acids to the dictionary, and initialize the count to 0
        for amino_acid in self.aa2mw.keys():
            self.protein_dict[amino_acid] = 0

        # add all the amino acids to the protein_string, and count the number of each amino acid
        for character in protein:
            if character.upper() in self.aa2mw.keys():
                self.protein_string += character.upper()
                self.protein_dict[character.upper()] += 1


    def aaCount (self):
        '''
        Return: the number of amino acids in the string
        '''
        return len(self.protein_string)

    def pI (self, precision = 2):
        '''
        precision: an integer, the number of decimal places to round to (default 2)

        Return: the theoretical isoelectric point of the protein
        '''
        

        # create a list of pH values
        pH_list = [0, 14]

        # find the pH value that gives the 0 net charge, binary search

        # while the difference between the two pH values is greater than 10 ^ (-precision)
        while pH_list[1] - pH_list[0] > 10 ** (-precision):

            # calculate the mid point
            mid = (pH_list[0] + pH_list[1]) / 2

            # determine which bound to replace, based on the sign of the product of the charge at the mid point and the charge at the bound
            if self._charge_(mid) * self._charge_(pH_list[0]) < 0:
                pH_list[1] = mid
            else:
                pH_list[0] = mid
        
        # return the pH value that gives the 0 net charge
        return round((pH_list[0] + pH_list[1]) / 2, precision)

    def aaComposition (self) :
        '''
        Return: a dictionary of each amino acid and its count
        '''
        return self.protein_dict

    def _charge_ (self, pH):
        '''
        _charge_: a private method that calculates the charge of the protein at a given pH
        pH: a float, the pH value

        Return: the net charge of the protein at the given pH
        '''
        # calculate the positive charge, add the charge of the N-terminus first
        positiveSum = (10 ** (self.aaNterm)) / (10 ** (self.aaNterm) + 10 ** (pH))

        # calculate the negative charge, add the charge of the C-terminus first
        negativeSum = (10 ** (pH)) / ((10 ** (self.aaCterm)) + (10 ** (pH)))

        # add the charge of each amino acid, only count those that affect the charge positively
        for amino_acid in self.aa2chargePos:
            positiveSum += ((self.protein_dict[amino_acid] * (10 ** self.aa2chargePos[amino_acid])) / ((10 ** (self.aa2chargePos[amino_acid])) + (10 ** (pH))))
        
        # add the charge of each amino acid, only count those that affect the charge negatively
        for amino_acid in self.aa2chargeNeg:
            negativeSum += ((self.protein_dict[amino_acid] * (10 ** (pH))) / (10 ** (self.aa2chargeNeg[amino_acid]) + (10 ** (pH))))
            
        # return the net charge
        return positiveSum - negativeSum


    def molarExtinction (self, Cystine = True):
        '''
        Cystine: a boolean, whether to include Cystine in the calculation

        Return: the molar extinction coefficient of the protein
        '''

        molar_sum = 0
        
        # add the molar extinction coefficient of each amino acid, only count those that have a molar extinction coefficient
        for character in self.protein_string:

            # if Cystine is False, do not count Cystine
            if character in self.aa2abs280.keys():
                if character == 'C' and Cystine == False:
                    continue
                molar_sum += self.aa2abs280[character]
        return molar_sum
            

    def massExtinction (self, Cystine = True):
        '''
        Cystine: a boolean, whether to include Cystine in the calculation
        
        Return: the mass extinction coefficient of the protein
        '''

        # calculate the molecular weight of the protein
        myMW =  self.molecularWeight()

        # return the mass extinction coefficient, divide by the molar extinction coefficient if the molecular weight is not 0
        return self.molarExtinction(Cystine) / myMW if myMW else 0.0

    def molecularWeight (self):
        '''
        Return: the total molecular weight of the protein
        '''

        myWeight = 0

        # add the molecular weight of each amino acid, minus the weight of H2O from peptide bond formation
        for character in self.protein_string:
            myWeight += (self.aa2mw[character] - self.mwH2O)
        
        # add the weight of H2O, return the value
        return myWeight + self.mwH2O


'''
ORFfinder will go through a DNA sequence and find all the open reading frames (ORFs) in the sequence.
In each frame
It will first go through the sequence and find all the start codons, then find the stop codons that follow the start codons.

'''
class ORFfinder:
    """ Class with method to return list of all ORFs from a DNA string. """

    def __init__(self, seq):
        """ 
        Constructor for ORFFinder. 
        Parameters:
            self.seq: DNA sequence
        """
        self.validStart = ['ATG']
        self.validStop = ['TAG', 'TGA', 'TAA']
        self.validORFs = []
        self.seq = seq

    ##Front overhangs always exist if there exist stop codons. 
    ##Back overhangs always exist if there exist start codons.
        
    def OrfFinder(self, validStartCodons = [], validStopCodons = [], minLen = 0, onlyLongest = True, rev=False):
        '''
        
        '''
        # If the validStartCodons and validStopCodons are not given, use the default values
        if validStartCodons:
            self.validStart = validStartCodons
        if validStopCodons:
            self.validStop = validStopCodons

        #List for stop and start codons
        start_codon_pos = []
        stop_codon_pos = []
        start_codon_pos_index = 0
        stop_codon_posIndex = 0

        #For each frame, find the start and stop codons
        for frame in range(3):
            #find if the start codon is found
            startCodonFound = False
            stopCodonFound = False
            for pos in range(frame, len(self.seq), 3):
                codon = self.seq[pos:pos + 3]
                if codon in self.validStart:
                    start_codon_pos.append(pos)
                    startCodonFound = True
                if codon in self.validStop:
                    stop_codon_pos.append(pos)
                    stopCodonFound = True

                #If the start codon is found, and the stop codon is found, add the ORF to the list