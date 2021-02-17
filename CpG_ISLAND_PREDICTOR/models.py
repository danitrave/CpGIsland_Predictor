import random

absolutDimerInside = {'A':0,'C':0,'G':0,'T':0}   #the following variables has been inizialized as global
absolutDimerOutside = {'A':0,'C':0,'G':0,'T':0}   #since they are used by different functions in the program
InsideMatrix = {}      #this variable will be provided of the inside model at the end
OutsideMatrix = {}    #this variable will be provided of the outside model at the end
seqForOutside = ''    #use to collaps sequences of chromosom 22 in a unique string, that are not CpG
AvLength = 0                #avarage length for sliding window

def Chr22_CpGLocFinder():     #this function opens a FASTA file containing the chromosome 22
                                                          #and return the locus of each CpG and a string ontaining non-CpG
      fileGenome = open("chr22.fa")     #open the FASTA file
      chr22 = ''.join(fileGenome.readlines()[1:]).replace('\n','') #insert it in a string      
      fileGenome.close()
      
      chr22 = chr22.upper()
      location = {}
      fileLocation = open('CpGLoc22.txt','r')  #open the file containing all staring points of CpG islands
      index = fileLocation.readlines()        #in the FASTA file
      e = 0

      global AvLength      #update its value with the avarage length
      allLengths = []
                  
      for i in index[1:]:           #from each line of the txt file containg the location of CpG
            boundaries = tuple()
            i = i.split()
            if i[0] == 'chr22':
                  AvLength += int(i[3])            #gradually add values at column 3 for the avarage length
                  allLengths += [int(i[3]),]     #in colum 3 we find the length of each CpG in chr22
                  boundaries += (int(i[1]),)+(int(i[2]),)   #add values at column 1 and 2 representing
                  location[e] = boundaries                  #starting end ending index of each CpG in chr22
                  e+=1
            
      AvLength =int( AvLength/len(location))    #get avarage length of CpG in chromosome 22
      
      global seqForOutside     #now built a sequence having length equal to a sequence containig all Cpg
                                               #but having NON-CpG sequences
      for j in allLengths:
            start = random.randint(0,len(chr22))     #select a random point in chr22 and take a sequence of the same length 
            seqForOutside += chr22[start:(start+j)+1].replace('N','')   #of each CpG found. Then collaps each in a unique string

      fileLocation.close()

      return chr22,location  

def CpG():   #this function builds a unique string containing all CpGs
      
      totCpG = Chr22_CpGLocFinder()   #consider the entire chromosome22 and the start-end idex of each CpG
      sequence = totCpG[0]    #chromosome 22
      loc = totCpG[1]          #loci of each CpG
      CpGseq = ''
      for index in loc:   #select each CpG from chromosome 22
            currentCpG = ''
            currentCpG = sequence[loc[index][0]:loc[index][1]+1]   #add it in a unique string
            CpGseq += currentCpG
      return CpGseq    #return the unique string with all CpGs

def absDimerInside(x):   #this function calculate the absolute frequency of alle A,C,G,T 
                                                          #in the sequence containing all CpG island sequences
      global absolutDimerInside
      if x[0] in absolutDimerInside:       #keep on adding 1 as one nucleotide is found
            absolutDimerInside[x[0]] += 1
      
      return absolutDimerInside         #return final absolute frequency
      
def InsideModel():          #this function create the inside model for CpG island

      global InsideMatrix    #dictionary having all possible dimer as keys and the relative frequency of each
                                               #of each dimer as values
      seq = CpG()                    #take the sequence containing all CpG
      
      dimers = {}             #build all possible dimers
      for x in 'ACGT':
            for y in 'ACGT':
                  dimers[x+y] = 0       #inizialize a relative frequence of 0 per each dimer
                  
      for i in range(len(seq)-2):   #scann all CpG sequence every 2 nucletides (dimer)
            absFrequencyIn = absDimerInside(seq[i:i+2])   #take each dimer and consider each nucleotide in the dimer. This is done by the AbsoluteFrequencyInside function
            if seq[i:i+2] in dimers:      #add 1 to the relative frequency of each dimer  
                  dimers[seq[i:i+2]] += 1

      InsideModel = []
      for dimer in dimers:     #calculate the conditional probability of each dime
            value = 0
            if dimer[0] in absFrequencyIn:
                  value = round(dimers[dimer]/absFrequencyIn[dimer[0]],2)  #calculate the conditional probability of the current dimer
                  InsideModel+= [value,]      #add result 

      i = 0
      for h in dimers:      #use the probabilities obtained before and the value of the dimers
            InsideMatrix[h] =InsideModel[i]   #to build the inside model "matrix"
            i+=1
            
      return InsideMatrix #the matrix is a python dictionary having {dimer:probability}

def absDimerOutside(x):   #this function calculate the absolute frequency of A,C,G,T
                                                   #in the sequence not containg CpG islands
      global absolutDimerOutside
      
      if x[0] in absolutDimerOutside:       #as a nucleotide is considered, add 1 to its ABS frequency value
            absolutDimerOutside[x[0]] += 1
      
      return absolutDimerOutside   #return a dictionary containing absolute frequencies
           
def OutsideModel():    #this function creates the outside model for sequences not being a CpG
      
      global seqForOutside
      global OutsideMatrix

      dimers = {}         #creates all dimer of nucleotides
      for x in 'ACGT':
            for y in 'ACGT':
                  dimers[x+y] = 0    #sets their relative frequency to 0

      for k in range(len(seqForOutside)-2):     #consider each dimer in a sequence having streches of DNA outside a CpG
            absFrequencyOut = absDimerOutside(seqForOutside[k:k+2])  #count the absolute frequency of each nucleotide in the sequence
            if seqForOutside[k:k+2] in dimers:       #update continusly the relative frequency of a dimer
                  dimers[seqForOutside[k:k+2]] += 1   #once you find a dimer, add 1

      OutsideModel = []
      for dimer in dimers:     #calculate the conditional probability of each dimer
            value = 0
            if dimer[0] in absFrequencyOut:     #add the relative conditional probability in a list
                  value = round(dimers[dimer]/absFrequencyOut[dimer[0]],2)
                  OutsideModel+= [value,]

      i = 0
      for h in dimers:      #use the value of the probabilities obtained before and
            OutsideMatrix[h] =OutsideModel[i]    #the dimers to create the outside "matrix"
            i+=1

      return OutsideMatrix  #return the outside matrix seen as a python dictionary {dimer:probability}

def models():
      
      inside = InsideModel()
      outside = OutsideModel()
      
      return inside,outside


      
function = models()
print('The inside model is:',function[0])
print('The outside model is:',function[1])

##NB: if the user wants to generate the models more the once, each time it has to initialize
##      6 GLOBAL VARIABLES that are initialized at the beginning of the program
      
            



