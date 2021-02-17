import random
import math

inside = {'AA': 0.19, 'AC': 0.28, 'AG': 0.4, 'AT': 0.14, 'CA': 0.19, 'CC': 0.36, 'CG': 0.25, 'CT': 0.2, 'GA': 0.17, 'GC': 0.33, 'GG': 0.36, 'GT': 0.14, 'TA': 0.09, 'TC': 0.34, 'TG': 0.38, 'TT': 0.19}
outside = {'AA': 0.29, 'AC': 0.2, 'AG': 0.29, 'AT': 0.22, 'CA': 0.32, 'CC': 0.3, 'CG': 0.07, 'CT': 0.31, 'GA': 0.25, 'GC': 0.24, 'GG': 0.3, 'GT': 0.21, 'TA': 0.18, 'TC': 0.23, 'TG': 0.3, 'TT': 0.29}

def checkModel(P):   #this function calculate the score for a unique sequence

      P = P.upper()         #makes everything upper

      global inside       #takes the inside and the outside model for CpG islands
      global outside     #calculated by the program models
      
      result = 0       #here the final score will be collected 

      d = 2
      n = 0
      probInside = math.log(0.25)    #initial probability for the inside score 
      probOutside = math.log(0.25)    #initial probability for the outside score 

      while d <= len(P)-1:                             #calculate the logarithm of the probability of each dimer
            ins = math.log(inside[P[n:d]])    #conditional probability of being inside a CpG
            probInside = probInside + ins
            out = math.log(outside[P[n:d]])   #conditional probability of being outside a CpG
            probOutside = probOutside + out
            n +=1       #increment these values to move into the next following dimer
            d+=1

      result = probInside - probOutside     #final score calculation

      
      if result > 0:
            return True,result      #TRUE means that the sequence is a CpG 
      else:
            return False,result    #FALSE means that the sequence is not a CpG

def slidingWindow(L,AvLength):   #this function calculates the score of streches of a sequence
                                                          #obtained using a sliding window of size AvLength
      L = L.upper()
      resultWindow = []   #list that will contain tuples having TRUE or FALSE and the score associated to each strech
      slides = []
      offsets = []

      for e in range(len(L)-AvLength):    #slide the sequence L with a window of size AvLength
            slides += [L[e:e+AvLength],]
            offsets += [int(e),]
            window = checkModel(L[e:e+AvLength])   #calculate the score of the current window 
            resultWindow += [window,]             #collect the result

      return resultWindow,slides,offsets    #the final output is a list conting the scores of each window 

def test():
      
      select = int(input('Insert 0 to generate a random sequence or 1 if you want to insert a sequence:'))
      if select == 0:
            length = int(input('Insert the length of the random sequence:'))
            genome = "".join(random.choice('AGCT') for i in range(length))
      if select == 1:
            genome = input('Insert the genome sequence for which you want to test the CpG island model:')
      choice = int(input('Insert 1 for single window check; insert 2 for a sliding window check:'))
      if choice == 1:
            program = checkModel(genome)
      if choice == 2:
            window = input('Insert a number of the window or press ENTER for the default value (567):')
            if window != '':
                   window = int(window)
                   program = slidingWindow(genome,window)
            else:
                  window = 567
                  program = slidingWindow(genome,567)
      print('Results')
      if choice == 1:
            print(genome.upper())
            if program[0] == True:
                  print('The sequence under analysis is a CpG island')
                  print('The score of the model evaluation is:',program[1])
            else:
                  print('The sequence under analysis is not a CpG island')
                  print('The score of the model evaluation is:',program[1])

      count = 0
      if choice == 2:
            output = program[0]
            i= 0
            for r in output:
                  if r[0] == True:
                        print(program[1][i].upper())
                        print('The sequence under analysis is a CpG island')
                        print('The score of the model evaluation is:',r[1])
                        print('The range of the CpG found in the genome sequence is:',[program[2][i],program[2][i]+int(window)+1])
                        i+=1
                        count += 1
                  else:
                        print(program[1][i].upper())
                        print('The sequence under analysis is not a CpG island')
                        print('The score of the model evaluation is:',r[1])
                        i+=1

            print('The number of CpG found is:',count)
                  
      return 'End of the program'

print(test())                 
      
