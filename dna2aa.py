# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 7:41:44 2015

@author: Mikalai Drabovich
"""

import unittest
import sys 
 
class dna2aa(object):

  rna_codon_table = \
       {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
        "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
        "UAU":"Y", "UAC":"Y", "UAA":".", "UAG":".",
        "UGU":"C", "UGC":"C", "UGA":".", "UGG":"W",
        "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
        "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
        "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
        "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
        "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
        "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
        "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
        "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
        "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
  def convert(self,rna):
    
    subseqs = []
    result = [] 
    
    startCodon = 'AUG'
    stopCodon0 = 'UAA'
    stopCodon1 = 'UAG'
    stopCodon2 = 'UGA'
    
    start_index = 0
    total_chars = len(rna)
  
    while start_index<total_chars-2:
      if rna[start_index:start_index+3] == startCodon:
         triplets = [rna[i:i+3] for i in range(start_index, len(rna), 3)]
         k=0
         while k<len(triplets):
           if triplets[k]==stopCodon0 or triplets[k]==stopCodon1 or triplets[k]==stopCodon2:
             break
           else:
             k=k+1
         
         temp_str = ''
         
         if k>=len(triplets): 
           # try new partitioning into triplets, startting from start_index + 1
           start_index = start_index + 1
           continue
       
         # case when string between start and stop contains noise:
         hasNoise = False
         numStarts = 0
         for triplet in triplets[0:k+1]:
           try:
               if triplet == startCodon:
                  numStarts = numStarts + 1
               if numStarts>1:
                  hasNoise = True 
                  break              
               temp_str = temp_str + self.rna_codon_table[triplet]
           except KeyError:
               hasNoise = True 
               break

         if hasNoise:
           start_index = start_index + 1
           continue
           
         # so, all the triplets has been successfully translated!
         result.append(temp_str)
         start_index = start_index + 3*(k+1)
         subseqs.append(''.join(triplets[0:k+1]))
         
      else:
         start_index=start_index+1
          
    return (subseqs, result) if len(result)>0 else (None, None)
    
  def dna2rna(self, dna):
    return ''.join(c if c!='T' else 'U' for c in dna)
 

  def rna2dna(self, rna):
    return ''.join(c if c!='U' else 'T' for c in rna)
        
class dna2aaTest(unittest.TestCase):
  """Unit tests for some typical and some edge cases """
  
  def test_typical_example(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGACGGCUUGUUUCUUUUCUGUGGCUGCGUGA"),['MTACFFSVAA.'])

  def test_typical_example_doubled(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGACGGCUUGUUUCUUUUCUGUGGCUGCGUGAAUGACGGCUUGUUUCUUUUCUGUGGCUGCGUGA"),['MTACFFSVAA.','MTACFFSVAA.'])
        
  def test_typical_example_doubled_with_noise(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("---SOME-NOISE---AUGACGGCUUGUUUCUUUUCUGUGGCUGCGUGA---SOME-NOISE---AUGACGGCUUGUUUCUUUUCUGUGGCUGCGUGA---SOME-NOISE---"),['MTACFFSVAA.','MTACFFSVAA.'])
  
  def test_typical_example_doubled_pseudo_stops(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("---PSEUDO-STOPS:UAAUGAUAG---AUGACGGCUUGUUUCUUUUCUGUGGCUGCGUGA---PSEUDO-STOPS:UAAUGAUAG---AUGACGGCUUGUUUCUUUUCUGUGGCUGCGUGA---PSEUDO-STOPS:UAAUGAUAG---"),['MTACFFSVAA.','MTACFFSVAA.'])
    
  def test_typical_example_doubled_pseudo_starts(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("---PSEUDO-STARTS:AUGAUGAUG---AUGACGGCUUGUUUCUUUUCUGUGGCUGCGUGA---PSEUDO-STARTS:AUGAUGAUG---AUGACGGCUUGUUUCUUUUCUGUGGCUGCGUGA---PSEUDO-STARTS:AUGAUGAUG---"),['MTACFFSVAA.','MTACFFSVAA.'])

  def test_choosing_the_nearest_stop1(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGUAAUAA"),['M.'])
    
  def test_choosing_the_nearest_stop2(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGAUGUAAUAA"),['M.']) 
 
  def test_noise(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("A"), (None, None))
    self.assertEqual(converter.convert("AA"),(None, None))
    self.assertEqual(converter.convert("AAAAAA"),(None, None))
    self.assertEqual(converter.convert("AAAAAAAAA"),(None, None))
    
  def test_only_start_stop(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGUAA"),['M.'])
    self.assertEqual(converter.convert("AUGUGA"),['M.'])
    self.assertEqual(converter.convert("AUGUAG"),['M.'])

  def test_start_noise_stopUAA(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGAUAA"),(None, None))
    self.assertEqual(converter.convert("AUGAAUAA"),(None, None))
    self.assertEqual(converter.convert("AUGABCUAA"),(None, None))
    self.assertEqual(converter.convert("AUGUUUAUAA"),(None, None))
    self.assertEqual(converter.convert("AUGXXXUAA"),(None, None))
  
  def test_start_noise_stopUGA(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGAUGA"),(None, None))
    self.assertEqual(converter.convert("AUGAAUGA"),(None, None))
    self.assertEqual(converter.convert("AUGABCUGA"),(None, None))
    self.assertEqual(converter.convert("AUGUUUAUGA"),(None, None))
    self.assertEqual(converter.convert("AUGXXXUGA"),(None, None))
     
  def test_start_noise_stopUAG(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGAUAG"),(None, None))
    self.assertEqual(converter.convert("AUGAAUAG"),(None, None))
    self.assertEqual(converter.convert("AUGABCUAG"),(None, None))
    self.assertEqual(converter.convert("AUGUUUAUAG"),(None, None))
    self.assertEqual(converter.convert("AUGXXXUAG"),(None, None))
     
  # --- checking all codon table ... ---
  def test_start_UUU_stop(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGUUUUAA"),['MF.'])
    self.assertEqual(converter.convert("AUGUUUUGA"),['MF.'])
    self.assertEqual(converter.convert("AUGUUUUAG"),['MF.'])

  def test_start_UUC_stop(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGUUUUAA"),['MF.'])
    self.assertEqual(converter.convert("AUGUUUUGA"),['MF.'])
    self.assertEqual(converter.convert("AUGUUUUAG"),['MF.'])
    
    '''
    ... here we would test all codon map - 60 + 4 entries ...
    '''
  def test_start_GGG_stop(self):
    converter = dna2aa()
    self.assertEqual(converter.convert("AUGUUUUAA"),['MF.'])
    self.assertEqual(converter.convert("AUGUUUUGA"),['MF.'])
    self.assertEqual(converter.convert("AUGUUUUAG"),['MF.'])
    
  # --- end of checking codon table ... ---  
    

if __name__ == '__main__':        
    
  #unittest.main()  

  if len(sys.argv)<3:
      print "Usage:\n --rna dna.txt : converts DNA to RNA; \n --genes dna.txt : finds sequence of amino acids which can be produced from the input DNA sequence"

  else:   
      with open(sys.argv[2], 'r') as f:
          contents = f.read()
          
      converter = dna2aa()   
      
      if sys.argv[1]=='--rna':
        print converter.dna2rna(contents)
      elif sys.argv[1]=='--genes':
        result = converter.convert(converter.dna2rna(contents))
        subseq = result[0]
        aminoacids = result[1]
        print "%d protein(s) can be produced from this sequence:\n" % (len(result[0]))
        for k in range(0,len(subseq)): 
          print "DNA: %s" % (converter.rna2dna(subseq[k]))
          print "RNA: %s" % (subseq[k])
          print "Amino acid sequence: %s\n" % (aminoacids[k])
          