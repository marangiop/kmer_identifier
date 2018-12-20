# -*- coding: utf-8 -*-
"""
Created on Sun May  6 23:24:28 2018
@author: maran
"""
from __future__ import division

#STEP 1: Extraction of DNA sequences from the input file

while True: #This while loop will continute infinetely unless stopped when a condition is met
    Actual_file_name = input("Enter file name including .fastq extension : \n ") #Allows user to enter the name of the file
    if Actual_file_name.endswith('.fastq'):  # Checks whether the file name entered by the user contains the string  'fastq' at the end
        break #If the above condition is fulfilled, the while loop will be exited
    else:
        print("Error! Please enter an input file in .fastq format!") #Prints this error statement and re-starts the loop.

def DNA_sequences_extractor(file_name):
    my_file = open(file_name)                            
    extracted = []                         
    for line in my_file: #This for loop treats the file as a list of lines and interates over each line                       
        extracted.append(line.rstrip("\n"))
    return extracted[1::4]             #This function returns every 4th line starting from element with index 1 from the list 'extracted'
extracted_DNA_sequences=DNA_sequences_extractor(file_name=Actual_file_name) 
#_____________________________________________________________________________________________________________________________

#STEP 2: Generation of a list of kmers for each sequence

while True: #This while loop will continute infinetely unless stopped when a condition is met
    try:
        k_value = int(input("Enter value of kmer length: \n ")) #Allows user to enter the kmer length and checks whether it is an integer
        break #If the above condition is fulfilled, the while loop will be exited
    except: #If there is an exception 
        print("That is not a valid number. Please insert an integer") #Prints this error statement and re-starts the loop.

def kmers_generator(dna_sequence_list, kmer_length):
    kmers = []
    for dna_sequence in dna_sequence_list:#This for loop iterates over a list of DNA sequences provided
        for start in range(len(dna_sequence) - kmer_length + 1): #Thi for loop iterates over every starting position in the DNA sequence. Subtracting k and adding 1 allows to stop the range() so that incomplete kmers at the end of the sequence are not generated.
            stop = start + kmer_length #This calculates the upper limit of the slice taken from the DNA sequence to generate the kmers.
            kmer = dna_sequence[start:stop] #This generates a kmer of desired kmer length at that position in the DNA sequence
            if kmer not in kmers: #This allows to append the kmers that are not already present in the list kmers.
                kmers.append(kmer) 
    return kmers #This will generate a list of non-redundant kmers for a given list of sequences 

all_kmers= kmers_generator(dna_sequence_list=extracted_DNA_sequences, kmer_length=k_value)
#_____________________________________________________________________________________________________________________________

#STEP 3:Generation of a collection of counts at each position across all the sequences for each kmer 

def count_kmer_at_positions_across_sequences(dna_list, kmer_length, kmer_1):
    kmers= {}
    for dna_sequence in dna_list: 
        for position in range(len(dna_sequence)-kmer_length+1): 
            kmer= dna_sequence[position:position+kmer_length]
            if kmer == kmer_1:   #This checks whether the kmer being generated at that position (kmer) is equal to any kmer that is fed into the function (kmer_1)
                kmers[position]= kmers.get(position, 0) +1  #This looks up the current count for that position and adds 1. The position acts as a key in the dict and the count is a stored as an integer value assigned to each position
    for position in range(len(dna_sequence) - kmer_length + 1): 
        if position not in kmers: #If a position has not been assigned any count and is therefore not present in the dict kmers, which holds the count at each position
            kmers[position] = 0 #The default value of such entry is set to 0
    return kmers #This returns a dict holding the count of a given kmer at all positions across a list of DNA sequences.

dict_of_dicts_counts_kmers={}
for kmer in all_kmers:   
    dict_of_dicts_counts_kmers[kmer]= count_kmer_at_positions_across_sequences(dna_list=extracted_DNA_sequences, kmer_length=k_value, kmer_1= kmer) #Each dict produced from feeding into the function each given kmer from the list of non-redundant kmers "all_kmers" is assigned to a new key, which corresponds to each given kmer. In this way a dict of dicts is generated.
#_____________________________________________________________________________________________________________________________

#STEP 4:Calculation of the expected count for each kmer

def kmer_expected_counts(dna_list, kmer_length):
    sequence_length_minus_kmer_length= len(dna_list[0])-kmer_length + 1 #This corresponds to the number of possible starting positions in a given sequence. In this case it is calculated based on the first sequence of the list of extracted sequence, but any other sequence could be used as it is assumed that all sequences have the same length
    f = {}
    for dna_sequence in dna_list:
        for position in range(len(dna_sequence)-kmer_length + 1):
            kmer = dna_sequence[position:position+kmer_length]
            f[kmer] = f.get(kmer, 0) + 1  #This looks up the current count for that kmer and adds 1.
    for kmer, count in f.items(): #This for loop iterates over every kmer and count in the dict f
        f[kmer]= count/sequence_length_minus_kmer_length # Each kmer is now assigned to the total count of that kmer devided by sequence_length_minus_kmer_length
    return f #This returns a dict that holds the expected count for each kmer across all the sequences
dict_expected_counts = kmer_expected_counts(dna_list = extracted_DNA_sequences, kmer_length= k_value) 
#_____________________________________________________________________________________________________________________________

#STEP 5: Comparing the expected count and the observed count for each kmer at each position

while True: #This while loop will continute infinetely unless stopped when a condition is met
    try:
        threshold_value= int(input("Enter value of reporting threshold: \n ")) #Allows user to enter the threshold value
        break #If the above condition is fulfilled, the while loop will be exited
    except: #If there is an exception 
        print ("That is not a valid number. Please insert an integer") #Prints this error statement and re-starts the loop.
        
def kmer_flagging(kmer, threshold, expected):
    for kmer1, dictionary in dict_of_dicts_counts_kmers.items():
        if kmer1 == kmer:     #This checks whether the kmer found in the list of observed counts (kmer) is equal to any kmer that is fed into the function (kmer_1)   
            for position, observed_count in dictionary.items(): #This for loop iterates over every position and observed count in each dictionary
                if (observed_count/expected) > threshold: #The following command is performed only for those kmers whose observed count devided by its respective expected count is greater than a value that is assigned by the user to the variable 'threshold' through raw_input function
                   return kmer, int(expected), dictionary.values() #The function returns the kmer sequence, the expected count value rounded up to a whole number and the counts for that kmer at each position

#_____________________________________________________________________________________________________________________________
     
#STEP 6: Generation of an output file that stores information about ‘suspicious’ kmers 

my_file = open(Actual_file_name.rstrip('.fastq') + "_" + "kmer_length" + "_" + str(k_value) + "_" + "reporting_threshold" + "_" + str(threshold_value) + '.txt', "w") # Automatically defines a file name that includes the input file name stripped of the .fastq extension and the kmer legnth and reporting threshold values. The '.txt' string is automatically added at the end of each output file name. This enables it to be opened with a word processor like Notepad or Microsoft Word
for kmer0, expected_count in dict_expected_counts.items(): #This for loop iterates over every kmer and expected count in the dict holding the expected counts 
    if kmer_flagging(kmer=kmer0,threshold=threshold_value, expected=expected_count) != None: #This checks whether a given kmer does not return None after it has been fed into the kmer_flagging() function. 
        output = str(kmer_flagging(kmer=kmer0,threshold=threshold_value, expected=expected_count)) #Every kmer is fed into function kmer_flagging and the returned information is stored as a string in the variable output
        my_file.write('%s \n' %output) #Writes the information stored in variable output to the output file and adds a new line in the file
my_file.close()
