"""
READING WEIGHTS FILE: reading the weights file line by line then add it's content in two different dictionaries,
'dict' in which the weight represents the key while the amino acid is the value 
and 'reverse_dic' the key here is the amino acid while the weight represents it's value.
"""
file = open("D:\\3rd year\\Second term\\Structural Bioinf\\task 1\\weight.txt")
dict={}
reverse_dic = {}
for f in file:
    line= f.rstrip()
    words=line.split(" ")
    aacid = words[0]
    weight = int(words[1])
    dict[weight]= aacid
    reverse_dic[aacid] = weight


"""
this should be the first function to be executed, it takes the input spectrum as it's input and determine from it all the amino acids composing the cyclic protien sequence.
as the largest amino acid in weights file is of weight 186 and all the k-mers  (where k>1) will exceed this values we can use this as a condition to detect the amino acid composing the sequence.
"""
def initial_list(spectrum):
    init_list = []
    for s in spectrum:
        if s< 186 and s!=0:
           init_list.append(dict[s])
    return init_list

"""
this function takes a sub-peptide sequence and calculate it's linear spectrum but without the circulation condition
"""

def linear_spectrum(sub_peptide):

    sublist = []
    for i in range(len(sub_peptide)):
        for j in range(len(sub_peptide)):
            if (i + j) < len(sub_peptide):
                sublist.append(sub_peptide[j:j + i + 1])        #generating all the possible sub_peptides composing the sub_peptide taken as a parameter.

    totals = [0]
    for i in sublist:
        totals.append(weight(i))                                #then calculate the weight for each of those sub_peptide to construct the linear spectrum.
    return sorted(totals)

"""
this function take a sub_peptide and calculate it's total weight by using the 'reverse_dic' dictionary.
#for each amino acid found in the sub_peptide, add it's weight to the final weight.
"""
def weight(sub_peptide):
    total=0
    for c in sub_peptide :
        total+= reverse_dic[c]
    return total

"""
This function check if a sub_peptide is consistent through:
1) calculate the linear spectrum for the sub_peptide taken as a parameter(which corresponds to a possible combination) by calling 'linear_spectrum' function.
2) check that each and every value in the returned linear spectrum list is exist in the input spectrum.
3) number of occurrences for every value in the returned linear spectrum can't exceed it's number of occurrences in the input spectrum.
in order to achieve that take a copy of the input spectrum and when a value in the returned linear spectrum is found in that copy change it's value to -1.
"""
def IsConsistent(sub_peptide, input_spectrum):

    sub_peptideSpectrum = linear_spectrum(sub_peptide)
    consistent = True
    input_spectrum_temp = input_spectrum.copy()
    for value in sub_peptideSpectrum:
        found = False
        for i in range(len(input_spectrum_temp)):
            if value == input_spectrum_temp[i]:
                found= True
                input_spectrum_temp[i] = -1
                break
        if found == False:                                   #once there's a weight that it's not in the input spectrum set 'consistent' flag to false.
            consistent = False
            break
    return consistent

#This function extend each amino acid and create all the possible combnations, if the created combination is consistent add it to a list of items thatwill be extended in the next iteration.

def get_Combinations(init_list ,spectrum):

    temp_list = init_list
    last_combinations = []
    for i in range(1,len(temp_list)):                        #this loop iterate k-1 times where k is the number of amino acid composing the protien sequence, each iteration we generate the consistent i+1-mers, why k-1? because we already have the 1-mers which are in the 'init_list'.
        for aa in temp_list:                                 #for each consistent i-mers generated in the previous iteration.
            for j in init_list:                              #concat it with all the amino acids
                temp = aa + j
                consistent = IsConsistent(temp, spectrum)
                if consistent == True:
                    if (temp in last_combinations) == False:
                        last_combinations.append(temp)      #if this combination is consistent and unique(not generated before) add it to the list that contains the current i-1-mers
        temp_list = last_combinations.copy()
        print('List of', i+1, '-mers :', temp_list)
        last_combinations.clear()
    return temp_list

#---the main function---

def cyclopeptide_sequencing(spectrum):
    init_list = initial_list(spectrum)
    print('List of 1-mers :', init_list)
    result = []
    result = get_Combinations(init_list, spectrum)
    return result

#main------------------------
spectrum = [0,97,97,99,101,103,196,198,198,200,202,295,297,299,299,301,394,396,398,400,400,497]
result = cyclopeptide_sequencing(spectrum)
print('All the final linear representation of the peptide sequence : ',result)