def ReverseComplement(Pattern):
    """Returns reverse complement of Pattern"""
    comp_bps = {"A":"T", "T":"A", "G":"C", "C":"G"}
    comp_seq = []
    
    for bp in Pattern:
        comp_seq.insert(0,comp_bps[bp])
    
    return ''.join(comp_seq)

def PatternToNumber(Pattern):
    """Returns a numeric code for Pattern according to base pair sequence in Pattern"""
    sym_dict = {"A":0, "C":1, "G":2, "T":3}
    symbol = Pattern[-1]
    
    if len(Pattern) == 1:
        return sym_dict[symbol]
    
    else:
        prefix = Pattern[:len(Pattern)-1]
        return 4 * PatternToNumber(prefix) + sym_dict[symbol]

def NumberToPattern(number, k):
    """Returns the base pair sequence of length k based on number"""
    from itertools import product
    
    bps = ["A", "C", "G", "T"]
    mer_dict = {}
    
    for index, perm in enumerate(itertools.product(bps, repeat = k)):
        mer_dict["".join(perm)] = index

    return [key for key, value in mer_dict.items() if value == number][0]    

def Neighbors(Pattern, d):
    """Generates a list of all kmers within a Hamming Distance <= d of Pattern"""
    nts = ["A", "C", "G", "T"]
    
    #Exact match only case
    if d is 0:
        return [Pattern] #Square braces ensures a list is returned, otherwise string is returned
    
    #Recursion base case: last bp in pattern, returns all nts
    if len(Pattern) is 1:
        return nts
    
    #neighborhood_list stores the neighbors with Hamming Distance <= d
    neighborhood_list = []
    
    #Create suffix of Pattern
    suffix = Pattern[1:]
    
    #Recursively generate all suffixes to build neighborhood_list in next steps
    suffix_neighbors = Neighbors(Pattern[1:], d)
    
    #Generate neighbors according to Hamming Distance logic
    for text in suffix_neighbors:
        
        #If Hamming Distance < d, we can add any nt to the suffix
        if HammingDistance(Pattern[1:], text) < d:
            for nt in nts:
                neighborhood_list.append(nt + text)
            
        #Otherwise, Hamming Distance = d, so we can only add back symbols from the pattern
        else:
            neighborhood_list.append(Pattern[0] + text)
        
    return neighborhood_list
    
def FrequentWordsWithMismatchesAndReverseComplementsSorting(Text, k, d):
    """Finds most frequent k-mers and their reverse complements in Text with d mismatches or less
    Implements a frequency array and neighbors sorting algorithm to achieve significant reduction in run time relative to nested loops"""
    
    import time
    
    FrequencyDictionary = {}
    FrequentPatterns = []
    Neighborhoods = []
    rc_Neighborhoods = []
    Index = []
    
    start = time.time()
    #Populate kmer_index_list with kmers in Text
    for index in range(len(Text)-k+1):
        kmer = Text[index:index+k]
        Neighborhoods.append(Neighbors(kmer, d))
    
    #Generate reverse complements of Neighborhoods and store in temporary list
    for i in range(len(Neighborhoods)):
        rc = [ReverseComplement(j) for j in Neighborhoods[i]]
        rc_Neighborhoods.append(rc)
    
    Neighborhoods.extend(rc_Neighborhoods)
    
    #return Neighborhoods
    finish = time.time()
    print("Creating Neighborhoods took", finish-start)
    
    start = time.time()
    #Converts contents of Neighborhoods to Number Codes and appends to Index
    for i in range(len(Neighborhoods)):
        ptn = [PatternToNumber(j) for j in Neighborhoods[i]]
        Index.append(ptn)
    
    #Flattens Index from a 2-D list to a 1-D list
    Index = [number for sublist in Index for number in sublist]
    
    #Sorts index for counting in next section
    SortedIndex = sorted(Index)
    finish = time.time()
    print("Creating SortedIndex took", finish-start)
    
    start = time.time()
    #Counts the number of identical, adjacent numbers. This is equal to frequency of the underlying pattern.
    for i in range(len(SortedIndex)-1):
        if SortedIndex[i] == SortedIndex[i+1]:
            FrequencyDictionary[SortedIndex[i]] = FrequencyDictionary.get(SortedIndex[i], 0) + 1
    
    #Finds maximum value from previous step by sorting FrequencyDictionary according to its values, and taking the last value
    max_freq = sorted([value for value in FrequencyDictionary.values()])
    max_freq = max_freq[-1]
    finish = time.time()
    print("Finding max_freq took", finish-start)
    
    start = time.time()
    FrequentNumbers = [key for key, value in FrequencyDictionary.items() if value == max_freq]
    
    FrequentPatterns = [NumberToPattern(number,k) for number in FrequentNumbers]
    finish = time.time()
    print("Finding Frequent Numbers and Patterns took", finish-start)
    
    return FrequentPatterns
