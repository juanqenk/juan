---
layout: post
title:  "Bioinfo"
date:   2023-01-15
desc: "Bioinfo"
keywords: "markdown,Jekyll,gh-pages,website,blog"
categories: [Linux]
tags: [Markdown,Github]
icon: fab fa-markdown
---


# A Rapid Introduction to Molecular Biology

Making up all living material, the cell is considered to be the building block of life. The nucleus, a component of most eukaryotic cells, was identified as the hub of cellular activity 150 years ago. Viewed under a light microscope, the nucleus appears only as a darker region of the cell, but as we increase magnification, we find that the nucleus is densely filled with a stew of macromolecules called chromatin. During mitosis (eukaryotic cell division), most of the chromatin condenses into long, thin strings called chromosomes.

One class of the macromolecules contained in chromatin are called nucleic acids. Early 20th century research into the chemical identity of nucleic acids culminated with the conclusion that nucleic acids are polymers, or repeating chains of smaller, similarly structured molecules known as monomers. Because of their tendency to be long and thin, nucleic acid polymers are commonly called strands.

The nucleic acid monomer is called a nucleotide and is used as a unit of strand length (abbreviated to nt). Each nucleotide is formed of three parts: a sugar molecule, a negatively charged ion called a phosphate, and a compound called a nucleobase ("base" for short). Polymerization is achieved as the sugar of one nucleotide bonds to the phosphate of the next nucleotide in the chain, which forms a sugar-phosphate backbone for the nucleic acid strand. A key point is that the nucleotides of a specific type of nucleic acid always contain the same sugar and phosphate molecules, and they differ only in their choice of base. Thus, one strand of a nucleic acid can be differentiated from another based solely on the order of its bases; this ordering of bases defines a nucleic acid's primary structure.

For reasons we will soon see, DNA is found in all living organisms on Earth, including bacteria; it is even found in many viruses (which are often considered to be nonliving). Because of its importance, we reserve the term genome to refer to the sum total of the DNA contained in an organism's chromosomes.

```{admonition} Assignment: **Counting DNA Nucleotides**
:class: important

A string is simply an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains. An example of a length 21 DNA string (whose alphabet contains the symbols ‘A’, ‘C’, ‘G’, and ‘T’) is “ATGCTTCAGAAAGGTCTTACG.”

**Given**: A DNA string s of length at most 1000 nt.

**Return**: Four integers (separated by spaces) counting the respective number of times that the symbols ‘A’, ‘C’, ‘G’, and ‘T’ occur in s.

    Sample Dataset

<code>AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC</code>

    Sample Output

<code>20 12 17 21</code>

```

We need an algorithm that counts the number of occurrences of the symbols 'A', 'C', 'G', and 'T' in a given string, represented by the variable "s". The algorithm starts by creating a dictionary with keys 'A', 'C', 'G', 'T' and initializing their values to 0. Then, the algorithm iterates through each symbol in the input string "s". For each symbol encountered, the algorithm checks if the symbol exists as a key in the dictionary, and if so, increments the value for that key by 1. After the iteration is complete, the algorithm returns a string containing the values for the keys 'A', 'C', 'G', and 'T' in the "counts" dictionary.


```{admonition} Algorithm: **Counting DNA Nucleotides**
:class: danger

_________________________________________________________________________________


<code>**create a dictionary** "counts" with keys "A", "C", "G", "T" and set all of them to 0     
**for each** symbol in the string "s"
    **if** symbol **exists** in the dictionary "counts"
        **increment the count for the key** "symbol" in the "counts" dictionary    
**return** a string containing the values for the keys "A", "C", "G", "T" in the "counts" dictionary</code>
_________________________________________________________________________________
```




```python
def count_symbols(s):
    """
    count_symbols
    ==============

    .. function:: count_symbols(s: str) -> str

    This function takes a DNA sequence string as input, and returns a string containing the count of each symbol ('A', 'C', 'G', and 'T') in the input string.

    :param s: The DNA sequence string.
    :returns: A string containing the count of each symbol in the input string, in the format "count_A count_C count_G count_T".

    Example:
        >>> count_symbols("AGCTAGCT")
        "2 2 2 2"

    """
    # Create a dictionary to store the count of each symbol
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    # Iterate through the string and increment the count for each symbol
    for symbol in s:
        if symbol in counts:
            counts[symbol] += 1

    # Return the count of each symbol as a string
    return f"{counts['A']} {counts['C']} {counts['G']} {counts['T']}"
```


```python
s = 'AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'
count_symbols(s)
```




    '20 12 17 21'




```python
s = 'AGGTAATTTTACTTCCCTTTGTCCAATGATGCTCATAAGGCCTCGGCGACTCAAAGTGAATAACTGCACGGTGAGTGGTCTCACCTTCAAATTATCTAGAGCTAGCTAGGTTACATAAGTATCCTTAGGGTTCGGTTGTATGAATGGGTGAAGATTTACATAGGCCATCGTGATGAAGAGTTACGCAGTAAGAAATTGGCACTAGCTTGCAACTGATCTAACGTGTAGAATCACATACACGTTCGTTGTAGCTCATCAATCCACTCGCCGTCGAGCGAGTGGGGTGTTCGCCAAGTTCGCTGGGTCCTTATCCGATGTACTTACACTAGTATTTACGGCTTATACGACCAATTAACACTAATGCAGATTAAGTAGGCCGTCTCCGGGGGAATTTCTCTGTGTCCTGCTCGGCTGACATATGTCACGCGCCCTTCTGTTGTGGTGCTAACCTAGTGAATTGGAGATGAATTATTGGCATTCGGCTGTGTCGCATCGTGACGGGCTGTATATCCACACCGCGGCTGAAATGTGCTGTAACTCCTATTAGCCAGTGCTCCAGGGGGCAGTGATGACCAACGTTGCCAGCGGTATCAATCCGTTCGTACTTATTTTTTAAGGAGCAAATGGGAACAATGTTTATTGTCGGGGTCTTGATCCGTCGGGTTCATTCTACGCTCGCTTAAGTATACATCCACTCCCAAGCCTAGTGCCAGACAACTCCGGCCCGATATAAGCTCGTGAGGACCCTTGTATTATCCATGTAGTTTTATCGGTGCGTTGAATATATGTTTAGAACAATTAAGAATGCTCGGTTGGGGGACGTAAAATCGGTAACGCCACTGG'
count_symbols(s)
```




    '204 186 208 245'




```python
help(count_symbols)
```

    Help on function count_symbols in module __main__:
    
    count_symbols(s)
            count_symbols
        ==============
        
        .. function:: count_symbols(s: str) -> str
        
        This function takes a DNA sequence string as input, and returns a string containing the count of each symbol ('A', 'C', 'G', and 'T') in the input string.
        
        :param s: The DNA sequence string.
        :returns: A string containing the count of each symbol in the input string, in the format "count_A count_C count_G count_T".
        
        Example:
            >>> count_symbols("AGCTAGCT")
            "2 2 2 2"
    



```python
def test_count_symbols():
    assert count_symbols("AGCTAGCT") == "2 2 2 2"
    assert count_symbols("AGCTAGCTAGCT") == "3 3 3 3"
    assert count_symbols("ACGTACGT") == "2 2 2 2"
    assert count_symbols("") == "0 0 0 0"
    assert count_symbols("ATCGATCG") == "2 2 2 2"
test_count_symbols()
```

The time complexity of this algorithm is O(n), where n is the number of symbols in the input string. The algorithm needs to iterate through the entire input string once, so the time it takes to run the algorithm grows linearly with the size of the input. Space complexity is also O(1) since it creates a constant number of fixed size dictionary and doesn't use any data structures that grow with input size.

It's a simple and easy to understand solution, but it also has a very low time complexity and space complexity, which makes it efficient and fast.

# The Second Nucleic Acid

A second nucleic acid exists alongside DNA in the chromatin; this molecule, which possesses a different sugar called ribose, came to be known as ribose nucleic acid, or RNA. RNA differs further from DNA in that it contains a base called uracil in place of thymine. Biologists initially believed that RNA was only contained in plant cells, whereas DNA was restricted to animal cells. However, this hypothesis dissipated as improved chemical methods discovered both nucleic acids in the cells of all life forms on Earth.

The primary structure of DNA and RNA is so similar because the former serves as a blueprint for the creation of a special kind of RNA molecule called messenger RNA, or mRNA. mRNA is created during RNA transcription, during which a strand of DNA is used as a template for constructing a strand of RNA by copying nucleotides one at a time, where uracil is used in place of thymine.

In eukaryotes, DNA remains in the nucleus, while RNA can enter the far reaches of the cell to carry out DNA's instructions. 

**Assignment: Transcribing DNA into RNA**

An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.

Given a DNA string t corresponding to a coding strand, its transcribed RNA string u is formed by replacing all occurrences of 'T' in t with 'U' in u.

Given: A DNA string t having length at most 1000 nt.

Return: The transcribed RNA string of t.

    Sample Dataset

<code>GATGGAACTTGACTACGTAAATT</code>

    Sample Output

<code>GAUGGAACUUGACUACGUAAAUU</code>

This code must be a function that takes in a string "t" which is a DNA sequence and replace the character T with U character, this is because in DNA the T is replaced with U in the RNA.

<code>**Replace** all occurrences of 'T' with 'U'
**Store** result of the replacing in a variable "u"
**return** "u" </code>


```python
def transcribe_rna(t):
    # Replace all occurrences of 'T' with 'U'
    u = t.replace('T', 'U')

    # Return the transcribed RNA string
    return u
```


```python
t = 'GATGGAACTTGACTACGTAAATT'
transcribe_rna(t)
```




    'GAUGGAACUUGACUACGUAAAUU'




```python
t = 'GATGTTTCCACTTATCCAGTTGGCTGGATGTGTTAGGCGCCCCTGATAATGAGCATCAACTGATCTGTGATGTAGTGGAGTAGCACTTTTCGATAAACTGTTTGAGCAGCCTGGACAGCTCGACGCAGACCCACATGCACACGGTATTTCTCCCGCGACCTGTCAGACTAACGAGTGCTAAAGCTTTAAGATGAGATAAAGTTGTCCTCTTGGGGACAAATAATGCGTCGAAATTGACGTCTATTACCCCTGAGCACCATTATAGAGCATTACCAGAGGAATGATACTCCGGGAAGCTATATAATTTCGAATGCGTCTGCCCTAGATACGGAATGAGTTTAGATGCGCTAGAGTGATGCCGATCTTTCCATGAAAGTGGGTGTAACTAGTAAAGGTGTGACTGGGCCTACGTCCTCGTATCGTGGACCTCCACACCATGACTGCGGCGCTCTTTGGCACCGCGCGACCCTCTTATACCTAGCTCGCGGGGTCGCTGGGGAAACAGCCCCTATGGGGGCAAGCGACTCAGTAAATCTCTGTTTACGTACCTAGTTATTATGACCAGATGCATTCTGTCATCGATCTAATCGCAGACAATTCTCCTCTCGACGGGGTCGTAAGTCCAGCGACTGGGGCAGACCGCGCTCAGTATACATCTAAGGGTATGGTGCGCTGTATAATACCGCCATAAAACCTTACGAAGTCCAGTGTCCGATAGAGACCTTATCCTTCCAAGAGCATTAGGATCGTATTAGGCCCGTAAGTTCTGGCTGCTGCACTATATATACCTGAATGCCAAAGGTTCTTCGCGAACGGGCCTAGGGTCGGGATCTAGGTTCGGATCTAAGGTTGAGTAGAGGGAGTCACTACGCGAATTGAGCCCGTCTTGCCAAACGTTTGTATAGAGCGACTAAACTATGGATGAGCCCACCACAACATGGCCACCGTCACCATAT'
transcribe_rna(t)
```




    'GAUGUUUCCACUUAUCCAGUUGGCUGGAUGUGUUAGGCGCCCCUGAUAAUGAGCAUCAACUGAUCUGUGAUGUAGUGGAGUAGCACUUUUCGAUAAACUGUUUGAGCAGCCUGGACAGCUCGACGCAGACCCACAUGCACACGGUAUUUCUCCCGCGACCUGUCAGACUAACGAGUGCUAAAGCUUUAAGAUGAGAUAAAGUUGUCCUCUUGGGGACAAAUAAUGCGUCGAAAUUGACGUCUAUUACCCCUGAGCACCAUUAUAGAGCAUUACCAGAGGAAUGAUACUCCGGGAAGCUAUAUAAUUUCGAAUGCGUCUGCCCUAGAUACGGAAUGAGUUUAGAUGCGCUAGAGUGAUGCCGAUCUUUCCAUGAAAGUGGGUGUAACUAGUAAAGGUGUGACUGGGCCUACGUCCUCGUAUCGUGGACCUCCACACCAUGACUGCGGCGCUCUUUGGCACCGCGCGACCCUCUUAUACCUAGCUCGCGGGGUCGCUGGGGAAACAGCCCCUAUGGGGGCAAGCGACUCAGUAAAUCUCUGUUUACGUACCUAGUUAUUAUGACCAGAUGCAUUCUGUCAUCGAUCUAAUCGCAGACAAUUCUCCUCUCGACGGGGUCGUAAGUCCAGCGACUGGGGCAGACCGCGCUCAGUAUACAUCUAAGGGUAUGGUGCGCUGUAUAAUACCGCCAUAAAACCUUACGAAGUCCAGUGUCCGAUAGAGACCUUAUCCUUCCAAGAGCAUUAGGAUCGUAUUAGGCCCGUAAGUUCUGGCUGCUGCACUAUAUAUACCUGAAUGCCAAAGGUUCUUCGCGAACGGGCCUAGGGUCGGGAUCUAGGUUCGGAUCUAAGGUUGAGUAGAGGGAGUCACUACGCGAAUUGAGCCCGUCUUGCCAAACGUUUGUAUAGAGCGACUAAACUAUGGAUGAGCCCACCACAACAUGGCCACCGUCACCAUAU'



The time complexity of this algorithm is O(n) where n is the number of symbols in the input string "t", the algorithm needs to iterate through the entire input string once, so the time it takes to run the algorithm grows linearly with the size of the input. Space complexity is also O(n) since the algorithm needs to keep all the n characters of the input DNA sequence in memory.

It's a simple and easy to understand solution, but it also has a very low time and space complexity which makes it efficient and fast.
This algorithm is a single step process, no loops or conditions, the string method replace make the replace in a single step.

# The Secondary and Tertiary Structures of DNA

We introduced nucleic acids, and we saw that the primary structure of a nucleic acid is determined by the ordering of its nucleobases along the sugar-phosphate backbone that constitutes the bonds of the nucleic acid polymer. Yet primary structure tells us nothing about the larger, 3-dimensional shape of the molecule, which is vital for a complete understanding of nucleic acids.

The search for a complete chemical structure of nucleic acids was central to molecular biology research in the mid-20th Century, culminating in 1953 with a publication in Nature of fewer than 800 words by James Watson and Francis Crick. Consolidating a high resolution X-ray image created by Rosalind Franklin and Raymond Gosling with a number of established chemical results, Watson and Crick proposed the following structure for DNA:

The DNA molecule is made up of two strands, running in opposite directions.Each base bonds to a base in the opposite strand. Adenine always bonds with thymine, and cytosine always bonds with guanine; the complement of a base is the base to which it always bonds. The two strands are twisted together into a long spiral staircase structure called a double helix.

In light of Watson and Crick's model, the bonding of two complementary bases is called a base pair (bp). Therefore, the length of a DNA molecule will commonly be given in bp instead of nt. By complementarity, once we know the order of bases on one strand, we can immediately deduce the sequence of bases in the complementary strand. These bases will run in the opposite order to match the fact that the two strands of DNA run in opposite directions.

**Problem**

In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C' and 'G'.

The reverse complement of a DNA string s
is the string sc formed by reversing the symbols of s

, then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

Given: A DNA string s

of length at most 1000 bp.

Return: The reverse complement sc
of s

.
Sample Dataset

AAAACCCGGT

Sample Output

ACCGGGTTTT



```python
def reverse_complement(s):
    # Create a dictionary to map each symbol to its complement
    complement_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    # Reverse the DNA string using slicing
    rev_s = s[::-1]

    # Initialize an empty string to store the reverse complement
    sc = ""

    # Iterate through the reversed string and replace each symbol with its complement
    for symbol in rev_s:
        sc += complement_map[symbol]

    # Return the reverse complement
    return sc
```


```python
# Test 
s = "AAAACCCGGT"
print(reverse_complement(s))  
```

    ACCGGGTTTT



```python
s = "ATGAATATCGAGGTTCGTCGTCTCAACACAAAGGTACGCAACCGGAGTGTTGCACGGCCGTATTGGTTTCAGCTATATTTGCATGGACCACTCTTTGTGGATCCGCAAGTATTACATACGCTCTCTAGTCTATCTATCTGTCATTATAAGCGCATTTCATGAATGGTGCTGACTCCGCGTGGTAGAGCATAGTATAGCTGAGTGACCCGTCACAGTCACCGCCAGACGGGTCCTTTAGACTATCATAATGAAGCCCACAGAGGTCAGACGTCTGTAGTAGGCAAGTACCCTCACGCGACCGTCACCGAATTGTACGAACCATGGGCTCTCATAGGTTCGCAACTAGCCGGCACGAGAATAACCCCCGAGTTGTGTTACATTCCGGGGGAATCAACAGACCCGGAGTGATACGTGGAAACGAATCCACGACATACGGAAGGCTTTAAGAACAAGCGGGTGCTTAATACCGCAGGTCTTCTGACTCAAGATAGTGGCATAAACGGACTCGCGGTAAGGTTCCTGCTCGCCTAGCTGATAGCAGCGTGTGGCTTCAATAGGACTATGCTCACACTCGTCCAAAGCTACTTCAAACCCTTGGCAGTGTCCACTTGGTCCTGAAGCCAACTTAAGCGAAAACTAATTGGAGGCGGACAGAAGACGGCGCGCTGGCCTTAGGGGAAACATTTGGCAGCAGTAGTACTCCGGAAACCATGCGTTGCCTAACAATATGTTCTGATTGCTTACGATGCCGGTCGGTCGAGTTCACCTACACGACATTAGCCTCAGGATTTCGCGCAATGACCCGTATGATGTCCCCCTGCGCAGACCAGTCGGG"
print(reverse_complement(s))  # Output: "ACCGGGTTTT"
```

    CCCGACTGGTCTGCGCAGGGGGACATCATACGGGTCATTGCGCGAAATCCTGAGGCTAATGTCGTGTAGGTGAACTCGACCGACCGGCATCGTAAGCAATCAGAACATATTGTTAGGCAACGCATGGTTTCCGGAGTACTACTGCTGCCAAATGTTTCCCCTAAGGCCAGCGCGCCGTCTTCTGTCCGCCTCCAATTAGTTTTCGCTTAAGTTGGCTTCAGGACCAAGTGGACACTGCCAAGGGTTTGAAGTAGCTTTGGACGAGTGTGAGCATAGTCCTATTGAAGCCACACGCTGCTATCAGCTAGGCGAGCAGGAACCTTACCGCGAGTCCGTTTATGCCACTATCTTGAGTCAGAAGACCTGCGGTATTAAGCACCCGCTTGTTCTTAAAGCCTTCCGTATGTCGTGGATTCGTTTCCACGTATCACTCCGGGTCTGTTGATTCCCCCGGAATGTAACACAACTCGGGGGTTATTCTCGTGCCGGCTAGTTGCGAACCTATGAGAGCCCATGGTTCGTACAATTCGGTGACGGTCGCGTGAGGGTACTTGCCTACTACAGACGTCTGACCTCTGTGGGCTTCATTATGATAGTCTAAAGGACCCGTCTGGCGGTGACTGTGACGGGTCACTCAGCTATACTATGCTCTACCACGCGGAGTCAGCACCATTCATGAAATGCGCTTATAATGACAGATAGATAGACTAGAGAGCGTATGTAATACTTGCGGATCCACAAAGAGTGGTCCATGCAAATATAGCTGAAACCAATACGGCCGTGCAACACTCCGGTTGCGTACCTTTGTGTTGAGACGACGAACCTCGATATTCAT



# Wascally Wabbits

In 1202, Leonardo of Pisa (commonly known as Fibonacci) considered a mathematical exercise regarding the reproduction of a population of rabbits. He made the following simplifying assumptions about the population:

        The population begins in the first month with a pair of newborn rabbits.
        Rabbits reach reproductive age after one month.
        In any given month, every rabbit of reproductive age mates with another rabbit of reproductive age.
        Exactly one month after two rabbits mate, they produce one male and one female rabbit.
        Rabbits never die or stop reproducing.

Fibonacci's exercise was to calculate how many pairs of rabbits would remain in one year. We can see that in the second month, the first pair of rabbits reach reproductive age and mate. In the third month, another pair of rabbits is born, and we have two rabbit pairs; our first pair of rabbits mates again. In the fourth month, another pair of rabbits is born to the original pair, while the second pair reach maturity and mate (with three total pairs). After a year, the rabbit population boasts 144 pairs.

Although Fibonacci's assumption of the rabbits' immortality may seem a bit farfetched, his model was not unrealistic for reproduction in a predator-free environment: European rabbits were introduced to Australia in the mid 19th Century, a place with no real indigenous predators for them. Within 50 years, the rabbits had already eradicated many plant species across the continent, leading to irreversible changes in the Australian ecosystem and turning much of its grasslands into eroded, practically uninhabitable parts of the modern Outback. 

**Problem**

A sequence is an ordered collection of objects (usually numbers), which are allowed to repeat. Sequences can be finite or infinite. Two examples are the finite sequence (π,−2–√,0,π)
and the infinite sequence of odd numbers (1,3,5,7,9,…). We use the notation an to represent the n-th term of a sequence.

A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms. In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal to the number of rabbits that were alive two months prior. As a result, if Fn represents the number of rabbit pairs alive after the n-th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence relation Fn=Fn−1+Fn−2 (with F1=F2=1 to initiate the sequence). Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago.

When finding the n-th term of a sequence defined by a recurrence relation, we can simply use the recurrence relation to generate terms for progressively larger values of n. This problem introduces us to the computational technique of dynamic programming, which successively builds up solutions by using the answers to smaller cases.

Given: Positive integers n≤40 and k≤5.

Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).

Sample Dataset

5 3

Sample Output

19



```python
def rabbit_pairs(n, k):
    # Initialize the array with the first two elements
    rabbit_pairs = [1, 1]

    # Iterate through the array, starting from the third element
    for i in range(2, n):
        # Use the recurrence relation to compute the value of the current element
        rabbit_pairs.append(rabbit_pairs[i-1] + k * rabbit_pairs[i-2])

    # Return the value of the n-th element
    return rabbit_pairs[n-1]

```


```python
print(rabbit_pairs(5, 3))
```

    19



```python
print(rabbit_pairs(31, 4))
```

    1117836738901



# Identifying Unknown DNA Quickly

A quick method used by early computer software to determine the language of a given piece of text was to analyze the frequency with which each letter appeared in the text. This strategy was used because each language tends to exhibit its own letter frequencies, and as long as the text under consideration is long enough, software will correctly recognize the language quickly and with a very low error rate. 

You may ask: what in the world does this linguistic problem have to do with biology? Although two members of the same species will have different genomes, they still share the vast percentage of their DNA; notably, 99.9% of the 3.2 billion base pairs in a human genome are common to almost all humans (i.e., excluding people having major genetic defects). For this reason, biologists will speak of the human genome, meaning an average-case genome derived from a collection of individuals. Such an average case genome can be assembled for any species, a challenge that we will soon discuss.

The biological analog of identifying unknown text arises when researchers encounter a molecule of DNA from an unknown species. Because of the base pairing relations of the two DNA strands, cytosine and guanine will always appear in equal amounts in a double-stranded DNA molecule. Thus, to analyze the symbol frequencies of DNA for comparison against a database, we compute the molecule's GC-content, or the percentage of its bases that are either cytosine or guanine.

In practice, the GC-content of most eukaryotic genomes hovers around 50%. However, because genomes are so long, we may be able to distinguish species based on very small discrepancies in GC-content; furthermore, most prokaryotes have a GC-content significantly higher than 50%, so that GC-content can be used to quickly differentiate many prokaryotes and eukaryotes by using relatively small DNA samples.

**Problem**

The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.

DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.

In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.

Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).

Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.

Sample Dataset

>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT

Sample Output

Rosalind_0808
60.919540



```python
def read_fasta(filename):
    # Initialize an empty dictionary to store the entries
    entries = {}

    # Open the FASTA file
    with open(filename, 'r') as f:
        # Initialize a variable to store the current ID
        current_id = None

        # Iterate through the lines in the file
        for line in f:
            # If the line starts with '>', it is an ID line
            if line[0] == '>':
                # Extract the ID from the line
                current_id = line[1:].strip()

                # Initialize an empty DNA string for the ID
                entries[current_id] = ''
            else:
                # If the line is not an ID line, it is a DNA string
                # Append it to the DNA string for the current ID
                entries[current_id] += line.strip()

    # Return the dictionary of entries
    return entries
```


```python
def highest_gc_content(fasta_file):
    fasta = read_fasta(fasta_file)
    gc_content = {}
    for id, dna in fasta.items():
        # Compute the GC-content for the DNA string
        gc = (dna.count('G') + dna.count('C'))*100 / len(dna)

        # Store the GC-content in the second dictionary
        gc_content[id] = gc

    # Find the ID with the highest GC-content
    highest_gc = 0
    highest_gc_id = None
    for id, gc in gc_content.items():
        if gc > highest_gc:
            highest_gc = gc
            highest_gc_id = id
    return f"{highest_gc_id}\n{highest_gc:.6f}"
```


```python
print(highest_gc_content('rosalind_gc.txt'))
```

    Rosalind_9397
    52.468427



# Evolution as a Sequence of Mistakes
   
A mutation is simply a mistake that occurs during the creation or copying of a nucleic acid, in particular DNA. Because nucleic acids are vital to cellular functions, mutations tend to cause a ripple effect throughout the cell. Although mutations are technically mistakes, a very rare mutation may equip the cell with a beneficial attribute. In fact, the macro effects of evolution are attributable by the accumulated result of beneficial microscopic mutations over many generations.

The simplest and most common type of nucleic acid mutation is a point mutation, which replaces one base with another at a single nucleotide. In the case of DNA, a point mutation must change the complementary base accordingly.

Two DNA strands taken from different organism or species genomes are homologous if they share a recent ancestor; thus, counting the number of bases at which homologous strands differ provides us with the minimum number of point mutations that could have occurred on the evolutionary path between the two strands.

We are interested in minimizing the number of (point) mutations separating two species because of the biological principle of parsimony, which demands that evolutionary histories should be as simply explained as possible.

**Problem**

Given two strings s and t of equal length, the Hamming distance between s and t, denoted dH(s,t), is the number of corresponding symbols that differ in s and t.

Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).

Return: The Hamming distance dH(s,t).



```python
def hamming_distance(s, t):
    return sum(c1!= c2 for c1, c2 in zip(s, t))
```


```python
s = 'GAGCCTACTAACGGGAT'
t = 'CATCGTAATGACGGCCT'
hamming_distance(s, t)
```




    7




```python
s = 'TGCCTAACCTGCGCCTTCTAGGTTTCAGTCCCACTCATGAACAAACCCCCGTAGCCTGATCAATTGATCAAACACTTTAACAGACTATTGGTGGATGCTAATAGGGTTCTGATCTTTTCTGCTGAGCTTTACAGGCTTCCGATTGCGGTTCAAGGCTAAGCATCGGGTTATATCTTTCGGGGTACAAACCGCATAGGGATAGAGCTGTGCCCGATGCTAGGCATTCGGCAAAGTGACGATTGGCTACTAACAAGAGGATCTGGGGAAGTGATACGCTGAGCGCGTATGGTCCCCGCTCCCGTCGTTGGGGCGCCCCCTTAGTTTACAGCTCGACCTGCAGCGGAGACCCTCGCTGAATGATTGCGGCGTAGGACGGTCAACGGAAGTGTACGTATATTCTATGGCGAGCCGAAGTCGGAGATGTAGCAGTTTTGTGTCTGGTGAATTGATATGCCACAAACAGCCTCACCATAACAGGGACCCCCAGTTCGTATAGGAGGTACTTTTCTGACATTAGTCTCCCAGTTACATATGCACTCCCCGAGCGGTAGGAGGGTTCAGAAGGGAACTAAGTAGGCCCCAAAAAAGAAATACCTAGTGACGTGACCGGTTCCAGACGGTGCGTTGCCCTATTGATTCGTCACCACCTGTCACCCTTCATATGCAGTCCTAGACGTGCAGGAGGATTCGGATTATGGGCTGTATAGACCCCTAAAGTTTCGCAAGGCCTGCGTTTGAGTCTGAGGGGTGCGTCAGCCTAACGCAGCGACCATCCGCGACACCGGTGTTATCCGGTCTATCTCCAACCCCCAAGAATCCTAACGGCAAATCCACCGGCACTTTTTGCCATTGGCAACAGACCCTCCTTTCGAGTTCATGCCGACCGGCAAGTGCAACCTCGAAGTCAAGTTAGACCTGAACGACGCAGTTGGCAGCCGCTTCCGGA'
```


```python
t = 'TGACTTTCCTTCGACATTTATGTGACTGTTTGGCTTTTGGGTACATCCCTGAAGAGCGAACGACTGCCTAGAATCGTATCCTAGTAAATGTTGCATGCTATTATAATAGTAATCGATCCTGCGGAGGTTACCAGTTTGCAGTCTGCGGGAAAGTGATAAGTGACCCGTGATGCCACTCAACGTACAACCTCGATTTGGATGACGCTGAGTGTGAGATAAGGAATCCACCCAAGTTTCCCTTATTTAGTAGCAAACGGAGCTGGAGCCGAGTAAAGCTAAGTGGTTAACCACTGAGTTCAGGCACTTGAGGCGAACCCTCAGGTCAAGGCTCGACTTTAACCGATGCCTCCTGAAGGCTGCCAGCTGACGGGTACGATAAGGGGATGGGTCTCTATTGCCTCTTGTGGGACGGGCTCAGAGCAATGTAAGTAGGGTATCATGTACATTCAGCTTGAACACCCAGCCCGTTCATACTGGGCCCGGCGAATATGTGCAGGAGCTACCCTGTTCTCGGGGGGCTACCTGACACATCTAAACACCGCGACCCGTAGCTGTGTTCCGAAAGATACCGAGTACTTGACCGGGGCGAGAACCCCACTGTAGCTATCCGGTCGTGCGCTTTTGTTGACTCAGATGGCAATCCCCATTAATCTCCCTGCATATCGATGCCTTCTTCCTAAGGAGGGTTCTGATTATGATCCGTATAGTACGATACAGCGCTTCACGGCCCGGGGTTGGGCCTAATGAGTCCGTGCTTCATGCGGAGCGCCCCGCAGCGACACTGGAGCTCTATGCTCTTACTCGAAACCCGAATACTCCCTATGGCGCAGAGTCCGAAAATGAGTCACATGGGGAATAGTTGCTCCACACGCGGTATCGCCGAACGACTCTGGGATCTTCGTCACCCGAAGCGCCGCCCGATCCGCACGTCGGTGCCGGCTCCGAA'
```


```python
hamming_distance(s, t)
```




    447



# Introduction to Mendelian Inheritance

Modern laws of inheritance were first described by Gregor Mendel (an Augustinian Friar) in 1865. The contemporary hereditary model, called blending inheritance, stated that an organism must exhibit a blend of its parent's traits. This rule is obviously violated both empirically (consider the huge number of people who are taller than both their parents) and statistically (over time, blended traits would simply blend into the average, severely limiting variation).

Mendel, working with thousands of pea plants, believed that rather than viewing traits as continuous processes, they should instead be divided into discrete building blocks called factors. Furthermore, he proposed that every factor possesses distinct forms, called alleles.

In what has come to be known as his first law (also known as the law of segregation), Mendel stated that every organism possesses a pair of alleles for a given factor. If an individual's two alleles for a given factor are the same, then it is homozygous for the factor; if the alleles differ, then the individual is heterozygous. The first law concludes that for any factor, an organism randomly passes one of its two alleles to each offspring, so that an individual receives one allele from each parent.

Mendel also believed that any factor corresponds to only two possible alleles, the dominant and recessive alleles. An organism only needs to possess one copy of the dominant allele to display the trait represented by the dominant allele. In other words, the only way that an organism can display a trait encoded by a recessive allele is if the individual is homozygous recessive for that factor.

We may encode the dominant allele of a factor by a capital letter (e.g., A) and the recessive allele by a lower case letter (e.g., a). Because a heterozygous organism can possess a recessive allele without displaying the recessive form of the trait, we henceforth define an organism's genotype to be its precise genetic makeup and its phenotype as the physical manifestation of its underlying traits.

The different possibilities describing an individual's inheritance of two alleles from its parents can be represented by a Punnett square.

**Problem**

The probability of any outcome in a probability tree diagram is given by the product of probabilities from the start of the tree to the outcome.

Probability is the mathematical study of randomly occurring phenomena. We will model such a phenomenon with a random variable, which is simply a variable that can take a number of different distinct outcomes depending on the result of an underlying random process.

For example, say that we have a bag containing 3 red balls and 2 blue balls. If we let X
represent the random variable corresponding to the color of a drawn ball, then the probability of each of the two outcomes is given by Pr(X=red)=35 and Pr(X=blue)=25.

Random variables can be combined to yield new random variables. Returning to the ball example, let Y model the color of a second ball drawn from the bag (without replacing the first ball). The probability of Y being red depends on whether the first ball was red or blue. To represent all outcomes of X and Y, we therefore use a probability tree diagram. This branching diagram represents all possible individual probabilities for X and Y, with outcomes at the endpoints of the tree. The probability of any outcome is given by the product of probabilities along the path from the beginning of the tree.

An event is simply a collection of outcomes. Because outcomes are distinct, the probability of an event can be written as the sum of the probabilities of its constituent outcomes. For our colored ball example, let A
be the event "Y is blue." Pr(A) is equal to the sum of the probabilities of two different outcomes: Pr(X=blue and Y=blue)+Pr(X=red and Y=blue), or 310+110=25.

Given: Three positive integers k, m, and n, representing a population containing k+m+n organisms: k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.

Sample Dataset

2 2 2

Sample Output

0.78333



```python
def dominant_probability(homoDominant, heterozygous, homoRecessive):
    total = homoDominant + heterozygous + homoRecessive
    dominantTotal = homoDominant/total
    heterozygousDominantTotal = (heterozygous/total) * ((homoDominant/(total-1)) + (heterozygous-1)/(total-1) * 0.75 + ((homoRecessive/(total-1)) * 0.5))
    domHeterRecesTotal = (homoRecessive/total) * (homoDominant/(total-1) + heterozygous/(total-1) * 0.5)
    completeProbablity = dominantTotal + heterozygousDominantTotal + domHeterRecesTotal

    return completeProbablity
```


```python
homoDominant = 26
heterozygous = 30
homoRecessive = 25
dominant_probability(homoDominant, heterozygous, homoRecessive)
```




    0.7581018518518519




```python
homoDominant = 2
heterozygous = 2
homoRecessive = 2
dominant_probability(homoDominant, heterozygous, homoRecessive)
```




    0.7833333333333332



# The Genetic Code

Just as nucleic acids are polymers of nucleotides, proteins are chains of smaller molecules called amino acids; 20 amino acids commonly appear in every species. Just as the primary structure of a nucleic acid is given by the order of its nucleotides, the primary structure of a protein is the order of its amino acids. Some proteins are composed of several subchains called polypeptides, while others are formed of a single polypeptide.

Proteins power every practical function carried out by the cell, and so presumably, the key to understanding life lies in interpreting the relationship between a chain of amino acids and the function of the protein that this chain of amino acids eventually constructs. Proteomics is the field devoted to the study of proteins.

How are proteins created? The genetic code, discovered throughout the course of a number of ingenious experiments in the late 1950s, details the translation of an RNA molecule called messenger RNA (mRNA) into amino acids for protein creation. The apparent difficulty in translation is that somehow 4 RNA bases must be translated into a language of 20 amino acids; in order for every possible amino acid to be created, we must translate 3-nucleobase strings (called codons) into amino acids. Note that there are 43=64 possible codons, so that multiple codons may encode the same amino acid. Two special types of codons are the start codon (AUG), which codes for the amino acid methionine always indicates the start of translation, and the three stop codons (UAA, UAG, UGA), which do not code for an amino acid and cause translation to end.

The notion that protein is always created from RNA, which in turn is always created from DNA, forms the central dogma of molecular biology. Like all dogmas, it does not always hold; however, it offers an excellent approximation of the truth.

An organelle called a ribosome creates peptides by using a helper molecule called transfer RNA (tRNA). A single tRNA molecule possesses a string of three RNA nucleotides on one end (called an anticodon) and an amino acid at the other end. The ribosome takes an RNA molecule transcribed from DNA, called messenger RNA (mRNA), and examines it one codon at a time. At each step, the tRNA possessing the complementary anticodon bonds to the mRNA at this location, and the amino acid found on the opposite end of the tRNA is added to the growing peptide chain before the remaining part of the tRNA is ejected into the cell, and the ribosome looks for the next tRNA molecule.

Not every RNA base eventually becomes translated into a protein, and so an interval of RNA (or an interval of DNA translated into RNA) that does code for a protein is of great biological interest; such an interval of DNA or RNA is called a gene. Because protein creation drives cellular processes, genes differentiate organisms and serve as a basis for heredity, or the process by which traits are inherited.

**Problem**

The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.

The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.

Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).

Return: The protein string encoded by s.

Sample Dataset

AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA

Sample Output

MAMAPRTEINSTRING



```python
def translate(s):
    equivalences = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    # Initialize the protein string to an empty string
    protein = ""
    # Divide the RNA string into codons
    codons = [s[i:i+3] for i in range(0, len(s), 3)]
    # Iterate over the codons
    for codon in codons:
        # Look up the codon in the table of equivalences
        amino_acid = equivalences[codon]
        # If the amino acid is "Stop", stop translation
        if amino_acid == "Stop":
            break
        # Otherwise, add the amino acid to the protein string
        protein += amino_acid
    # Return the protein string
    return protein

```


```python
s = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
translate(s)
```




    'MAMAPRTEINSTRING'




```python
s = 'AUGGUGAUCGAUUGCUGGAACCCUUUCUAUAUAAGUCGAUAUACUUUAAAGCUCACAGCAUCAUACAGUACGCUAACGGUGGGGCGGACCCUCAUCCCGCUAUCUGGCUCCCCCGGAGCGGAUAAUGAUAGACUACGCAUACAUGUUAAAAGACGCUCCAUUUUUAGACCAUUGUCCAGUGUCGCGACACGCAUCCCUAAGUCGAGUUCCAGCCUACUCUUUCCGUUAACGACGGACUGUUUCUGCGGCAAGUUUACGUCAAACCUCGCAUGGACGGUACAGCUGGCCUGUCAUGAGUUAGCAACACCUCUCUCAUACCCACCUCAAUCUGAAGCCCUAGCAGGGCUUCUACCCCCAUUUAAUCAUAACUAUUUGACGUCGCGCGACACGUAUAUUAAAAGUACCCAUUUAUGUACACUCACCCCAUGGCGAAAGGGCAAAGGUUACGACCAAUGUCUGAUAAAAAUUCGUUGGAUCGUCGUGCGUAGUAGAGCCCUUCUCGGUUACUUGGCCGUUCUGACACGUAACGCAACCACGCACGGGUCCGCGUACAAGAGCCGACGCACUGUGCGAGAACGAUCGCCCGGACUGAAUCUCGCUGCCCAGGCGCGUUCCGGCUUCACAAACGGACGUCAAUGCAACAGCUGUGACAUUCGGUCUCCGCCACCUGUUCCCGAUCAUUUCCCUGUUCAGCGUACGGAGGGUCGUUCGUUGAAAUUAAUUCAGGGCAAGGAAGAAUUACGGGGGGGGCGGGUAAGAACCAACUGGCAAAACGGAGUCGAAAAGAUGGCCACUGAUAUGAUGAGGUUCGUAAGGCUUCACAUUACUAAGAGUCGGUAUCAUGGUGGACAGGGGUAUUUACCGCUCCAAGUUGUUAGCGUAUUAGUGGCCAAUUUGGGUGAUGCUGUUAUGCAGUGUUUCUUGCCUGUGGUAGCGGUCCACGGUGAAAGCGCGAUAGCAGCUAGCCUUACUUCACCGAUAUCGCCGUGGGUCACUGGCAUGCCGACGUUGACGGCGCGCCGAGAGCUCCACGGGACCUUAUGUCGCGCUCAGAGGGUACACGACAUAGCUUGCUUCACUUGUAAACCCCGCCUUUGUGAGCCGAACAGAGGUCGCCUAAUCGUCACCGCGACUGCUGGACCCGCGGGCCAACCUUCGAGCAAAGAUGUCACUGUGGCGCCAGGCCUCUGUAAUCGUUGGACAAUGAUACGCUUAUAUCGAUUAUACAGGAGCCCGUUAGGUACGCACAGCGUAGUCUGUGUAACGAUCUGGCGUACCUGUUUGCCGGCCUUGGGAUUAGUUCACGCGUAUCUGUCGACACCUUUAAUCCGAAACCAUCGUUCUAAACCUUAUCCCACCGGUCUGUAUCCUUCGGGCAACGGUCUGUGUCGGACGCCCGUGGCAACGUUACAGCAUGCUUACAGUACUUUAUACCUAGUUUCUGCAAUACAAGUAGACCUUAUUUUUCCUCAUCUCGCAAUAUGGAACGUUGAGACGGCACUACUCAUGUACGCCCCAUAUACCGCCAAGGGAGUCCCCUAUGAACCCGAGUGCCUACGGGAUAUCCCUAACAGGGGUCCAACUGACGAGAAACGACCAAGCUACGGUAGCUCCAAAGCUCCGUUUUCACUACGACGAGUCGUAGUUAAGCACUUAGACAGGCCCAGAGCUACUUUACCUUUCAUUUACCGGGUCAAUCCAUUUUCAGAACCAUCAGUUUGUUAUUCGUGGAGCUUGGGCAGUUACGCCGUAAUCGUUUUCCGCCACUCAGCUUUAGACACAAGGGCGCUACAAGGGCGUGCUGGACAUACCUCGUCUCUGGUCUCCAUCCGCAGCCUCGGUGACGAAGUGGAUCUUUGUCCGCGUGCAUUGCCAGGCGUAACAGGAGAAUGGUGGUUACAAAGGGCUUUCGUACUCACACCGGUAACUAUCGACCUGUUAGCCCGCAAGCGCAACCCUGCAUGGUGGAACGGCAGAUCACACCUAACCCGGAUGGUGCCCUACACGACCUCGUCCGAUGUGUUGAUUUGCAGAGGCGCGUCCUCCAAAUCAAAGUUUACUAUUAGUUGCCCGAGGCGGGCACAUGCAGCGAACGUGCUCAUAUAUCCCAAGUACGUCCGUAGGUACACGAUCUUUAUUGCUAUCAGAGACGGUAUGGGCGGCUGCAAUCUAGCAAGAGGGCCUCUAGGAGAAGAAUCGAACUUACGAAGGUCAGUCGAAAACACUUAUCAUAUCUUACACGCGUUCAAUUGCCUUUCCAAGGAUGACCGGAUCACGUCACAGCGGUGGUGUCCGAUAAGGACCGCUUACUGUAUAAACCAUCGUAGCAGCUGUCAAAGUCGUAGUAGGAGACGGUGGAUUCAUCCUAUCUCCUACGAUUCGUGUUUACUCCGCGGUCGCUAUUCAGUGCGUCUAUGUCAGGAACGGAUUUUUUGCACUGCAAUUUCUCCGGCAGUGCACAGUGUUGUAAGCGAUUACGGGAAAUUUGUAUGCCAGCCCGUAGGUGGGGCAUGUUUGACGGCUGGGGCAUGCGUUGGCCGCAAUUUAACAGUGGAAUUGUGCCAUAGGGUCUUUCUAAUUUGUCCGGUGCCCUUUCGUUUAUCCGCCCAGGUGCUCGUUGAUAUGUGGACACAUCGCGCCGAACCCUGCGUUUUCCUUUGCUGGCUCUUGCGCCAUCACAAUCAGCUAAAAACUGGUUUCUCUGAAAAGCGAUGCCUCGGAGACUGGAAUUCUAGGGUCGCUCAGCUACGGUUACUAGCAUGGGACGAUCGCAACCAGGCUGGCAUACGCUUCACGUACCAGAGGACGCCAGAAACCUGCGUUUCGAGGGGACACCGAUCUGGCCUGAAUGUGUACCCAGGAGAAACCCGUGUCCACGUAAUACUCCCCCCGAAGUAUAAGCGCUGUGGUUGGUUUGAGGAGAGGAAUUUCCUUACCAACCAGUGCCAGGAAACAAUCCACUCAUUUGACUAUACACAUCUGGUAGGGAUACCUUCACCAAGCAUAAGUGUGCUUUCCCUUCUCCGAACUUUUGGCAUGUUACCAGUUCAGAAAUGGAACUACCACGCCGACACCCCCAUGUUUCUAGCUGCGGUGUGGAUCCAAGUGCCAUGGGUACAGAACGAAGAAAUUUGUAGCACUUUCUGUGAGAGCCGGUCCUUAUAUGUGACCUUCCAUGUGGGAAAGAGGCUAGGACAACAUUCCCCAACUAGACUCCAUCCGAGGGUAUAUCAUCCGUGCGUGUGUCGGCAUAGGAGAAUCAGCAGGCUGCAAAGGGACAUUUGCUCCAUGAGACAAGAGUAUGACGCAUGUGUGGCCCCAACCGCUGGGAUAGGUAAGUUACUCACACAAUGCUUGAGACUAAAAAAACCAAACAAUGCCGUGCGAGUGGAUUGCGUCUCGUCACUACCUAACUCCACGAGUCUAAUCUUAUUGGAGGGUAUCUUCUUCUUCCUGUGUGAGAACGGGUCACCGGGACGAAAAAUCGAGCUUGAGUUAUUUGAAAAAUCCACGGCCUUUAGUUUGAAUCCACGCGACUCGGCCAUGAUUUUCGGCCCUCCUCUACGGAUAUCUGAUAUUCCCCAAUUAGUGUGGGACACCAAAACGGGUAAUAGAGAUCAUAGGGAAAUUCCCUGUGUGACAUUAUGUCUUUUGAAGCAAAAAACAACCAGACCUUGGGUGCGUCGUUGGGAUCACGCGUCAGAACUUAACAUCCGUCAUCCGGGGGGCGCCUGCGUAGACGCUAAGCCGGUCCGUGCCCGCACCUUUGCAGGUGUUGGCUCCCCGGUAGGCCAUGGUUCACGACGCUUCGGUAGCCUACACUUCGAAAAGAACUUCAUAUGGAUCGGGGAUUCCCCGCGGGAUCUGAGGAGUGGUGUAGUGAAAUGUAACGCCAGCCCCAUUGCUCUGAAGGCAGUGGAAAGUGCCUUAGAUGGGUCACGCGCGGAUGGCGGUUACGGGCGAGAGGUGUUCGCGGGAUGGCCUCUUCGUCCCCGCGUCACUACGCGACGCUAUACUGCACAUUCGCAUAAUCCUCGUUUACGCUCCCAGUACGCUAACCGUUUCCUCUUUUGCUAUAGUUAUUCGCAGAAAACAGGACAUACCCGGGGUCAGAUUCGGUGGCGAAAGGAGGGAUUCCUCAACGGGCCCCGCCUGUGUGCUAGCCGUUUUGGCUUGAUGGGAUUAGGAGGAAAUUUGUUUAGGGCCGGAGCUAAAAUUGGUGCAUACAAAUCGGGUUAUUGCUAUGUAUCUUGUAAUUCCUACGGUCCCCCGUGCGACGAAUGGGAGUUACUCCUGGGGACUUUCUAUCAAGUUGCCGGCGGAACAAAAAUCGGAGUAUCGAUUAGCAACCCGGUUGGGCAGCGAGCUCCAAAGGCGAGAGUGUGCAACUCCAGAUUGCGGAGGACGAUAAGACUCCCUUCUGGUUGCCCCGGAAAGUACUGGGGUUGUGCUACGGGUCGCCACUGGAGCGAUUUGGAUCAGUUGACAACCUGGGCGCAAGUAAUACCCGUGUUGCAAGCAUGGCGGGCUAUUGUACAGCGGGUUUUCACGGGAACCACUGCGUUGACGCUUUUACGAUCUCAGAGAGCUUGCCGAGAGAUACCUCUGCGCAAUCGCGCCGUGAAGGUGUGCGGCUCUUCCUGGGAGAUGAAAUCGACAGAUACGGUGCUCCAUAUCUUCUGGACAGUCAUGGCACUUAUCCUGACUUGGUUUCAUUUUCGAAUUGCGACCUACAGGCUGGGAGCUUCAACGGCGCCCUCCCCAUGUUCAGUCCCGUUCCGGGGACGAGUCCAGCCAGGGGAACUAUGUCCAGCAGUGACACGCGAUCAAUGGCAUCUGGUUAGGAGGCGCAAGAAGUGUGACUUAGCGUUUGAUAAUAGCGAACUUGGCAAAGCCGGGCUUAAACCCUCAGCAAGAGUUCCAUACCGUGCAAUGUCUACCGGCGUCGCAUUGUCCAAGGUGGGUUUGAUCCGUGCGAGUGCGGUGCCCAAGUUUUAUGCCUAUCGAAGAGACAGCCGCCUUCGCGUUACCAGGCUCACUACAUACCAUACACUAGAGACGGUUUUGUUCUACUGCGUAGAGUGCGUUCCUAAGCUUCUCCUCUUGAAUUGCGCAAUACAUCCCGGGUACUGCACGUUAACUGCUAUAUCGCCGGCGUUUAGACCUUGGAGACUUUCGGCAGAGCAUCUUCAAGGUCCAUACGAGCACUGCACGCGUGCGGCAGAAGGGUGCAGGGGUCAAGCGAGUCACGUGAUGCAACCGUUUCCCAAAAUGGAACGCGAUUACCCUAAUCUGGAGAGGGCGCACCCACGUACUCCUCGCGGAGACCUGGCGCAAGUCUCGGUCCACCAUCUACUUCGCGGACUCGAGCUCAAAAGCGAUAACUCUCUCACUGAGACAUCUAGUACUAAUCUAAGGCCGGCGCAUAUGCAAAACAUUAUCACACUCCAAGCCUGCUGCCCUAGUCGCUUCCCUAGGAAUGGAUGGCCUGGCACGAAUACCACCGGUUGUAACUACUCGCAUCAAAGAUCUGAAAACCACGGAUCUAGCUCCUGUAGCCUGGCCCUCACCGGUGGCGGGACAGUAGGGAGCAUAUCAUGCUUGACCACUAGAGGUAUACCGCGGUUUCGACUAAGAAAGACAGCAGAAGGACGAGAGCACAUAAUGAACGAUUCGGGAGAAGAGGUGACCGCUGAAACUGCUCCAACAAAUUACUUAAGACCCAGAGAGUUGCAAACUCAUUACUAUCUAGUCCUCCGUCACGGUGUCGCUACUAACAAUAAGUCAUCGGGCUGUCCUCCUCGGUGGGUUUACGUGGUUCAUGUGCCGGACUCUAUGUAUAUGGCUGGGCCACAACCAAGGAGUCAAUUAUGGAACUCGGCGCCAAACAGCGCUCUACAACCGGGCGUCACCCCUUUGCUGGGGUGUCAUCAAAACGUGAGGUUCAACUGUUCUGCCGACGCAGGGCGUCAUGCCCGUUUAACUACAGCGCGAUAUGUCCACAUAAUAAACUGCGCUGACAGUGAGUGGAAGCUGUGCGAACGGAUAUUGAUGGAUAGCCCGGCUAUGGAGCGCGUUGCGGCUUCGGUCAAUAUUACAAUCUAUAGCUGCCACAAAAGACCAUCACGAUAUGAGGAUUUGGAAGUAAACGUUCUCCGUGAACCGUCGUAUUACUUGGUGCAGUGCGUCGCUAUAGGUGAUAUUUACAUUUCAGAAAUCGAGCCACACAUUCUUCGACUCUGCCUUCAGAAAUGCCAAUGCGCUGCGCAGGUCCAUAGACGUGCAUCCGCUGUACCAAUCAGUCACCGGAAGCAGCUCCAUCUGGGUAGUAAUUACACUUCCUCACUACUAACGGGGUUCCCCCCUAGCUCGCUUCACAGUCCGUUGCGCUGUACGAGCUUAUUGCUUACCCCCAUAUUAGGAACUGUCCUGAAUCCCUCCUGUUCCCAUUGUAAGGAACAAACACUCUUGCUGUGUCAGAGCACUCUGAUGUCACGGUCUGUUUUAAUAGCAUCAUGUCGCAUGACUUGUAGCACCGGUUGUCGAGCUCAUGACGGUUUCGCUGACGUCAAUCGUAUAGCGAUCGUCGUCCGUGCAAGAUUCGAGCAUUUCAUGCGCCAAUCAGGGACAAUUGACGCUGAGCUAGGAGUUGAGGCGUCAGGUUACGUACCGCUCUUAAGAAAAAAGUGGUGCAUCCGAGUGUUACAGUGGAAUAACGAAAAGGGUGCCCCGCCUACCUGUGUAUGGAUUGUUCGGAAAAGAGGUGACUGCACUGCAAAAAAGACGAACGAGGAUACUGGUUCUAUUCCAAAGACAGCCUGCUCGGUUAAGGCAUUUGCUAAGUGUACAUUUACACGGUACGGACUCAAACCAAUGCAUGUUCACUGUAUAGGUAAGUGGACUGAUAGCUGGACGCGUGAACCGCAAUCCGGCCCAAAACACGGGCACUCACGCGACUGCCGAUUACUUUCGAUAUGGUGGGUGGACCGACUGCCCACGUCGGUAGCACUACGAGGUCCCGCUUCCUGCAACCGGGUCCGUCUAAUUCUCUCUAUGGUCUGCCUAAUAGGUUUCACUGGACGCGACCACGGUCUGCCUGCGCGUCAUAAUUGGCGCAUUGGGGCCAAUGCUUUUCGCCACGAUGGUGACUCUUCUGGAGGACACACCACUGCCAAGGCAUACGACGAUUAUCGAUAUCCAAUUUACUGCUGUAGUCUAAGUUAUUUCCGGAACUCGACAGCCUGCCGCUUAGAGUUACUGCGGGACGAUCUUGGCAUGACGAUCACUCUCGAGAUGGUUAUGAAAGUAAUUGAGGUCCCCGAGACCCUUACGUGUGAAGAGCUAGACUCACUAAAGACUAUGGGUCGUUAUGGUUAUAAAACCAUUGAAAUGCAACACCGCGAUCAUGUCCGAGCGCCCCAUCGACGCAAAGAACAACGCACCACCACUCAUAGGGUGCCUAGACUGGGACUAUAUUCGAUUGGACCAAUUCGACAUACACGAUUACGAGCUCAUAGACUAGUGCUCAAAAGUGCCUCUAUCGCAUAUUCCGCGAGGGAGACACUUUCCCAGUGGUGUAAGUCGGGCCGCUUAGAGCCAUACUGUUUGCCCCUGUUUUGCUCUGGAAACCUGAACGAGUAUGCAUGCCCAAGGCGUGUGUUUCACGGCAUUACAAAUAGAACGACUACUUCCGUCGAUACUGAUGGCACAGACUAUUAUGGACUCGUUUUUUCACAGAAUUCCAUUACCUCGUAUAAGUUAUUACUCCCUCAUCCCUGUGACGCAGGAGCGACGUGUAUAAUACAGCUGGGAUUAUGCCCUCUGACUGGGUUCCCCAUACCGUGGCUGCCCACUGGUCACGAGCCCAAGGAUAUUGACCUUGCCAAUGCGAGGGCUGAUCCAGCCAGAGAUCGCGACGAGUUCUCCCCGUUUACAACCAUACUGUCGGUAUUCAAUAGUGGGACGUCGAGCUCGCGGUCCGAUACAUAUAAAGCUUAUGCUACUUCCAAUAGGCGUCCGCCAUACCCUGGCCUCGCUUAUUCCCACCUUGUCAUCUGUCCAAAACACGGGUCGAAUCUACAGUCAAACGUCAUUAGGUUUGAUGACGUCAAGAGUGGCCUUCAUUAUGGUUUCGCUCCAGGGGUCAUCAACAUUCAUCCCCUAUACCACAGAAGUAAAGAAAGCCCUGUUCUGAGAACCCACCGAGCCACCACCGACCCGGUAACGUUAUUACGAUCCACUGGCAGCUCGUCCACUAUGAUAUGCCAUUAUAGCGACGUAGGCUACCGCACUAUGAAGUUAGGAUUCGCAUGCCCGUCCCCCUAUUGGUUGGUACCACAAGUCGAGAAUUGGUAUGUGAUGAUACGGAGAGAACUGUACUUCGGCAGCGAUAGUCACUGCAACCCGUUACUCCUUGAUAGCUUUAUAGGGUCAAUAAACCUACAAGAAUCAGUAGUUUACGUUGUGAAGCGGGAGUCAUUGGAUAGAUACGACCUGUCUCAGCUGAGGCUCCCCUCGGCAAAGUACCGUGUUAUCGCGACAAACCAAGCGAGGCCAGGCUGUGUUACGAAGACGGACCGUGUGCGGGCUAGUCUGCGCGCAGAAACGCUAGACUCCAUUGAACGGUUGAGAAGUCGUCAAGCCGUGUUCGAACUUACAAUUCACCUUCACCGUUAUGAAAUUCAACACCCCACCGCUACCUUCAAGACGUGUUACUGGUUCAGACGGGACGGCCCUGCUGUAAUGCAGGCCACAGUCGGCGUUACACCAAGUUGCAACAUCCAACCUACUAAGCGAAGACCACAUCGUGAAAGCGUCUACAAAUCUGUUAACGUUGUUACUGUAAUAACACAGGCAAGUGCUCGACCCAUUCAAGCCAAAAUCUCGAAAAGUGGCGUCGCCACAGCGCCCGGCAAACAAUUAGCGGCGAUACGCAUAUCUCAUAGGGAAAGAACAUUAGUGGCCGGUCUACACACGCUGAAGUGUUGUCCCUGGGUGCCAUUGGGUCAUCACGUGGAGCCUCUCCACUUCGGAUACCGGCGAGCCUCGAUCAAGCGUAGAAGAAAGCAACGUCGGUGGGACCGAGGGAGACUGUUCCAUUGCGGGGGUGCCAGCCAUACGAGUGCCUGGUGGCGGGGAGGACGCUGCACGCUUCUGCUACUGCCUGCUUCGCUACCUGAGCUAGUGAGAAGGUUGCUUUUUUUUCCCACGGUUAAACAGGUGAUUUUAUCACGUAAACGUUUACUUACGGGCCCAAGUCCCAGAUUGUUGCGUGACGGCGGUGUCCUGCUCGUGUCGGGGUUCGGAUUUGUCCCACCAGGGCUUAGAUCCGUGAUGAAUAAUUUCCGUCCUGUACCCCGUUUCUGCGGAGAUACUGUGGACCUCGUUAGCGACUUAGUAAGAGUGGUUCCUGGGACAUGGAUACUCCAACAAGAUGACAGAGGGCUGAUAAGACUGGAACGCGCCACCUGUCACCUGUAUCACAUUAGGAGAUUACAAGAAAAUUCAAGCUGUCUUAGCGUAGAUUAUAAGCUUACUUGGCAGAAGGAAUGUUUAACAAUUGUGUCAGCGCGGGAAUCAACCAAUAACGGUUUCGGGAGCGACCCUAGAUUCAGUACCUUUGCAGUGGGUCGGGCGGGGCGGGGUUGGCUACAACCUGUUGCAGCUGAACCACAUGUGAACUUUGCAACGCGACUAACGCCAAGUGAACAGGUGGACACGAAGUUAAAUUCGACCGUUGGUGGGUUAGGGAUACAUGCUAAGUCCUGA'
translate(s)
```




    'MVIDCWNPFYISRYTLKLTASYSTLTVGRTLIPLSGSPGADNDRLRIHVKRRSIFRPLSSVATRIPKSSSSLLFPLTTDCFCGKFTSNLAWTVQLACHELATPLSYPPQSEALAGLLPPFNHNYLTSRDTYIKSTHLCTLTPWRKGKGYDQCLIKIRWIVVRSRALLGYLAVLTRNATTHGSAYKSRRTVRERSPGLNLAAQARSGFTNGRQCNSCDIRSPPPVPDHFPVQRTEGRSLKLIQGKEELRGGRVRTNWQNGVEKMATDMMRFVRLHITKSRYHGGQGYLPLQVVSVLVANLGDAVMQCFLPVVAVHGESAIAASLTSPISPWVTGMPTLTARRELHGTLCRAQRVHDIACFTCKPRLCEPNRGRLIVTATAGPAGQPSSKDVTVAPGLCNRWTMIRLYRLYRSPLGTHSVVCVTIWRTCLPALGLVHAYLSTPLIRNHRSKPYPTGLYPSGNGLCRTPVATLQHAYSTLYLVSAIQVDLIFPHLAIWNVETALLMYAPYTAKGVPYEPECLRDIPNRGPTDEKRPSYGSSKAPFSLRRVVVKHLDRPRATLPFIYRVNPFSEPSVCYSWSLGSYAVIVFRHSALDTRALQGRAGHTSSLVSIRSLGDEVDLCPRALPGVTGEWWLQRAFVLTPVTIDLLARKRNPAWWNGRSHLTRMVPYTTSSDVLICRGASSKSKFTISCPRRAHAANVLIYPKYVRRYTIFIAIRDGMGGCNLARGPLGEESNLRRSVENTYHILHAFNCLSKDDRITSQRWCPIRTAYCINHRSSCQSRSRRRWIHPISYDSCLLRGRYSVRLCQERIFCTAISPAVHSVVSDYGKFVCQPVGGACLTAGACVGRNLTVELCHRVFLICPVPFRLSAQVLVDMWTHRAEPCVFLCWLLRHHNQLKTGFSEKRCLGDWNSRVAQLRLLAWDDRNQAGIRFTYQRTPETCVSRGHRSGLNVYPGETRVHVILPPKYKRCGWFEERNFLTNQCQETIHSFDYTHLVGIPSPSISVLSLLRTFGMLPVQKWNYHADTPMFLAAVWIQVPWVQNEEICSTFCESRSLYVTFHVGKRLGQHSPTRLHPRVYHPCVCRHRRISRLQRDICSMRQEYDACVAPTAGIGKLLTQCLRLKKPNNAVRVDCVSSLPNSTSLILLEGIFFFLCENGSPGRKIELELFEKSTAFSLNPRDSAMIFGPPLRISDIPQLVWDTKTGNRDHREIPCVTLCLLKQKTTRPWVRRWDHASELNIRHPGGACVDAKPVRARTFAGVGSPVGHGSRRFGSLHFEKNFIWIGDSPRDLRSGVVKCNASPIALKAVESALDGSRADGGYGREVFAGWPLRPRVTTRRYTAHSHNPRLRSQYANRFLFCYSYSQKTGHTRGQIRWRKEGFLNGPRLCASRFGLMGLGGNLFRAGAKIGAYKSGYCYVSCNSYGPPCDEWELLLGTFYQVAGGTKIGVSISNPVGQRAPKARVCNSRLRRTIRLPSGCPGKYWGCATGRHWSDLDQLTTWAQVIPVLQAWRAIVQRVFTGTTALTLLRSQRACREIPLRNRAVKVCGSSWEMKSTDTVLHIFWTVMALILTWFHFRIATYRLGASTAPSPCSVPFRGRVQPGELCPAVTRDQWHLVRRRKKCDLAFDNSELGKAGLKPSARVPYRAMSTGVALSKVGLIRASAVPKFYAYRRDSRLRVTRLTTYHTLETVLFYCVECVPKLLLLNCAIHPGYCTLTAISPAFRPWRLSAEHLQGPYEHCTRAAEGCRGQASHVMQPFPKMERDYPNLERAHPRTPRGDLAQVSVHHLLRGLELKSDNSLTETSSTNLRPAHMQNIITLQACCPSRFPRNGWPGTNTTGCNYSHQRSENHGSSSCSLALTGGGTVGSISCLTTRGIPRFRLRKTAEGREHIMNDSGEEVTAETAPTNYLRPRELQTHYYLVLRHGVATNNKSSGCPPRWVYVVHVPDSMYMAGPQPRSQLWNSAPNSALQPGVTPLLGCHQNVRFNCSADAGRHARLTTARYVHIINCADSEWKLCERILMDSPAMERVAASVNITIYSCHKRPSRYEDLEVNVLREPSYYLVQCVAIGDIYISEIEPHILRLCLQKCQCAAQVHRRASAVPISHRKQLHLGSNYTSSLLTGFPPSSLHSPLRCTSLLLTPILGTVLNPSCSHCKEQTLLLCQSTLMSRSVLIASCRMTCSTGCRAHDGFADVNRIAIVVRARFEHFMRQSGTIDAELGVEASGYVPLLRKKWCIRVLQWNNEKGAPPTCVWIVRKRGDCTAKKTNEDTGSIPKTACSVKAFAKCTFTRYGLKPMHVHCIGKWTDSWTREPQSGPKHGHSRDCRLLSIWWVDRLPTSVALRGPASCNRVRLILSMVCLIGFTGRDHGLPARHNWRIGANAFRHDGDSSGGHTTAKAYDDYRYPIYCCSLSYFRNSTACRLELLRDDLGMTITLEMVMKVIEVPETLTCEELDSLKTMGRYGYKTIEMQHRDHVRAPHRRKEQRTTTHRVPRLGLYSIGPIRHTRLRAHRLVLKSASIAYSARETLSQWCKSGRLEPYCLPLFCSGNLNEYACPRRVFHGITNRTTTSVDTDGTDYYGLVFSQNSITSYKLLLPHPCDAGATCIIQLGLCPLTGFPIPWLPTGHEPKDIDLANARADPARDRDEFSPFTTILSVFNSGTSSSRSDTYKAYATSNRRPPYPGLAYSHLVICPKHGSNLQSNVIRFDDVKSGLHYGFAPGVINIHPLYHRSKESPVLRTHRATTDPVTLLRSTGSSSTMICHYSDVGYRTMKLGFACPSPYWLVPQVENWYVMIRRELYFGSDSHCNPLLLDSFIGSINLQESVVYVVKRESLDRYDLSQLRLPSAKYRVIATNQARPGCVTKTDRVRASLRAETLDSIERLRSRQAVFELTIHLHRYEIQHPTATFKTCYWFRRDGPAVMQATVGVTPSCNIQPTKRRPHRESVYKSVNVVTVITQASARPIQAKISKSGVATAPGKQLAAIRISHRERTLVAGLHTLKCCPWVPLGHHVEPLHFGYRRASIKRRRKQRRWDRGRLFHCGGASHTSAWWRGGRCTLLLLPASLPELVRRLLFFPTVKQVILSRKRLLTGPSPRLLRDGGVLLVSGFGFVPPGLRSVMNNFRPVPRFCGDTVDLVSDLVRVVPGTWILQQDDRGLIRLERATCHLYHIRRLQENSSCLSVDYKLTWQKECLTIVSARESTNNGFGSDPRFSTFAVGRAGRGWLQPVAAEPHVNFATRLTPSEQVDTKLNSTVGGLGIHAKS'



# Combing Through the Haystac
   
Finding the same interval of DNA in the genomes of two different organisms (often taken from different species) is highly suggestive that the interval has the same function in both organisms.

We define a motif as such a commonly shared interval of DNA. A common task in molecular biology is to search an organism's genome for a known motif.

The situation is complicated by the fact that genomes are riddled with intervals of DNA that occur multiple times (possibly with slight modifications), called repeats. These repeats occur far more often than would be dictated by random chance, indicating that genomes are anything but random and in fact illustrate that the language of DNA must be very powerful (compare with the frequent reuse of common words in any human language).

The most common repeat in humans is the Alu repeat, which is approximately 300 bp long and recurs around a million times throughout every human genome. However, Alu has not been found to serve a positive purpose, and appears in fact to be parasitic: when a new Alu repeat is inserted into a genome, it frequently causes genetic disorders.

**Problem**

Given two strings s and t, t is a substring of s if t is contained as a contiguous collection of symbols in s (as a result, t must be no longer than s).

The position of a symbol in a string is the total number of symbols found to its left, including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position i
of s is denoted by s[i].

A substring of s can be represented as s[j:k], where j and k represent the starting and ending positions of the substring in s; for example, if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5]= "UGCU".

The location of a substring s[j:k] is its beginning position j; note that t will have multiple locations in s if it occurs more than once as a substring of s.

Given: Two DNA strings s and t (each of length at most 1 kbp).

Return: All locations of t as a substring of s.
Sample Dataset

GATATATGCATATACTT
ATAT

Sample Output

2 4 10

Different programming languages use different notations for positions of symbols in strings. Above, we use 1-based numbering, as opposed to 0-based numbering, which is used in Python. For s = "AUGCUUCAGAAAGGUCUUACG", 1-based numbering would state that s[1] = 'A' is the first symbol of the string, whereas this symbol is represented by s[0] in 0-based numbering. The idea of 0-based numbering propagates to substring indexing, so that s[2:5] becomes "GCUU" instead of "UGCU".

Note that in some programming languages, such as Python, s[j:k] returns only fragment from index j up to but not including index k, so that s[2:5] actually becomes "UGC", not "UGCU".


```python
def find_substring(s, t):
    # Initialize a list to store the locations
    locations = []
    # Initialize the start index to 0
    start = 0
    # Find the first occurrence of the substring
    index = s.find(t, start)
    # While the substring is found
    while index != -1:
        # Add the index to the list of locations
        locations.append(index + 1)
        # Set the start index to the index + 1
        start = index + 1
        # Find the next occurrence of the substring
        index = s.find(t, start)
    # Return the list of locations
    return locations
```


```python
s = 'TGTAAATTTATCGCACGTTCGAAATTTAAAATTTAAAATTTAAAAAAATTTAAAATTTAGGAAATTTACGAAATTTAGCGGAAATTTATGAATAAATTTAACGGTGCACTCGTTGAAATTTAAAATTTATAAATTTAAAAATTTAACAAAATTTATAAATTTAAAATTTAGGCGCGCTCTGAGAAAATTTATAAATTTATAAATTTAAAATTTAAAATTTAAAATTTAAAATTTAAGAGAAATTTATAAATTTATGCTGAATGAAGGGAACTAAAATTTAGCAAATTTAGGAAATTTATAAATTTAGCTCTCTAAATTTAAATCCGTTAAATTTACAAATTTAATGATAAATTTATAGTGTAACTTAAATTTAAAATTTAGAAATTTAATGTGTAAATTTAAAAATTTAAAGTAAATTTAAAATTTAAAATTTAAAATTTACATCGAAATTTACAAATTTATATGTAAATTTATTAAATTTAAGGTCGACTAAATTTATAAATTTACGATCGAGAGTCGTAAAATTTAAAATTTACGGAAAATTTAAAGTAAAATTTACCGAAATTTAGCGGAAAAATTTAAGGCCAAATTTAGAAATTTAGAAATTTACAAAATTTAACCCTACATGCTAAATTTAAAAATTTAAAATTTATTCTACCATAAATTTAGTGCGGAAATTTAAAATTTAAAAAAATTTATGTCAAATTTACGAAATTTAAAAATTTAAAATTTATGAAATTTAAAATTTAAAATTTAAAAATTTACAGGAAAATTTAAAATTTAAGCTTTTAAATTTAACAAATTTAACTGTAAAATTTACGAAATTTAAAAATTTAACAAATTTAAAATTTAGCTATCAAATTTAAAATTTACTGGTTAAAATTTAATCTCTCGTAAATTTATAAATTTACAAAAATTTAAAAAATTTATAAATTTAAAATTTAAAGAAATTTAAACAAATTTAAAATTTAAAATTTACGACAAATTTACC'
t = 'AAATTTAAA'
elements = find_substring(s, t)
for element in elements:
    # Print the element followed by a space
    print(str(element) + " ", end="")
```

    22 29 36 46 116 131 157 201 208 215 222 314 367 395 403 414 421 428 522 540 631 639 675 682 712 720 736 743 750 770 822 839 859 914 931 938 948 958 965 

# Finding a Most Likely Common Ancestor

In “Counting Point Mutations”, we calculated the minimum number of symbol mismatches between two strings of equal length to model the problem of finding the minimum number of point mutations occurring on the evolutionary path between two homologous strands of DNA. If we instead have several homologous strands that we wish to analyze simultaneously, then the natural problem is to find an average-case strand to represent the most likely common ancestor of the given strands.

**Problem**

A matrix is a rectangular table of values divided into rows and columns. An m×n matrix has m rows and n columns. Given a matrix A, we write Ai,j to indicate the value found at the intersection of row i and column j.

Say that we have a collection of DNA strings, all having the same length n. Their profile matrix is a 4×n matrix P in which P1,j represents the number of times that 'A' occurs in the jth position of one of the strings, P2,j represents the number of times that C occurs in the j th position, P3,j represents the number of times that G occurs in the j th position, P4,j represents the number of times that T occurs in the j th position.

A consensus string c is a string of length n formed from our collection by taking the most common symbol at each position; the jth symbol of c therefore corresponds to the symbol having the maximum value in the j-th column of the profile matrix. Of course, there may be more than one most common symbol, leading to multiple possible consensus strings.

Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.

Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)

Sample Dataset

>Rosalind_1
ATCCAGCT
>Rosalind_2
GGGCAACT
>Rosalind_3
ATGGATCT
>Rosalind_4
AAGCAACC
>Rosalind_5
TTGGAACT
>Rosalind_6
ATGCCATT
>Rosalind_7
ATGGCACT

Sample Output

ATGCAACT
A: 5 1 0 0 5 5 0 0
C: 0 0 1 4 2 0 6 1
G: 1 1 6 3 0 1 0 0
T: 1 5 0 0 0 1 1 6



```python
def find_consensus(fasta):
    import numpy as np
    import pandas as pd
    # Parse the FASTA format and extract the DNA strings
    strings = read_fasta(fasta)
    
    n = len(list(strings.values())[0])
    profile = np.zeros((4, n), dtype=int)
    
    # Count the number of occurrences of each nucleotide in each position
    for s in strings.values():
        for j in range(len(s)):
            for i in range(4):
                profile[i][j] += (s[j] == "ACGT"[i])
    df = pd.DataFrame(profile, index = ['A','C','G', 'T'])
    
    consensus_strings = []
    consense = ""
    # Iterate over the columns of the matrix
    for col in df.columns:
        # Find the row with the maximum value
        max_row = df[col].idxmax()
        consense += max_row
        # If there are multiple rows with the maximum value, add all of them to the list of consensus strings
        if df[col].isin([df[col].max()]).sum() > 1:
            new_string = ""
            for i in df[col][df[col] == df[col].max()].index.tolist():
                new_string += i
                
            #print('['+new_string+']')
            consensus_strings.extend('['+new_string+']')
            #print('['+str(df[col][df[col] == df[col].max()].index.tolist())+']')
        else:
            consensus_strings.append(max_row)
    consense2 = ""        
    for string in consensus_strings:
        consense2 += string
        
    
    return consense, consense2, df
```


```python
fasta = 'f.txt'
c, con, mat = find_consensus(fasta)
print(c)
print(con)
for index, row in mat.iterrows():
    print(index + ":", *row)

```

    ATGCAACT
    ATG[CG]AACT
    A: 5 1 0 1 5 5 0 0
    C: 0 0 1 3 2 0 6 1
    G: 1 1 6 3 0 1 0 0
    T: 1 5 0 0 0 1 1 6



```python
fasta = 'rosalind_cons.txt'
c, con, mat = find_consensus(fasta)
print(c)
print(con)
for index, row in mat.iterrows():
    print(index + ":", *row)
```

    CGCAAATCGCATTACACTGATCCTGAAAGCCCCACCCACGGACGGAAAACCGCTAACATTGAAGACGACATCTACCGCATCAAAATCCTAAGTGGGCTGCCGAGTGTACCAAAGCAGACCGTAGACCTCCACAGATGACCGCAGGACCAATCAGGAACTGAATTTCACGTCGAAGGGACACGAGGCTGTGGTTATCACGCCAGGAGCAAATAACCACCGAGCGAGCGAAGCAGTAAAATCCAAACCAAGCGGAAAAGTTTCATGCGCCAAACAGCCCGGACTCCACCAGTTTACACCACAAATGCCGCCGGTGTAAGAAGGCACGCATTAATCAAGCCCAGAACCACGTCACGCCTATAAAGTCACCGAGCAAAGCCAAGAGACGACAATCAGCCATAAGAACCCAACGTACTTATCGGGTCCGGAAACCGATCAATCGCGAGGAAGCAGCTCCACAAAGCTTTTCCCCATTCACACCAAATGAACCATGACCAACACAGCATCTAGAGCTGGGCGGAGCAGAGTTCTTACAACCAATTCTTGAGCTGAAACTAAAAGCATCTATGTAAAGACATAACCGAGTAGTGCGAGAGAACGGATGTATAGTGTGGGGCGCAGAGTACCAGCGAGCGACAACGAAGGAGACTCTCCCTAGCCAAGGGAGGCACATAAATCGTATGACAACATGGACGTAAAAACCCGATCCTCGACCTATTCACGAGTACGTAGACTCAAAAGAGCAAGAACCCTCGCTCCCATATACTAGTTGTCACTGTCGGAACTCCTACTACCAGCTCTTGGTCCGTCGAGAAGAAAAGTTTGTAGACCGAGGCAGTTAAAGAGGACTCGCCAGACTAGCAGACCTATCAGAGCAGAAGAACACATTACGCCGCGGATGTCAAGCAAGTGGGGACGGCCGTTAAGCGCCACG
    CGCAA[AG]TCGCATTA[CG][ACT]CTG[AGT]TCCTGA[ACG][AC]GCCCC[ACT]CC[CT][AG]CGGACG[GT]AA[AT][AC]CC[GT]CTAA[CG][AC]TTG[ACG][AT]G[AC]C[GT]ACAT[CGT]TACCG[CG]ATCA[AT][AC]AT[CG]CTA[ACT]GTG[GT]GCTG[CG]CGAGTGT[AT][CG]C[AT][AC]AGCA[GT]AC[CT]GT[AG]GA[CT]CTC[CG]AC[ACG]GATGA[CT]CG[CG]AGGA[CT]CAATCAGG[AG]ACTGA[ACG]TTTC[AG][CG]GTCGA[AG]G[GT]GACACG[AG]GG[CG]T[GT]TGGTT[AC]T[CG][AT]C[GT]CC[AG]GGAGC[AC][AGT]AT[AT]AC[CG]AC[CT]GAGCGA[GT]CGA[AGT]GCAGT[AGT][AC]AATCCA[AG][AT]C[CG][AG][AT]GCGGA[AC][AG][AT]GTTT[CT]ATG[CG]GC[CT]A[AC][ACG]C[AG]GCCCGGACTCCAC[CG][AT]GTTTAC[AC]C[CGT][AT][CG][ACT]A[AT]TGCCGCCGGT[GT]TA[AC][GT][AC][AGT][GT]GCA[CT][GT]CATTA[AC]TCA[AC]GC[CG]CA[GT]AA[CG]C[AG][CG]GTC[AC][CG]GCCTATA[AC][ACT][GT]TC[ACT]CCGA[GT]CA[ACT][AC]G[CG][CT]A[ACT]GAGACG[AC][CT][AGT][AT]T[CG]A[GT][CT]C[AG]T[AGT]AGA[AGT][CT]CC[AGT]AC[GT]TA[CGT]TT[AGT]TCG[GT]GTCC[GT][GT][AC][AG][AG]C[CG]GATCAAT[CT]GCGAGGAAGCAGCTC[CT][AT][CG][AC][AGT]AGCTTTT[CG][CT]C[CG]ATTCACACCAAAT[GT][ACG][AC]C[CT][AGT]TG[AC]C[CT][AT]A[CG]AC[ACG]GCAT[CT]T[AG]GAGCTGG[GT][CG][GT]G[AG]G[CG][AT]GA[GT]TTCTTAC[ACG][AGT]CC[ACT][AC]TTCTTG[AG]GCTG[AC][AT]A[CT]T[ACG]A[AT]AGCAT[CGT]TATGT[AT]A[ACG][GT][AC][CGT]AT[AT]ACCG[ACG]GT[AG]GTG[CGT]G[AC]G[AGT][GT][ACT][AT]C[GT][GT][AGT]TGT[ACG]T[AC][GT]T[GT]TGGGGCGCAGAGTA[CT]C[AG]GC[GT]AGCG[AC]CAACGA[AGT]GGAG[ACG]CTCTC[CG]CT[AGT]G[CGT]CA[AG]GGG[AC]G[GT]C[ACG][CG]AT[AT]AATCGT[ACT]T[GT]A[CT]AAC[AC]TGGACGT[ACG][AT][ACT]AA[CG]CCGATCCT[CG]GACCT[AG]TTC[AC]CG[AG]GT[AG]CGTA[GT]ACT[CT]AAAAGAGC[AC]AGAAC[CT][CT]TCG[CT]TCCCATATA[CG]TA[GT]TT[GT]TC[AC][CT]TGT[CG]G[GT][AC][AT]CT[CGT]CTACTACC[AG]GCT[CG]TT[GT]GTCCGT[CT]GAG[AG]AGAAA[AGT]GTTTGTAG[AG][CGT][CT]G[AT]GGCA[GT]TT[AG][AC][ACG]GAGGACTCG[CT]C[ACT]G[ACG][CG]T[AC]GCAG[AC]CCT[AC]T[CG][AT][GT]AG[CT][ACT]GAAG[AG]A[CG][AT]CATTACGCC[GT]CG[GT]ATGTCA[ACG][GT]C[AG][AG][GT]TGGGG[AGT]CGG[CG]CGTTAA[GT][CT]GCCACG
    A: 2 0 2 4 5 4 0 2 3 3 6 0 3 4 2 3 2 2 2 3 2 0 1 3 3 5 3 4 2 2 3 1 3 3 2 1 2 3 1 1 1 6 2 0 2 4 4 3 3 1 1 1 1 1 4 4 2 3 1 2 2 3 3 1 3 1 2 4 2 4 0 1 3 4 3 2 3 2 5 1 3 6 3 3 5 2 2 3 1 4 3 1 3 2 2 1 3 1 3 2 3 1 5 1 0 2 3 4 2 1 3 4 5 2 1 6 1 5 1 1 3 2 3 1 5 2 4 3 3 2 7 2 3 1 5 2 1 4 2 2 2 1 6 1 0 4 2 3 4 6 0 4 4 2 2 3 4 2 1 1 4 3 1 2 3 2 3 1 3 1 2 2 4 4 1 0 3 4 4 5 2 2 4 1 3 2 2 2 3 0 0 1 3 3 1 2 3 2 1 1 1 3 1 2 5 3 1 3 3 5 2 4 4 3 2 4 2 2 3 4 2 1 1 6 2 2 1 5 3 2 1 6 1 2 3 3 4 5 3 1 2 4 3 4 0 2 3 4 3 1 3 1 5 3 3 3 1 1 1 3 2 4 1 2 2 2 1 2 6 3 3 1 3 2 2 1 3 0 1 4 2 2 2 2 4 3 2 3 3 1 2 2 7 3 3 0 1 3 1 3 4 3 3 2 3 2 2 1 2 3 1 1 2 1 4 3 2 4 3 2 2 3 5 2 2 1 4 2 0 5 3 1 2 5 3 2 1 1 1 4 0 4 4 1 1 3 2 1 1 1 3 1 2 1 2 3 4 1 4 4 3 2 2 3 3 1 0 1 4 2 2 5 3 4 2 2 2 7 3 2 4 2 4 4 2 3 2 3 3 1 2 4 1 2 2 3 2 3 4 2 4 3 2 4 3 3 5 2 2 0 6 1 2 2 3 2 0 2 1 1 1 2 3 2 2 3 3 3 2 2 2 5 1 2 5 5 2 2 1 3 0 4 1 2 5 5 3 0 5 3 1 2 2 2 3 0 3 3 6 2 3 2 3 1 2 2 1 1 2 5 2 1 3 5 2 4 1 1 4 5 4 1 2 3 3 3 2 3 4 1 3 3 0 3 4 2 4 2 3 2 3 4 1 2 2 3 2 4 4 1 1 1 3 2 2 2 2 3 3 2 3 1 4 2 1 1 2 1 2 4 2 3 3 2 2 3 3 1 3 1 2 1 3 3 1 2 2 1 3 3 5 2 1 3 4 3 6 2 1 5 1 1 3 5 2 3 2 3 4 3 2 3 1 4 3 3 5 2 2 2 3 2 1 4 3 3 3 1 2 3 2 3 2 3 3 1 0 2 3 2 3 2 3 1 3 2 2 2 3 0 1 3 0 2 3 2 5 1 4 1 2 4 0 1 3 2 0 2 4 2 2 2 4 1 4 4 2 3 5 3 0 3 4 1 3 3 2 2 3 2 2 2 1 3 1 1 0 4 3 0 1 1 4 2 2 1 3 2 5 3 4 4 4 3 1 3 2 3 1 2 5 2 5 5 3 3 2 2 3 4 1 2 3 3 3 3 4 4 2 2 3 0 4 2 0 2 2 2 1 4 2 2 3 3 3 1 3 3 1 3 3 2 2 3 2 1 3 4 2 5 2 1 2 4 4 4 4 3 4 3 3 4 4 3 6 4 2 2 2 2 3 2 2 2 3 1 1 4 1 4 1 4 2 3 5 2 1 1 2 2 2 3 2 1 2 3 2 3 2 3 3 1 1 1 1 1 5 1 1 6 2 0 5 1 3 2 2 2 3 1 3 1 3 2 3 0 2 3 5 0 4 4 3 6 5 5 3 1 2 2 3 3 3 5 4 3 1 2 1 4 1 1 1 4 2 1 2 3 3 3 2 5 2 2 4 2 2 3 2 0 3 3 0 3 2 3 3 1 2 4 4 4 1 0 3 3 2 1 3 2 4 2 2 3 1 4 5 3 4 5 0 3 3 4 3 2 4 1 2 0 2 2 3 0 1 5 3 1 2 2 4 3 1 2 3 3 0 2 1 1 2 1 3 1 0 1 2 1 1 1 3 6 4 2 1 2 3 1 5 3 1
    C: 4 3 4 1 0 1 2 4 1 4 2 2 1 1 3 3 4 1 3 1 3 4 5 1 0 2 3 4 2 4 5 4 6 3 5 5 3 2 4 3 2 3 4 3 2 3 3 2 3 4 4 1 4 3 1 0 3 3 3 2 0 3 2 3 3 5 2 3 6 1 2 3 1 2 4 4 1 3 2 1 4 2 2 3 0 1 3 4 3 3 3 2 1 1 0 1 4 1 3 3 4 2 2 3 1 1 0 1 3 4 2 4 3 0 5 2 1 1 5 4 1 2 2 2 2 3 5 3 4 3 1 5 3 2 3 0 3 3 3 4 0 4 1 2 3 1 3 6 2 0 3 5 2 3 1 2 1 4 3 3 2 3 3 1 3 4 2 4 1 2 6 1 3 2 2 2 1 3 5 3 4 3 1 3 0 3 2 2 1 2 3 3 0 3 3 3 2 4 1 5 5 2 2 2 1 2 5 3 1 3 3 1 1 4 3 2 5 3 3 3 2 4 2 1 0 4 2 0 1 3 4 3 2 1 1 3 3 2 0 5 5 1 2 1 6 3 2 2 1 5 1 2 1 3 2 2 2 3 1 0 3 1 1 3 4 2 4 3 0 3 3 4 2 1 5 4 4 3 3 2 6 1 5 4 2 4 3 2 1 3 1 1 1 4 3 5 3 2 4 3 2 2 1 2 4 4 2 6 4 1 3 2 2 2 0 3 2 4 1 2 2 5 1 3 2 4 2 3 3 3 3 3 5 3 3 1 4 4 4 1 2 2 2 4 4 2 3 3 2 5 3 4 3 4 4 2 3 3 3 4 3 2 0 4 3 5 6 3 1 2 7 1 3 4 1 3 3 2 3 2 1 3 2 5 3 3 3 1 2 3 3 2 1 3 4 2 3 1 1 0 2 1 3 5 4 1 1 4 2 3 2 3 1 1 1 1 6 1 1 3 2 5 4 2 2 3 2 2 4 3 2 3 1 6 2 2 2 3 1 4 3 3 2 3 1 1 1 5 0 1 5 3 5 3 2 4 3 1 1 3 6 1 0 3 0 3 4 4 3 2 1 3 4 0 5 0 4 5 1 2 1 3 2 3 3 4 3 1 1 2 3 4 4 2 2 3 2 6 3 1 4 1 1 3 2 2 2 3 0 4 2 3 1 2 4 2 1 2 1 4 2 3 3 2 1 3 4 1 0 3 4 3 1 5 5 3 3 3 2 5 1 3 0 2 3 6 2 0 3 2 1 3 2 3 3 2 0 1 5 2 1 3 0 2 2 1 3 2 0 3 2 3 3 2 2 2 2 5 5 2 3 1 0 2 0 0 1 3 2 3 1 1 2 3 2 4 2 2 1 2 2 1 3 4 3 0 1 2 1 2 2 0 1 4 2 5 1 3 2 3 1 2 4 4 2 1 5 2 3 3 4 1 4 6 1 2 4 2 4 1 1 1 2 2 3 4 2 4 2 5 3 5 2 1 2 3 5 3 2 1 2 2 4 2 2 5 3 3 0 2 1 2 2 2 4 0 2 3 3 2 2 4 0 2 4 3 3 2 1 3 5 2 2 3 2 3 1 3 3 4 6 1 1 2 4 4 3 3 3 2 5 5 2 2 2 1 4 3 5 0 2 1 2 2 5 0 0 2 2 1 5 1 3 1 2 3 2 1 0 2 4 4 2 1 1 1 5 3 4 2 4 1 3 1 5 4 4 3 2 2 2 3 3 2 2 2 3 1 0 2 5 3 3 3 1 1 3 1 2 3 2 4 3 3 5 2 2 6 2 0 4 4 0 0 5 1 3 3 1 1 2 1 4 6 2 3 3 3 1 3 0 2 2 1 3 3 1 2 1 2 1 0 0 3 1 2 3 3 2 0 2 2 4 1 2 2 2 2 3 3 3 1 2 3 3 6 0 5 2 4 5 3 3 3 3 1 3 2 4 0 1 4 4 6 2 3 1 4 2 2 1 3 3 3 2 3 1 2 0 3 4 2 4 2 2 3 2 5 2 6 4 2 4 1 1 4 0 3 1 4 1 3 1 6 2 2 2 2 2 3 3 2 1 4 3 1 3 4 2 2 0 0 0 0 4 3 6 4 2 4 1
    G: 3 5 2 3 1 4 2 1 4 2 2 3 1 3 3 1 2 2 5 3 1 3 4 2 4 1 3 1 4 3 2 3 1 1 1 2 2 3 3 4 4 1 1 5 3 1 2 2 2 2 3 4 2 2 3 3 3 2 1 1 5 3 2 5 2 3 3 2 0 2 3 3 1 3 2 2 5 3 1 3 0 2 2 2 3 3 3 1 1 2 1 5 1 5 4 5 1 3 4 3 2 4 1 5 1 4 3 1 3 3 2 0 0 7 1 1 4 3 3 1 4 2 3 5 1 2 0 0 1 3 1 1 3 4 1 3 4 2 2 3 5 4 0 5 7 2 2 0 1 1 3 1 1 4 5 3 2 3 2 4 1 3 2 2 0 1 3 4 4 3 2 5 2 4 5 4 4 3 1 1 2 4 4 4 6 3 1 3 2 6 4 2 2 2 2 3 2 1 4 1 2 3 6 5 2 4 1 2 3 0 1 1 2 1 3 2 2 2 4 2 5 3 6 1 4 2 4 3 3 5 2 0 4 1 3 2 1 1 3 1 3 3 3 1 1 3 3 0 4 3 6 5 2 2 3 2 4 1 3 3 2 2 3 4 4 5 3 2 1 2 3 2 3 5 0 2 0 4 4 2 1 3 1 3 3 0 3 2 5 2 2 2 1 1 2 3 3 2 4 1 2 2 2 5 0 2 5 2 1 4 4 1 3 3 3 2 3 1 3 3 4 1 3 2 3 3 3 0 2 1 2 2 0 1 2 4 3 4 3 2 4 2 1 4 3 3 3 5 3 1 2 4 4 2 3 1 2 1 1 0 1 3 2 2 1 1 2 4 3 3 1 3 1 1 4 3 2 0 1 4 2 4 3 1 4 2 2 3 2 2 3 2 4 2 2 3 1 3 3 5 3 3 2 0 1 3 2 3 3 3 2 3 3 2 3 3 0 4 4 5 3 2 2 3 3 2 3 3 2 3 5 2 3 0 1 1 2 2 6 2 4 2 4 4 2 2 4 2 4 4 3 1 2 2 2 4 2 3 1 4 1 2 3 2 3 3 1 3 3 1 3 2 2 3 3 3 3 2 3 3 3 2 3 3 2 1 2 3 0 6 2 2 2 2 3 3 1 0 3 4 2 3 3 2 2 3 5 2 5 3 2 4 4 3 4 3 5 3 4 4 2 4 2 3 3 2 2 3 3 2 1 3 3 2 1 1 2 1 1 4 3 2 6 3 5 1 2 5 2 2 2 2 3 3 1 2 1 4 2 2 3 3 3 2 1 4 1 2 3 3 3 2 3 2 1 2 1 0 2 4 3 5 4 4 4 3 4 3 4 2 5 3 3 1 2 3 4 3 3 2 4 2 3 0 2 4 3 3 1 6 4 4 5 1 4 2 0 4 1 5 2 2 2 2 3 6 3 3 3 4 2 4 1 1 2 1 3 4 1 3 6 4 3 5 3 1 1 3 1 2 3 1 2 3 4 3 2 2 3 5 5 4 2 5 3 2 3 3 3 1 1 3 3 1 2 4 2 1 2 3 2 0 4 2 2 2 1 5 5 2 3 5 1 3 2 1 3 2 3 1 0 5 2 2 3 2 0 3 4 2 2 2 1 3 1 3 3 2 3 5 3 4 0 3 1 6 3 2 3 3 0 2 2 3 3 1 2 4 3 5 2 1 2 4 1 2 1 2 0 1 1 4 2 3 1 3 3 2 3 2 3 1 3 1 2 3 2 1 4 1 1 2 2 2 6 1 3 4 3 2 2 2 1 3 2 1 2 2 2 2 3 3 5 6 1 3 3 0 2 4 4 2 2 0 4 3 2 4 2 4 4 2 4 1 1 0 3 4 2 1 1 5 3 0 5 3 3 2 4 2 5 5 2 2 3 3 2 3 2 3 4 2 5 4 2 0 3 2 5 2 0 1 5 3 3 2 2 4 1 3 5 0 3 0 1 2 3 4 2 3 3 4 2 1 4 2 1 4 4 1 4 2 2 1 0 1 3 2 5 1 2 3 1 5 4 0 2 4 3 1 2 3 4 2 3 3 4 1 5 5 4 4 3 2 4 5 3 3 4 2 2 3 3 4 1 4 0 2 1 2 7
    T: 1 2 2 2 4 1 6 3 2 1 0 5 5 2 2 3 2 5 0 3 4 3 0 4 3 2 1 1 2 1 0 2 0 3 2 2 3 2 2 2 3 0 3 2 3 2 1 3 2 3 2 4 3 4 2 3 2 2 5 5 3 1 3 1 2 1 3 1 2 3 5 3 5 1 1 2 1 2 2 5 3 0 3 2 2 4 2 2 5 1 3 2 5 2 4 3 2 5 0 2 1 3 2 1 8 3 4 4 2 2 3 2 2 1 3 1 4 1 1 4 2 4 2 2 2 3 1 4 2 2 1 2 1 3 1 5 2 1 3 1 3 1 3 2 0 3 3 1 3 3 4 0 3 1 2 2 3 1 4 2 3 1 4 5 4 3 2 1 2 4 0 2 1 0 2 4 2 0 0 1 2 1 1 2 1 2 5 3 4 2 3 4 5 2 4 2 3 3 4 3 2 2 1 1 2 1 3 2 3 2 4 4 3 2 2 2 1 3 0 1 1 2 1 2 4 2 3 2 3 0 3 1 3 6 3 2 2 2 4 3 0 2 2 4 3 2 2 4 2 1 0 2 2 2 2 3 3 5 5 4 3 3 5 1 0 1 2 3 3 2 1 3 2 2 3 3 3 3 2 2 1 4 2 1 1 3 2 3 1 4 5 5 1 2 2 2 3 3 1 3 2 3 4 1 3 2 1 1 3 2 2 6 3 4 3 2 3 1 3 3 2 1 1 3 3 2 1 5 5 1 2 4 3 1 2 3 2 1 2 3 4 2 3 1 2 2 2 1 4 3 2 1 1 3 1 4 1 5 2 2 3 3 6 1 3 3 2 2 2 3 0 1 3 1 3 2 3 1 3 2 3 1 1 0 1 2 3 3 3 4 2 2 4 3 2 2 4 3 2 3 1 3 3 1 2 3 2 1 3 4 0 3 4 5 3 4 4 3 4 1 4 1 1 3 3 2 2 2 2 2 1 0 5 2 2 2 4 3 2 1 3 1 3 1 2 2 2 3 1 2 1 4 1 3 3 2 2 3 2 1 0 5 4 4 5 2 4 2 2 2 4 4 1 2 0 3 2 2 2 0 2 4 3 1 2 2 3 3 5 1 2 1 4 3 1 2 3 2 1 3 1 2 5 3 4 2 1 1 1 2 5 2 2 3 0 3 2 2 2 0 3 2 1 3 5 4 2 5 5 1 3 1 3 1 2 3 2 5 4 0 4 4 1 2 1 1 4 4 2 3 2 3 4 1 2 3 3 3 2 1 5 3 4 1 5 2 4 3 3 1 3 2 3 2 4 3 2 3 1 2 1 2 5 0 3 4 2 3 2 2 2 3 3 3 3 2 4 3 3 4 1 5 1 5 2 4 4 3 5 2 3 3 4 3 1 1 4 2 3 1 5 2 4 3 2 1 2 3 0 1 2 3 1 2 3 3 1 1 0 3 3 2 1 2 1 2 5 1 4 1 2 2 5 3 3 3 3 1 2 4 2 3 0 1 3 2 1 2 2 4 4 1 1 4 3 3 4 3 4 3 1 4 1 1 1 2 4 1 1 1 1 1 4 1 3 3 2 1 2 3 1 4 3 4 3 2 5 2 2 2 1 1 4 2 4 5 0 2 1 2 2 3 6 2 2 3 4 2 3 1 3 6 3 2 1 2 2 2 3 0 1 1 2 2 2 3 2 3 4 5 2 3 3 4 1 2 2 1 4 2 4 2 2 4 1 3 4 7 4 5 2 2 3 4 1 5 2 2 3 2 3 3 5 3 2 6 1 1 5 2 1 3 0 3 1 4 2 5 4 4 1 6 1 2 1 4 3 0 2 3 2 2 1 2 1 2 3 3 5 5 5 2 4 2 0 2 3 3 3 4 2 2 3 3 3 4 4 2 2 1 1 2 1 1 1 2 5 0 1 4 2 3 2 1 2 4 2 3 3 3 0 2 2 4 4 2 4 1 3 3 2 1 3 3 3 1 3 1 2 1 2 3 1 3 5 4 1 2 1 3 2 3 2 4 4 1 5 2 4 3 3 1 4 0 2 2 4 5 2 1 1 3 3 3 3 3 2 2 3 5 5 1 3 4 4 1 1 3 2 1 1


# Wabbit Season
    
In “Rabbits and Recurrence Relations”, we mentioned the disaster caused by introducing European rabbits into Australia. By the turn of the 20th Century, the situation was so out of control that the creatures could not be killed fast enough to slow their spread.

The conclusion? Build a fence! The fence, intended to preserve the sanctity of Western Australia, was completed in 1907 after undergoing revisions to push it back as the bunnies pushed their frontier ever westward. If it sounds like a crazy plan, the Australians at the time seem to have concurred.

By 1950, Australian rabbits numbered 600 million, causing the government to decide to release a virus (called myxoma) into the wild, which cut down the rabbits until they acquired resistance. In a final Hollywood twist, another experimental rabbit virus escaped in 1991, and some resistance has already been observed.

The bunnies will not be stopped, but they don't live forever, and so in this problem, our aim is to expand Fibonacci's rabbit population model to allow for mortal rabbits.

**Problem**

Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation Fn=Fn−1+Fn−2 and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed number of months, for exemple in which rabbits live for three months (meaning that they reproduce only twice before dying).

Given: Positive integers n≤100 and m≤20.

Return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.

Sample Dataset

6 3

Sample Output

4



```python
def MortalFibonacci(n, m):
    living = [1, 1]
    for i in range(2, n):
        # first reproduction
        tmp = living[i - 1] + living[i - 2]
        # then death
        if i == m:
            tmp = tmp - 1
        if i > m:
            tmp = tmp - living[i - m - 1]
        living.append(tmp)
    return living[-1]
```


```python
# months/generations
n = 6
# survival time
m = 3

print(MortalFibonacci(n, m))
```

    4



```python
# months/generations
n = 84
# survival time
m = 16

print(MortalFibonacci(n, m))
```

    158240333533794817


# A Brief Introduction to Graph Theory

Networks arise everywhere in the practical world, especially in biology. Networks are prevalent in popular applications such as modeling the spread of disease, but the extent of network applications spreads far beyond popular science. Our first question asks how to computationally model a network without actually needing to render a picture of the network.

First, some terminology: graph is the technical term for a network; a graph is made up of hubs called nodes (or vertices), pairs of which are connected via segments/curves called edges. If an edge connects nodes v and w, then it is denoted by v,w (or equivalently w,v).

- an edge v,w is incident to nodes v and w; we say that v and w are adjacent to each other;
- the degree of v is the number of edges incident to it;
- a walk is an ordered collection of edges for which the ending node of one edge is the starting node of the next (e.g., {v1,v2}, {v2,v3}, {v3,v4}, etc.);
- a path is a walk in which every node appears in at most two edges;
- path length is the number of edges in the path;
- a cycle is a path whose final node is equal to its first node (so that every node is incident to exactly two edges in the cycle); 
- and the distance between two vertices is the length of the shortest path connecting them.

Graph theory is the abstract mathematical study of graphs and their properties.

**Problem**

A graph whose nodes have all been labeled can be represented by an adjacency list, in which each row of the list contains the two node labels corresponding to a unique edge.

A directed graph (or digraph) is a graph containing directed edges, each of which has an orientation. That is, a directed edge is represented by an arrow instead of a line segment; the starting and ending nodes of an edge form its tail and head, respectively. The directed edge with tail v and head w is represented by (v,w) (but not by (w,v)). A directed loop is a directed edge of the form (v,v).

For a collection of strings and a positive integer k, the overlap graph for the strings is a directed graph Ok in which each string is represented by a node, and string s is connected to string t with a directed edge when there is a length k suffix of s that matches a length k prefix of t, as long as s≠t; we demand s≠t to prevent directed loops in the overlap graph (although directed cycles may be present).

Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.

Return: The adjacency list corresponding to O3. You may return edges in any order.

Sample Dataset

>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG

Sample Output

Rosalind_0498 Rosalind_2391
Rosalind_0498 Rosalind_0442
Rosalind_2391 Rosalind_2323



```python
def overlap_graph(fasta, k):
    # Initialize an empty list to store the edges in the overlap graph
    edges = []
    # Iterate over the DNA strings
    for label1, dna1 in fasta.items():
        # Compare the string with every other string
        for label2, dna2 in fasta.items():
            # If the strings are different and the end of the first string matches the beginning of the second string
            if label1 != label2 and dna1[-k:] == dna2[:k]:
                # Add an edge from the first string to the second string to the list of edges
                edges.append(label1 + " " + label2)
    return edges
```


```python
filename = 'rosalind_grph.txt'
fasta = read_fasta(filename)
k = 3
for i in overlap_graph(fasta, k):
    print(i)

```

    Rosalind_8634 Rosalind_3905
    Rosalind_1760 Rosalind_7686
    Rosalind_1760 Rosalind_0110
    Rosalind_1760 Rosalind_1200
    Rosalind_2783 Rosalind_8499
    Rosalind_2783 Rosalind_8392
    Rosalind_2783 Rosalind_0807
    Rosalind_6549 Rosalind_7913
    Rosalind_6549 Rosalind_1664
    Rosalind_6108 Rosalind_2231
    Rosalind_6108 Rosalind_1408
    Rosalind_6108 Rosalind_9717
    Rosalind_6282 Rosalind_1222
    Rosalind_4474 Rosalind_4631
    Rosalind_4474 Rosalind_4070
    Rosalind_4474 Rosalind_6434
    Rosalind_4474 Rosalind_6477
    Rosalind_3966 Rosalind_9499
    Rosalind_4631 Rosalind_3438
    Rosalind_2552 Rosalind_5174
    Rosalind_2552 Rosalind_0054
    Rosalind_7227 Rosalind_2552
    Rosalind_7227 Rosalind_0594
    Rosalind_7227 Rosalind_9979
    Rosalind_3350 Rosalind_3438
    Rosalind_9986 Rosalind_6282
    Rosalind_9986 Rosalind_3394
    Rosalind_3905 Rosalind_3147
    Rosalind_3905 Rosalind_8253
    Rosalind_1894 Rosalind_2723
    Rosalind_0182 Rosalind_2723
    Rosalind_6784 Rosalind_4474
    Rosalind_5174 Rosalind_2600
    Rosalind_5174 Rosalind_7055
    Rosalind_4143 Rosalind_9696
    Rosalind_0647 Rosalind_6549
    Rosalind_0647 Rosalind_6562
    Rosalind_0647 Rosalind_8658
    Rosalind_0041 Rosalind_9199
    Rosalind_0041 Rosalind_6858
    Rosalind_0041 Rosalind_8039
    Rosalind_3396 Rosalind_3198
    Rosalind_3396 Rosalind_3404
    Rosalind_3394 Rosalind_6985
    Rosalind_3394 Rosalind_6082
    Rosalind_3072 Rosalind_0484
    Rosalind_3072 Rosalind_4963
    Rosalind_3072 Rosalind_0232
    Rosalind_4909 Rosalind_9199
    Rosalind_4909 Rosalind_6858
    Rosalind_4909 Rosalind_8039
    Rosalind_1408 Rosalind_2723
    Rosalind_9849 Rosalind_4474
    Rosalind_3613 Rosalind_3350
    Rosalind_3613 Rosalind_7687
    Rosalind_2615 Rosalind_0182
    Rosalind_2615 Rosalind_1257
    Rosalind_2615 Rosalind_1634
    Rosalind_2615 Rosalind_5009
    Rosalind_4687 Rosalind_0484
    Rosalind_4687 Rosalind_4963
    Rosalind_4687 Rosalind_0232
    Rosalind_1413 Rosalind_3350
    Rosalind_1413 Rosalind_7687
    Rosalind_7686 Rosalind_5174
    Rosalind_7686 Rosalind_0054
    Rosalind_0594 Rosalind_2600
    Rosalind_0594 Rosalind_7055
    Rosalind_5659 Rosalind_9696
    Rosalind_4070 Rosalind_1894
    Rosalind_4070 Rosalind_6489
    Rosalind_2600 Rosalind_6190
    Rosalind_2600 Rosalind_5135
    Rosalind_0054 Rosalind_9849
    Rosalind_0054 Rosalind_3909
    Rosalind_7687 Rosalind_3905
    Rosalind_4748 Rosalind_7913
    Rosalind_4748 Rosalind_1664
    Rosalind_6434 Rosalind_3905
    Rosalind_0319 Rosalind_2723
    Rosalind_3198 Rosalind_0182
    Rosalind_3198 Rosalind_1257
    Rosalind_3198 Rosalind_1634
    Rosalind_3198 Rosalind_5009
    Rosalind_1253 Rosalind_9199
    Rosalind_1253 Rosalind_6858
    Rosalind_1253 Rosalind_8039
    Rosalind_3404 Rosalind_0613
    Rosalind_3404 Rosalind_0319
    Rosalind_3382 Rosalind_4143
    Rosalind_3382 Rosalind_1413
    Rosalind_7055 Rosalind_8499
    Rosalind_7055 Rosalind_8392
    Rosalind_7055 Rosalind_0807
    Rosalind_2084 Rosalind_3905
    Rosalind_0372 Rosalind_1894
    Rosalind_0372 Rosalind_6489
    Rosalind_0484 Rosalind_0647
    Rosalind_6985 Rosalind_4474
    Rosalind_0110 Rosalind_6190
    Rosalind_0110 Rosalind_5135
    Rosalind_7913 Rosalind_2552
    Rosalind_7913 Rosalind_0594
    Rosalind_7913 Rosalind_9979
    Rosalind_6562 Rosalind_0041
    Rosalind_6190 Rosalind_3905
    Rosalind_9955 Rosalind_6282
    Rosalind_9955 Rosalind_3394
    Rosalind_3565 Rosalind_8499
    Rosalind_3565 Rosalind_8392
    Rosalind_3565 Rosalind_0807
    Rosalind_6788 Rosalind_1253
    Rosalind_6082 Rosalind_1347
    Rosalind_5135 Rosalind_2783
    Rosalind_5135 Rosalind_3072
    Rosalind_4589 Rosalind_7102
    Rosalind_8392 Rosalind_4589
    Rosalind_6477 Rosalind_6108
    Rosalind_6477 Rosalind_4524
    Rosalind_6477 Rosalind_0372
    Rosalind_7102 Rosalind_1894
    Rosalind_7102 Rosalind_6489
    Rosalind_4963 Rosalind_2231
    Rosalind_4963 Rosalind_1408
    Rosalind_4963 Rosalind_9717
    Rosalind_1347 Rosalind_9849
    Rosalind_1347 Rosalind_3909
    Rosalind_1222 Rosalind_1347
    Rosalind_8658 Rosalind_2783
    Rosalind_8658 Rosalind_3072
    Rosalind_6919 Rosalind_6784
    Rosalind_1664 Rosalind_1760
    Rosalind_1664 Rosalind_4687
    Rosalind_0906 Rosalind_7102
    Rosalind_7371 Rosalind_4631
    Rosalind_7371 Rosalind_4070
    Rosalind_7371 Rosalind_6434
    Rosalind_7371 Rosalind_6477
    Rosalind_3909 Rosalind_3396
    Rosalind_3909 Rosalind_9955
    Rosalind_3909 Rosalind_6788
    Rosalind_0232 Rosalind_6190
    Rosalind_0232 Rosalind_5135
    Rosalind_1634 Rosalind_9986
    Rosalind_9717 Rosalind_2723
    Rosalind_6489 Rosalind_7686
    Rosalind_6489 Rosalind_0110
    Rosalind_6489 Rosalind_1200
    Rosalind_6858 Rosalind_3396
    Rosalind_6858 Rosalind_9955
    Rosalind_6858 Rosalind_6788
    Rosalind_0807 Rosalind_1222
    Rosalind_3438 Rosalind_7227
    Rosalind_3438 Rosalind_3613
    Rosalind_1200 Rosalind_3905
    Rosalind_8039 Rosalind_0695
    Rosalind_8039 Rosalind_2615
    Rosalind_5009 Rosalind_4474
    Rosalind_8253 Rosalind_0695
    Rosalind_8253 Rosalind_2615


```{note}
This is a note
```

:::{important}
We believe in a community-driven approach of open-source tools that are...
:::

:::{figure} https://source.unsplash.com/random/400x200?beach,ocean
:name: my-fig
:alt: Random image of the beach or ocean!

Relaxing at the beach 🏝 🌊 😎
:::

:::{list-table} This is a nice table!
:header-rows: 1
:name: example-table

* - Training
  - Validation
* - 0
  - 5
* - 13720
  - 2744
:::

```{warning}
This is a warning block.
```

```{admonition} Take this as a warning!
:class: warning

This is a warning block.
```

```{danger}
This is a danger block.
```

```{seealso} 
This is a see also block.
```

```{admonition} Algorithm: **Counting DNA Nucleotides**
:class: danger

_________________________________________________________________________________


<code>**create a dictionary** "counts" with keys "A", "C", "G", "T" and set all of them to 0     
**for each** symbol in the string "s"
    **if** symbol **exists** in the dictionary "counts"
        **increment the count for the key** "symbol" in the "counts" dictionary    
**return** a string containing the values for the keys "A", "C", "G", "T" in the "counts" dictionary</code>
_________________________________________________________________________________
```


```python

```


```python
[1
 2
 3
 4
 5
 6]
```
