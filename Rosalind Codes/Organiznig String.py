#!/usr/bin/env python
# coding: utf-8

# In[1]:


rosalind_lexf = open('/Users/mac/Desktop/courses/Bioinformatic/rosalind_lexf.txt').read()
Dataset = "A C G T\n2"

def alpha_combination(alphabet, n, s = '', x = []):
    if n == 0:
        x.append(s)
    else:
        for c in alphabet:
            alpha_combination(alphabet, n - 1, s + c, x)
    return x

letter = rosalind_lexf.split()

def reverse(letter):
    s = ''
    for c in letter:
        s = c + s
    return s

alphabet = letter
organized_length = int(letter[-1])

for output in alpha_combination(alphabet[:(len(alphabet)-1)], organized_length):
    print(output) 


# In[ ]:




