#!/usr/bin/env python
# coding: utf-8

# In[1]:


rosalind_gc = open('/Users/mac/Desktop/courses/Bioinformatic/rosalind_gc.txt')
content = rosalind_gc.readlines()

def add_letter(content):
    x = []
    for i in range (len(content)):
        content[i] = str(content[i])
        if line.startswith('>'):
            id = line[1:]
            sequence = ''
        else:
            sequence_new = line.strip()
            sequence = sequence + sequence_new
            if i==len(content)-1 or content[i+1].startswith('>'):
                GC = 100 * (sequence.count('G') + sequence.count('C')) / len(sequence)
                # string that are 'G' or 'C'
                if GC > max_GC:
                    max_GC = gc
                    max_id = id
                    x.append(max_GC)
                    x.append(max_id)
                    
    return x
id = add_letter(content)
print(id[1], end='')
print(id[0])

