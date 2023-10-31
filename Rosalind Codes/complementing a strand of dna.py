#!/usr/bin/env python
# coding: utf-8

# In[7]:


DNA = input()
DNA = list(DNA)
for i in range(len(DNA)):
    if DNA[i] == "A":
        DNA[i] = "T" 
    elif DNA[i] == "T":
        DNA[i] = "A" 
    elif DNA[i] == "C":
        DNA[i] = "G"
    elif DNA[i] == "G":
        DNA[i] = "C"
        
def reverse(string):
    s = ''
    for c in string:
        s = c + s
    return s

print(reverse(DNA))


# In[ ]:




