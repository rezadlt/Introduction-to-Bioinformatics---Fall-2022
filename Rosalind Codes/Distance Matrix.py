#!/usr/bin/env python
# coding: utf-8

# In[10]:


def p_distance(s1, s2):
    assert len(s1) == len(s2)
    hamming_distance = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            hamming_distance += 1
    return hamming_distance / len(s1)

def main():
    file = open("rosalind_pdst.txt", "r")
    sequences = []
    file_list = []
    temp = 0
    for line in file:
        line = line.strip()
        if line[0] == '>':
            temp = temp +1
            if line and temp != 1:
                sequences.append(" ".join(file_list))
                file_list = []
        else:
            file_list.append(line)
    sequences.append(" ".join(file_list))
    distance_matrix = [[0 for i in range(len(sequences))] for j in range(len(sequences))]
    for x in range(len(sequences)):
        for y in range(len(sequences)):
            if x != y:
                distance_matrix[x][y] = p_distance(sequences[x], sequences[y])
    for row in distance_matrix:
        for data in row:
            print("{:.5f}".format(round(data,3)), end=' ')
        print()
    return

if __name__ == "__main__":
    main()


# In[ ]:




