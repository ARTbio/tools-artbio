#!/usr/bin/python

import sys

# python script to perform global searches and replacements of multiple words
# The script takes in input a tabular file for correspondances between searched items and item replacements
# a file with the text to replace
# and outputs the new text with replacements

# Usage:
# searchReplace.py <correspondance_table> <file_to_edit> <edited_file>

# define the search replacement method
def replace_all(text, dic):
    for i, j in dic.iteritems():
        text = text.replace(i, j)
    return text

# Create the dictionary
table = open (sys.argv[1], "r")
Dic = {}
for line in table:
    fields = line.split()
    Dic[fields[0]] = fields[1]
table.close()

# takes the input file to edit
input_file = open (sys.argv[2], "r")
input_text = input_file.read()
input_file.close()

# bind the returned text of the method to a variable
new_text = replace_all(input_text, Dic)

# outputs the new_text in a new file

Output = open (sys.argv[3], "w")
print >> Output, new_text[:-1]

