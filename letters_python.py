# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 16:53:59 2020

@author: Jess
"""
import pandas as pandas
inputfile = open("three_column_significance_letters_K0.txt", "r")

list_of_interesting_COGs = []
for line in inputfile:
    line = line.rstrip().split(" ")
    cog = line[0]
    hapletters = line[1]
    napletters = line[2]
    if any(x in hapletters for x in napletters) is False:
        list_of_interesting_COGs.append([cog, hapletters, napletters])
outputdf = pandas.DataFrame(list_of_interesting_COGs)
outputdf.to_csv(r'interesting_letters.csv')

inputfile.close()
