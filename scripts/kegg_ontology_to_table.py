__author__ = 'vetrot'

import sys
import os

kegg_ontology = sys.argv[1]
#A09100 Metabolism
#B
#B  09101 Carbohydrate metabolism
#C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
#D      K00844  HK; hexokinase [EC:2.7.1.1]
output_tab = sys.argv[2]

A_class = ''
B_class = ''
C_class = ''
D_class = ''

fp = open(output_tab, 'w')
rec = "level1\tlevel2\tlevel3\tfunction\tid"
fp.write(rec + "\n")
for n, line in enumerate(open(kegg_ontology)):
    char = line.rstrip()[0]
    nl = ' '.join(line.rstrip()[1:].split())
    if char == 'A':
        vals = nl.split(" ")
        v =  vals[1]
        if len(vals)>2:
            for i in range(2, len(vals)):
                v = v + ' '+ vals[i]
            # print nl
            # print v
        A_class = v
    if char == 'B':
        if nl != '':
            vals = nl.split(" ")
            v =  vals[1]
            if len(vals)>2:
                for i in range(2, len(vals)):
                    v = v + ' '+ vals[i]
                #print nl
                #print v
            B_class = v
    if char == 'C':
        vals = nl.split(" ")
        v =  vals[1]
        if len(vals)>2:
            for i in range(2, len(vals)):
                v = v + ' '+ vals[i]
            #print nl
            #print v
        C_class = v
    if char == 'D':
        vals = nl.split(" ")
        v =  vals[1]
        if len(vals)>2:
            for i in range(2, len(vals)):
                v = v + ' '+ vals[i]
            #print nl
            #print v
        D_class = v
        rec = A_class+"\t"+B_class+"\t"+C_class+"\t"+D_class+"\t["+vals[0]+"]"
        fp.write(rec + "\n")
fp.close()

print ("done...")



