from glob import glob
from pandas import DataFrame, read_csv, Series, to_datetime
from pickle import dump, load
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from random import randint 
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord
annotatedRegions = open("assemblies/AnnotatedContigs.fa","w")
geneDescripts = open("annotations/AnnotationDescripts.bed","w")
annotations = glob("mags/final.contigs.*/*.gbk")
seqs = set()
hypoCounter = 0
for fname in annotations:
    print(fname)    
    for rec in parse(fname,'genbank'):
        for feature in rec.features:
            try:
                product = feature.qualifiers['product'][0]
                hypoCounter += int(product == "hypothetical protein")
                geneDescripts.write("%s\t%i\t%i\t%s\n" % (rec.id,feature.location.start,feature.location.end,product))
                if rec.id in seqs: continue
                rec.seq = rec.seq.upper()
                write(rec,annotatedRegions,"fasta")
                seqs.add(rec.id)
            except:pass
annotatedRegions.close()
print(len(seqs),"Seqs with hypo prots:",len(hypoCounter))