from Bio.SeqIO import parse, write
from glob import glob
annotatedRegions = open("assemblies/annotatedContigs.fa","w")
annotations = glob("annotations/prokka/final.contigs.*/*.gbk")
for fname in annotations:
	    print(fname)    
	    for rec in parse(fname,'genbank'):
		    for feature in rec.features:
		        if "product" in feature.qualifiers and feature.qualifiers['product'][0] != "hypothetical protein":
		            rec.description = ""
		            write(rec,annotatedRegions,"fasta")
		            break
annotatedRegions.close()
