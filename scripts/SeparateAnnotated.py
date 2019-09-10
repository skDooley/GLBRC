from glob import glob
from Bio.SeqIO import parse, write
from Bio.SeqRecord import SeqRecord
annotatedRegions = open("assemblies/assemblyGeneSeqs.fa","w")
geneDescripts = open("annotations/annotationDescripts.tsv","w")
annotations = glob("annotations/prokka/final.contigs.*/*.gbk")
for fname in annotations:
	print(fname)    
	for rec in parse(fname,'genbank'):
		for feature in rec.features:
			try:
				product = feature.qualifiers['product'][0]
				if product == "hypothetical protein":continue
				subSeq = rec.seq[feature.location.start:feature.location.end]
				seqID = "%s_%i_%i" % (rec.id,feature.location.start,feature.location.end)
				newSeq = SeqRecord(subSeq,seqID,'','')
				geneDescripts.write("%s\t%s\n" % (seqID,product))
				write(newSeq,annotatedRegions,"fasta")
			except:pass
annotatedRegions.close()
