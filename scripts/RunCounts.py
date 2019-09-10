from pickle import load,dump
annoMap,revMap = {},{}
for line in open("annotations/annotationDescripts.tsv"):
    rec = line.strip().split("\t")
    annoMap[rec[0]]= rec[1]
    try:revMap[rec[1]].add(rec[0])
    except:revMap[rec[1]] = set([rec[0]])
excCounter = 0
geneCounts = {}
allCounts = load(open("pickles/allCounts_genes.p","rb"))
for index,sample in enumerate(allCounts.columns):
    geneCounts[sample]={}
    sampleCounts = allCounts[sample]
    print("%i. %s" % (index+1,sample))
    for function,contigs in revMap.items():
        functCounts = sampleCounts[sampleCounts.index.isin(contigs)]
        nReads = float(functCounts.sum(axis = 0))
        try: geneCounts[sample][function] = nReads/len(functCounts)
        except:
            excCounter +=1
            geneCounts[sample][function] = 0.0
dump(geneCounts,open("pickles/pooledGeneCounts.p","wb"))
print("There were %i genes with no representative contigs" % (excCounter))