'''
Takes the KEGG mapping files and connects each level of KO mapping pipeline
'''


from bs4 import BeautifulSoup
from urllib.request import urlopen
from time import sleep

class KO_Mapper:
    def __init__(self):
        self.koToMap = {}
        self.descripts = {}
        self.mapToKO = {}
        self.levelToMap = {}
        self.ko = ""
        self.KEGG_URL = "https://www.kegg.jp/dbget-bin/www_bget?%s"

    def setItem(self, level, mapNum, descrip):
        # print(self.ko, level, mapNum, descrip)
        # Map descriptions
        self.descripts[mapNum] = descrip
        
        # Map to KO
        try: self.mapToKO[mapNum].add(self.ko)
        except: self.mapToKO[mapNum] = set([self.ko])
        
        # KO to level to mapNum
        try: self.koToMap[self.ko][level].add(mapNum)
        except: 
            try: self.koToMap[self.ko][level] = set([mapNum])
            except: self.koToMap[self.ko] = { level:set([mapNum]) }
        
        # Level to map
        try: self.levelToMap[level].add(mapNum)
        except: self.levelToMap[level] = set([mapNum])
                
    def _processKOInfo(self, koText):
        rec = koText.strip().split('\n')
        while 'KEGG Orthology' not in rec[0] and len(rec)>0: rec = rec[1:]
        if len(rec)==0:return
        for line in rec[1:]:
            level = line.count('\xa0')
            if (level == 0) or self.ko in str(line): break #Don't record any ribosome info
            info = line.strip('\xa0')
            bIndex = info.find(" ")
            self.setItem(level,info[:bIndex],info[bIndex+1:])
    
    def mapKO(self,ko):
        siteContent = BeautifulSoup(urlopen(self.KEGG_URL%(ko)).read(),features="lxml")
        self.ko = ko
        found = False
        for i, elm in enumerate(siteContent.find_all("td", {"class": "td41"})):
            if "KEGG Orthology" in str(elm.contents[0]):
                pathway = elm.find("nobr")
                self._processKOInfo(pathway.text)
                found = True
                break
        if not found: 
            for i, elm in enumerate(siteContent.find_all("td", {"class": "td40"})):
                if "KEGG Orthology" in str(elm.contents[0]):
                    pathway = elm.find("nobr")
                    self._processKOInfo(pathway.text)
                    found = True
                    break
        if not found: raise Exception("Unable to find " + self.ko)

if __name__ == '__main__':
    from os import chdir
    chdir("/mnt/research/ShadeLab/GLBRC")
    ko_s = set() 
    failed = set()
    for line in open("annotations/KEGG_tools_out/diamondAnnotations_0_KOtable.txt"): ko_s.add(line.split('\t')[2])
    print(len(ko_s))
    logFileName = "AnnotationProgress.txt"
    with open(logFileName,'w') as fh: fh.write("Progress of annotation mapping\n")

    ko_mapper = KO_Mapper()
    for i,ko in enumerate(ko_s): 
        try: ko_mapper.mapKO(ko)
        except:
            sleep(3)
            try: ko_mapper.mapKO(ko)
            except: failed.add(ko)
        if i % 100 == 0:
            with open(logFileName,'a') as fh: fh.write("%i " % (i))

    print(ko_mapper.koToMap.keys())
    print("Failed %i KOs" % (len(failed)))
    print("# of KOs:",len(ko_mapper.koToMap) )
    from pickle import dump
    dump(failed,open("pickles/koMapFails.p","wb"))
    dump(ko_mapper,open("pickles/koMap.p","wb"))

