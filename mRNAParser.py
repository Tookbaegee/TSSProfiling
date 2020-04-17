import sys
  
def parsemRNA():
    with open("hg19.ncbiRefSeq.gff3", "r") as gff3:
        result = []
        s = gff3.readlines()
        for line in s:
            if "mRNA" in line:
                result.append(line)
        
        with open("output.gff3", "w") as output:
            output.writelines(result)  

def calcRange():
    with open("output.gff3", "r") as mRNAs:
        s = mRNAs.readlines()
        prevChr = ""
        genes = dict()
        genesInt = dict()
        for line in s:
            parsed = line.split()
            currChr = parsed[0]
            startInd = int(parsed[3])
            endInd = int(parsed[4])
            currGene = parsed[8].split(';')[1][7:]
            #print("currChr: {}\ncurrGene: {}".format(currChr, currGene))
                        
            if (currChr, currGene) not in genes:
                genes[(currChr, currGene)] = []
            genes[(currChr, currGene)].append((startInd, endInd))
 
        
            

        for gene in genes:
            minStart = sys.maxsize
            maxStart = 0
            for start, _ in genes[gene]:
                if start < minStart:
                    minStart = start
                if start > maxStart:
                    maxStart = start
            genesInt[gene] = (minStart - 20000, maxStart)

            
        with open("promoterRegions.txt", "w") as promoterfile:
            for gene in genesInt:
                promoterfile.write("{}\t{}\t{}\t{}\n".format(gene[0], gene[1], genesInt[gene][0], genesInt[gene][1]))
            
        

            
              

            
        

if __name__ == "__main__":
    calcRange()