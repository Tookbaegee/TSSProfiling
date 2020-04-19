import sys

def parsemRNA():
    with open("hg19.ncbiRefSeq.gff3", "r") as gff3, open('outputNeg.gff3', 'w') as writer:
        result = []
        r = []
        s = gff3.readlines()
        i = 0
        for line in s:
            # if i > 1:
            #     p = line.split()
            #     if p[2] == 'mRNA' and p[6] == '-':
            #         r.append(line)
            #     if p[2] == 'mRNA' and p[6] == '+':
            #         result.append(line)
            # i = i  + 1
            if "mRNA" in line:
                result.append(line)


        with open("output.gff3", "w") as output:
            output.writelines(result)
        # writer.writelines(r)
def checkIntervals():
    with open("output.gff3", "r") as gff3, open('check.gff3', 'w') as writer:
        s = gff3.readlines()
        linesToKeep = []
        for line in s:
            p = line.split()
            flag = True
            for xline in s:
                t = xline.split()

                if t[6] == p[6]:
                    if p[3] > t[3] and p[4] < t[4]:
                        flag = False
                        break
                if p[0] != t[0]:
                    break
            if flag == True:
                linesToKeep.append(line)

        writer.writelines(linesToKeep)

def calcRange():
    with open("mRNARegions.gff3", "r") as mRNAs:
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
            currStrand = parsed[6]
            #print("currChr: {}\ncurrGene: {}".format(currChr, currGene))

            if (currChr, currGene) not in genes:
                genes[(currChr, currGene, currStrand)] = []
            genes[(currChr, currGene, currStrand)].append((startInd, endInd, currStrand))




        for gene in genes:
            minStart = sys.maxsize
            maxStart = 0
            for start, _ , _ in genes[gene]:
                if start < minStart:
                    minStart = start
                if start > maxStart:
                    maxStart = start
            genesInt[gene] = (minStart - 20000, maxStart)



        with open("promoterRegions.txt", "w") as promoterfile:
            for gene in genesInt:
                promoterfile.write("{}\t{}\t{}\t{}\t{}\n".format(gene[0], gene[1], genesInt[gene][0], genesInt[gene][1], gene[2]))









if __name__ == "__main__":
    calcRange()
    # parsemRNA()
    # checkIntervals()
