def comparingFiles(mrnaRegions, tsrfile1, tsrfile2, tsrfile3, tsrfile4):
    with open(mrnaRegions, 'r') as mrna, open(tsrfile1, 'r') as tsr1, open('peaksInGenes.txt', 'w') as fileWrite, open(tsrfile2, 'r') as tsr2, open(tsrfile3, 'r') as tsr3, open(tsrfile4, 'r') as tsr4:

        mrnaLines = mrna.readlines()
        tsrLines1 = tsr1.readlines()
        tsrLines2 = tsr2.readlines()
        tsrLines3 = tsr3.readlines()
        tsrLines4 = tsr4.readlines()
        genePeakDict = dict()
        tsr1_peaks = dict()
        tsr2_peaks = dict()
        tsr3_peaks = dict()
        tsr4_peaks = dict()

        for line in mrnaLines:
            chrom, ID, start, end, strand = line.split()
            fileWrite.write("{} {}:\n".format(chrom, ID))

            if ID not in genePeakDict:
                genePeakDict[ID] = {'+': [], '-': []}
                tsr1_peaks[ID] = {'+': [], '-': []}
                tsr2_peaks[ID] = {'+': [], '-': []}
                tsr3_peaks[ID] = {'+': [], '-': []}
                tsr4_peaks[ID] = {'+': [], '-': []}


            # genePeakDict[chrom][strand].append('here')
            for t in tsrLines1:
                parse = t.split()

                if parse[0] != chrom:
                    break

                if int(parse[1]) >= int(start) and int(parse[2]) <= int(end):
                    genePeakDict[ID][parse[3]].append(parse)
                    tsr1_peaks[ID][parse[3]].append(parse)

            for t in tsrLines2:
                parse = t.split()

                if parse[0] != chrom:
                    break

                if int(parse[1]) >= int(start) and int(parse[2]) <= int(end):
                    genePeakDict[ID][parse[3]].append(parse)
                    tsr2_peaks[ID][parse[3]].append(parse)

            for t in tsrLines3:
                parse = t.split()

                if parse[0] != chrom:
                    break

                if int(parse[1]) >= int(start) and int(parse[2]) <= int(end):
                    genePeakDict[ID][parse[3]].append(parse)
                    tsr3_peaks[ID][parse[3]].append(parse)

            for t in tsrLines4:
                parse = t.split()

                if parse[0] != chrom:
                    break

                if int(parse[1]) >= int(start) and int(parse[2]) <= int(end):
                    genePeakDict[ID][parse[3]].append(parse)
                    tsr4_peaks[ID][parse[3]].append(parse)

            fileWrite.write("\t{}:\n".format('+'))
            for value in genePeakDict[ID]['+']:
                fileWrite.write("\t\t{}\n".format(value))

            fileWrite.write("\t{}:\n".format('-'))
            for value in genePeakDict[ID]['-']:
                fileWrite.write("\t\t{}\n".format(value))


if __name__ == '__main__':
    mrnaRegions = 'promoterRegions.txt'
    tsrfile1 = 'TSRset-1_AML.tab'
    tsrfile2 = 'TSRset-2_B20.tab'
    tsrfile3 = 'TSRset-3_NEM.tab'
    tsrfile4 = 'TSRset-4_PC.tab'

    comparingFiles(mrnaRegions, tsrfile1, tsrfile2, tsrfile3, tsrfile4)
