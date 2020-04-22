import csv

def comparingFiles(mrnaRegions, tsrfile1, tsrfile2, tsrfile3, tsrfile4):
    with open(mrnaRegions, 'r') as mrna, open(tsrfile1, 'r') as tsr1, open('peaksInGenes.txt', 'w') as fileWrite, open(tsrfile2, 'r') as tsr2, open(tsrfile3, 'r') as tsr3, open(tsrfile4, 'r') as tsr4, open('tsr1_peaks.csv', 'w') as tsr1_write, open('tsr2_peaks.csv', 'w') as tsr2_write, open('tsr3_peaks.csv', 'w') as tsr3_write, open('tsr4_peaks.csv', 'w') as tsr4_write:

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
        flag = 0
        offset = 0
        len_offset1 = {}
        len_offset2 = {}
        len_offset3 = {}
        len_offset4 = {}

        for t in tsrLines1:
            parse = t.split()
            if parse[0] not in len_offset1:
                len_offset1[parse[0]] = offset
            offset += len(t)

        offset = 0
        for t in tsrLines2:
            parse = t.split()
            if parse[0] not in len_offset2:
                len_offset2[parse[0]] = offset
            offset += len(t)

        offset = 0
        for t in tsrLines3:
            parse = t.split()
            if parse[0] not in len_offset3:
                len_offset3[parse[0]] = offset
            offset += len(t)

        offset = 0
        for t in tsrLines4:
            parse = t.split()
            if parse[0] not in len_offset4:
                len_offset4[parse[0]] = offset
            offset += len(t)

        tsr1.seek(0)
        # print(tsr1.readline())

        fieldNames= ['chrom', 'ID', 'start', 'end', 'strand', 'nTSS', 'nTags', 'tsrwidth']
        writer1 = csv.DictWriter(tsr1_write, fieldnames=fieldNames)
        writer1.writeheader()
        for line in mrnaLines:
            chrom, ID, start, end, strand = line.split()
            fileWrite.write("{} {}:\n".format(chrom, ID))

            if ID not in genePeakDict:
                genePeakDict[ID] = {'+': [], '-': []}
                tsr1_peaks[ID] = {'+': [], '-': []}


            # print(chrom)
            # genePeakDict[chrom][strand].append('here')
            prev = ""
            tsr1.seek(0)
            if len_offset1.get(chrom, "") == "":
                continue
            else:
                tsr1.seek(len_offset1[chrom])
                tsrLines1 = tsr1.readlines()
            for t in tsrLines1:
                parse = t.split()


                if prev == chrom and parse[0] != chrom:
                    break

                if parse[0] != chrom:
                    continue

                if int(parse[1]) >= int(start) and int(parse[2]) <= int(end):
                    genePeakDict[ID][parse[3]].append(parse)
                    tsr1_peaks[ID][parse[3]].append(parse)
                    writer1.writerow({'chrom': chrom, 'ID': ID, 'start': start, 'end': end, 'strand': strand, 'nTSS': parse[4], 'nTags': parse[5] , 'tsrwidth': parse[7] })
                prev = parse[0]


            # fileWrite.write("\t{}:\n".format('+'))
            # for value in genePeakDict[ID]['+']:
            #     fileWrite.write("\t\t{}\n".format(value))

            # fileWrite.write("\t{}:\n".format('-'))
            # for value in genePeakDict[ID]['-']:
            #     fileWrite.write("\t\t{}\n".format(value))

        writer2 = csv.DictWriter(tsr2_write, fieldnames=fieldNames)
        writer2.writeheader()
        for line in mrnaLines:
            chrom, ID, start, end, strand = line.split()

            if tsr2_peaks.get(ID, "") == "":
                # genePeakDict[ID] = {'+': [], '-': []}
                tsr2_peaks[ID] = {'+': [], '-': []}


            # print(chrom)
            # genePeakDict[chrom][strand].append('here')
            prev = ""
            tsr2.seek(0)
            if len_offset2.get(chrom, "") == "":
                continue
            else:
                tsr2.seek(len_offset2[chrom])
                tsrLines2 = tsr2.readlines()
            for t in tsrLines2:
                parse = t.split()


                if prev == chrom and parse[0] != chrom:
                    break

                if parse[0] != chrom:
                    continue

                if int(parse[1]) >= int(start) and int(parse[2]) <= int(end):
                    genePeakDict[ID][parse[3]].append(parse)
                    tsr2_peaks[ID][parse[3]].append(parse)
                    writer2.writerow({'chrom': chrom, 'ID': ID, 'start': start, 'end': end, 'strand': strand, 'nTSS': parse[4], 'nTags': parse[5] , 'tsrwidth': parse[7] })
                prev = parse[0]

        writer3 = csv.DictWriter(tsr3_write, fieldnames=fieldNames)
        writer3.writeheader()
        for line in mrnaLines:
            chrom, ID, start, end, strand = line.split()

            if tsr3_peaks.get(ID, "") == "":
                # genePeakDict[ID] = {'+': [], '-': []}
                tsr3_peaks[ID] = {'+': [], '-': []}


            # print(chrom)
            # genePeakDict[chrom][strand].append('here')
            prev = ""
            tsr3.seek(0)
            if len_offset3.get(chrom, "") == "":
                continue
            else:
                tsr3.seek(len_offset3[chrom])
                tsrLines3 = tsr3.readlines()
            for t in tsrLines3:
                parse = t.split()


                if prev == chrom and parse[0] != chrom:
                    break

                if parse[0] != chrom:
                    continue

                if int(parse[1]) >= int(start) and int(parse[2]) <= int(end):
                    genePeakDict[ID][parse[3]].append(parse)
                    tsr3_peaks[ID][parse[3]].append(parse)
                    writer3.writerow({'chrom': chrom, 'ID': ID, 'start': start, 'end': end, 'strand': strand, 'nTSS': parse[4], 'nTags': parse[5] , 'tsrwidth': parse[7] })
                prev = parse[0]

        writer4 = csv.DictWriter(tsr4_write, fieldnames=fieldNames)
        writer4.writeheader()
        for line in mrnaLines:
            chrom, ID, start, end, strand = line.split()

            if tsr4_peaks.get(ID, "") == "":
                # genePeakDict[ID] = {'+': [], '-': []}
                tsr4_peaks[ID] = {'+': [], '-': []}


            # print(chrom)
            # genePeakDict[chrom][strand].append('here')
            prev = ""
            tsr4.seek(0)
            if len_offset4.get(chrom, "") == "":
                continue
            else:
                tsr4.seek(len_offset4[chrom])
                tsrLines4 = tsr4.readlines()
            for t in tsrLines4:
                parse = t.split()


                if prev == chrom and parse[0] != chrom:
                    break

                if parse[0] != chrom:
                    continue

                if int(parse[1]) >= int(start) and int(parse[2]) <= int(end):
                    genePeakDict[ID][parse[3]].append(parse)
                    tsr4_peaks[ID][parse[3]].append(parse)
                    writer4.writerow({'chrom': chrom, 'ID': ID, 'start': start, 'end': end, 'strand': strand, 'nTSS': parse[4], 'nTags': parse[5] , 'tsrwidth': parse[7] })
                prev = parse[0]

        return [tsr1_peaks, tsr2_peaks, tsr3_peaks, tsr4_peaks]


if __name__ == '__main__':
    mrnaRegions = 'promoterRegions.txt'
    tsrfile1 = 'TSRset-1_AML.tab'
    tsrfile2 = 'TSRset-2_B20.tab'
    tsrfile3 = 'TSRset-3_NEM.tab'
    tsrfile4 = 'TSRset-4_PC.tab'

    peaks = comparingFiles(mrnaRegions, tsrfile1, tsrfile2, tsrfile3, tsrfile4)


    # print(peaks)
