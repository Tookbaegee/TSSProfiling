'''
Requires scikit-learn, scipy, numpy, and pandas
'''

import pandas as pd
import numpy as np
import csv
from scipy.stats import wilcoxon

def normalize(ts1_csv, ts2_csv, ts3_csv, ts4_csv):
    # only use columns we want and read the csv for each file
    cols = ['nTSS', 'nTags', 'tsrwidth']
    ts1_data = pd.read_csv(ts1_csv, usecols = cols)
    ts2_data = pd.read_csv(ts2_csv, usecols = cols)
    ts3_data = pd.read_csv(ts3_csv, usecols = cols)
    ts4_data = pd.read_csv(ts4_csv, usecols = cols)

    # get data values to manipulate later
    data_val1 = ts1_data.values
    data_val2 = ts2_data.values
    data_val3 = ts3_data.values
    data_val4 = ts4_data.values

    # now use MinMaxScaler to scale each data set for range (0, 1)
    from scipy import stats
    from sklearn.preprocessing import MinMaxScaler
    scaler1 = MinMaxScaler(feature_range=(0,1))
    scaler1 = scaler1.fit(data_val1)
    scaler2 = MinMaxScaler(feature_range=(0,1))
    scaler2 = scaler2.fit(data_val2)
    scaler3 = MinMaxScaler(feature_range=(0,1))
    scaler3 = scaler3.fit(data_val3)
    scaler4 = MinMaxScaler(feature_range=(0,1))
    scaler4 = scaler4.fit(data_val4)

    # finally normalize the data
    norm1 = scaler1.transform(data_val1)
    norm2 = scaler2.transform(data_val2)
    norm3 = scaler3.transform(data_val3)
    norm4 = scaler4.transform(data_val4)

    print(norm1)
    print(norm2)
    with open(ts1_csv, 'r') as tsr1, open(ts2_csv, 'r') as tsr2, open(ts3_csv, 'r') as tsr3, open(ts4_csv, 'r') as tsr4:
        tsr1_lines = csv.reader(tsr1)
        tsr2_lines = csv.reader(tsr2)
        tsr3_lines = csv.reader(tsr3)
        tsr4_lines = csv.reader(tsr4)
        next(tsr1_lines)
        next(tsr2_lines)
        next(tsr3_lines)
        next(tsr4_lines)
        tsr1_data = dict()
        tsr2_data = dict()
        tsr3_data = dict()
        tsr4_data = dict()
        i = 0

        for line in tsr1_lines:
            chrom = line[0]
            ID = line[1]
            start = line[2]
            end = line[3]
            strand = line[4]
            ntss = line[5]
            ntags = line[6]
            tsrwidth = line[7]
            if ID not in tsr1_data:
                tsr1_data[ID] = { '+': [], '-':[]}
            normal = norm1[i][0] + norm1[i][1] + norm1[i][2]
            tsr1_data[ID][strand].append([chrom, ID, start, end, strand, ntss, ntags, tsrwidth, normal])
            i+=1

        i=0
        for line in tsr2_lines:
            chrom = line[0]
            ID = line[1]
            start = line[2]
            end = line[3]
            strand = line[4]
            ntss = line[5]
            ntags = line[6]
            tsrwidth = line[7]
            if ID not in tsr2_data:
                tsr2_data[ID] = { '+': [], '-':[]}
            normal = norm2[i][0] + norm2[i][1] + norm2[i][2]
            tsr2_data[ID][strand].append([chrom, ID, start, end, strand, ntss, ntags, tsrwidth, normal])
            i+=1

        i=0
        for line in tsr3_lines:
            chrom = line[0]
            ID = line[1]
            start = line[2]
            end = line[3]
            strand = line[4]
            ntss = line[5]
            ntags = line[6]
            tsrwidth = line[7]
            if ID not in tsr3_data:
                tsr3_data[ID] = { '+': [], '-':[]}
            normal = norm3[i][0] + norm3[i][1] + norm3[i][2]
            tsr3_data[ID][strand].append([chrom, ID, start, end, strand, ntss, ntags, tsrwidth, normal])
            i+=1

        i=0
        for line in tsr4_lines:
            chrom = line[0]
            ID = line[1]
            start = line[2]
            end = line[3]
            strand = line[4]
            ntss = line[5]
            ntags = line[6]
            tsrwidth = line[7]
            if ID not in tsr4_data:
                tsr4_data[ID] = { '+': [], '-':[]}
            normal = norm4[i][0] + norm4[i][1] + norm4[i][2]
            tsr4_data[ID][strand].append([chrom, ID, start, end, strand, ntss, ntags, tsrwidth, normal])
            i+=1

    regions = []
    with open('promoterRegions.txt', 'r') as txt:
        lines = txt.readlines()
        for line in lines:
            t = line.split()
            regions.append(t[1])
    test_ts1_tsr2_ps, test_ts1_tsr2_ns = coxonTest(tsr1_data, tsr2_data, regions)
    test_ts1_tsr3_ps, test_ts1_tsr3_ns = coxonTest(tsr1_data, tsr3_data, regions)
    test_ts1_tsr4_ps, test_ts1_tsr4_ns = coxonTest(tsr1_data, tsr4_data, regions)
    test_ts2_tsr3_ps, test_ts2_tsr3_ns = coxonTest(tsr2_data, tsr3_data, regions)
    test_ts2_tsr4_ps, test_ts2_tsr4_ns = coxonTest(tsr2_data, tsr4_data, regions)
    test_ts3_tsr4_ps, test_ts3_tsr4_ns = coxonTest(tsr3_data, tsr4_data, regions)

    count1 = 0
    count2 = 0
    for x in test_ts3_tsr4_ps:
        if test_ts3_tsr4_ps[x] < .5:
            count1+=1
    for x in test_ts2_tsr4_ps:
        if test_ts2_tsr4_ps[x] < .5:
            count2+=1
    print(count1)
    print(count2)

def coxonTest(tsr1_data, tsr2_data, regions):
    positiveT = {}
    negativeT = {}


    for r in regions:
        flag1 = False
        flag2 = False

        if tsr1_data.get(r, "") != "":
            if len(tsr1_data[r]['+']) == 0:
                narrp1 = np.array([[0.0]], dtype=float)
            else:
                narrp1 = np.array([[tsr1_data[r]['+'][0][8]]], dtype=float)

            if len(tsr1_data[r]['-']) == 0:
                narrn1 = np.array([[0.0]], dtype=float)
            else:
                narrn1 = np.array([[tsr1_data[r]['-'][0][8] ]], dtype=float)

            flag1 = True


        if tsr2_data.get(r, "") != "":
            if len(tsr2_data[r]['+']) == 0:
                narrp2 = np.array([[0.0]], dtype=float)
            else:
                narrp2 = np.array([[tsr2_data[r]['+'][0][8]]], dtype=float)

            if len(tsr2_data[r]['-']) == 0:
                narrn2 = np.array([[0.0]], dtype=float)
            else:
                narrn2 = np.array([[tsr2_data[r]['-'][0][8]]], dtype=float)

            flag2 = True


        # if flag1 == False:
        #     narrp1 = np.array([[0.00001]], dtype=float)
        #     narrn1 = np.array([[0.00001]], dtype=float)
        #
        # if flag2 == False:
        #     narrp2 = np.array([[0.0]], dtype=float)
        #     narrn2 = np.array([[0.0]], dtype=float)

        if flag1 == False and flag2 == False or narrp1[0][0] == 0 and narrp2[0][0] == 0 or narrn1[0][0] == 0 and narrn2[0][0] == 0:
            positiveT[r] = 100
            negativeT[r] = 100
        else:
            w,p = wilcoxon(narrp1[:,0], y=narrp2[:,0])
            positiveT[r] = p

            w,p = wilcoxon(narrn1[:,0], y=narrn2[:,0])
            negativeT[r] = p

    return positiveT, negativeT

    '''
    Later we will be able to use these normalized values to run paired t test
        stats.ttest_rel(norm1[:,0], norm1[:,1])
        ^^ this would run paired t test on column 1 and column 2

    Professor said can use the un-normalized data for the wilcoxon test in the email
        from scipy.stats import wilcoxon
        w, p = wilcoxon(data_val1[:,0], y=data_val2[:,1])
        ^^
    '''

if __name__ == "__main__":
    ts1_csv = "tsr1_peaks.csv"
    ts2_csv = "tsr2_peaks.csv"
    ts3_csv = "tsr3_peaks.csv"
    ts4_csv = "tsr4_peaks.csv"
    normalize(ts1_csv, ts2_csv, ts3_csv, ts4_csv)
