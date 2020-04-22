'''
Requires scikit-learn, scipy, numpy, and pandas
'''

import pandas as pd
import numpy as np 


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
    scaler = MinMaxScaler(feature_range=(0,1))
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