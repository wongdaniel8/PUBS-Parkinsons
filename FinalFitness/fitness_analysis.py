#!/usr/bin/env python3
# deltaFitness = Drug - NoDrug

import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()
from collections import Counter
from collections import defaultdict
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os
from sys import argv

def createHeatmap(df, name):
    df = pd.DataFrame(df)
    plt.figure(figsize=(20, 8))
    ax = sns.heatmap(df, center = 0)

    plt.setp(ax.get_xticklabels(), rotation=90)
    plt.setp(ax.get_yticklabels(), rotation=0)
    plt.title("Frequency of Barcodes")
    plt.show()
    plt.savefig("heatmap_%s.png" %name)

#Unpickle
R1 = pkl.load(open("anm_fit_df1.pkl", "rb"))
R2 = pkl.load(open("anm_fit_df2.pkl", "rb"))

C1 = pkl.load(open("ctrl_fit_df1.pkl", "rb"))
C2 = pkl.load(open("ctrl_fit_df2.pkl", "rb"))
C1.replace(-10, np.nan, inplace=True)
print(C1)

# Create heatmaps
#createHeatmap(R1, "R1")
#createHeatmap(R2, "R2")
#createHeatmap(C1, "C1")
#createHeatmap(C2, "C2")

# Combine control datasets
#dfC = pd.concat((C1, C2))
#dfC.groupby(dfC.index).mean()

# Subtract each replicate Drug - NoDrug
#rd1 = R1.subtract(dfC)
