import numpy as np
import scanpy as sc
import scipy as sci
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import autogenes as ag

from sklearn.svm import NuSVR
import pickle

import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts informative genes with autogeneS')
    parser.add_argument('--input_dir',default = 'autogeneS_input', required=True)
    parser.add_argument('--ngen',default = 5000, type =int, help='Number of generations',required=False)
    parser.add_argument('--seed',default = 0, type =int, help='a random seed value',required=False)
    parser.add_argument('--nfeatures',default = 400, type =int, help='Number of informative genes',required=False)
    parser.add_argument('--mode',default = 'fixed', help='autogeneS running mode',required=True)


    args = parser.parse_args()


    autogeneS_input_files=os.listdir(args.input_dir)

    for input_file in autogeneS_input_files:
        centroids_sc_hv=pd.read_csv(args.input_dir+'/'+input_file,sep=',')
        
        ag.init(centroids_sc_hv.T)

        if args.mode == 'fixed':
            ag.optimize(ngen=args.ngen, seed=args.seed, nfeatures=args.nfeatures, mode= 'fixed', offspring_size=100, verbose=False) 
        elif args.mode == 'standard':
            ag.optimize(ngen=args.ngen, seed=args.seed, nfeatures=None, mode= 'standard', offspring_size=100, verbose=False) 

        index = ag.select(index=0) # Pareto index
        centroids_sc_pareto = centroids_sc_hv[index]

        if os.path.exists('autogeneS_output'):
            centroids_sc_pareto.to_csv('autogeneS_output/' + input_file)
        else:
            os.mkdir('autogeneS_output')
            centroids_sc_pareto.to_csv('autogeneS_output/' + input_file)


