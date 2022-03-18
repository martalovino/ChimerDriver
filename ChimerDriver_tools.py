# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:04:32 2020

@author: vener
"""

import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel
import csv
from sklearn import svm
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras import optimizers
from keras.models import load_model
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt
import random
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_predict



class DataIntegration:
    def __init__(self, filename1, filename2, filename3, filename4):
        self.filename1 = filename1
        self.filename2 = filename2
        self.filename3 = filename3
        self.filename4 = filename4

    def read_data(self):
        X = pd.read_csv(self.filename1, sep='\t')
        for col in ['strand5p','strand3p','name5p','name3p']:
            if col in X.columns:
                X = X.drop(columns=[col])        
        # load GO matrix
        GO = pd.read_csv(self.filename2, sep='\t')
        GO = GO.drop(columns=['FusionPair','Label'])
        # load TF matrix
        TF_matrix = pd.read_csv(self.filename3, sep='\t')
        TF_matrix = TF_matrix.drop(['FusionPair', 'Label'], axis=1)
        # load miRNA matrix
        miRNA_matrix = pd.read_csv(self.filename4, sep='\t')
        miRNA_matrix = miRNA_matrix.drop(['FusionPair','Label'], axis=1)
        # concatanate all of the features in one dataframe
        X_TF_GO_miRNA = pd.concat([X, TF_matrix, GO, miRNA_matrix], axis=1, sort=False).set_index("FusionPair")
        Y = X_TF_GO_miRNA['Label']
        X_TF_GO_miRNA = X_TF_GO_miRNA.drop(columns=['Label'])
        # drop inconsistent data:
        try:
            ind_0 = X_TF_GO_miRNA[(X_TF_GO_miRNA['perc5p']<0) |(X_TF_GO_miRNA['perc3p']<0) ].index.to_list()
            ind_100 =X_TF_GO_miRNA[(X_TF_GO_miRNA['perc5p']>100) |(X_TF_GO_miRNA['perc3p']>100)].index.to_list()
            ind_nan=X_TF_GO_miRNA[(np.isnan(X_TF_GO_miRNA['perc5p'])) |(np.isnan(X_TF_GO_miRNA['perc3p']))].index.to_list()    
            X_TF_GO_miRNA = X_TF_GO_miRNA.drop(ind_0+ind_100+ind_nan)
            Y = Y.drop(ind_0+ind_100+ind_nan)
        except:
            pass
        X_TF_GO_miRNA = X_TF_GO_miRNA.fillna(0)
        return X_TF_GO_miRNA, Y

class JoinDatasets():
    def __init__(self, filename):
        self.filename = filename
        self.filenames=[]
    def isToMerge(self):
        if '+' in self.filename:
            return True
    def MergeDatasets(self, folder_name):        
        for times in range(self.filename.count('+')+1):

            if '/' in self.filename:
                filename_clean = self.filename.split('+')[times].split('/')[-1] # if filename includes a folder
            else:
                filename_clean = self.filename.split('+')[times]

            self.filenames.append(('.').join(filename_clean.split('.')[0:-1])) #separate from extension, if multiple dots appear a part from the extension then reunite them

        filename_out_base = ('_').join(self.filenames)        

        try: # try if file has already been created
            df_dummy = pd.read_csv(folder_name+'/'+filename_out_base+'_structfeat.csv')
        except:
            featsets = ['structfeat','TF','GO','miRNA']
            for featset in featsets:
                mergedDatasets = pd.DataFrame()
                for times in range(self.filename.count('+')+1):
                    filename1 = self.filename.split('+')[times].split('/')[-1].split('.')[0]
                    df = pd.read_csv(folder_name+'/'+filename1+'_'+featset+'.csv',sep='\t')
                    mergedDatasets = pd.concat([mergedDatasets,df]).reset_index(drop=True)
                mergedDatasets.to_csv(folder_name+'/'+filename_out_base+'_'+featset+'.csv', sep='\t',index=False)
        return filename_out_base + '_'

class FeatureSelection(DataIntegration):
    def __init__(self, filename_out):
        self.filename_out = filename_out

    def Select_with_trees(self, filename1, filename2, filename3, filename4, thresh=0.001, n_est=50):
        data = DataIntegration(filename1, filename2, filename3, filename4)
        X, Y = data.read_data()
        clf = ExtraTreesClassifier(n_estimators=n_est)
        clf = clf.fit(X, Y)
        SelectFromModel(clf, prefit=True)
        # X_new = model.transform(X)
        # find the features with tree value > thresh
        features = X.columns.to_list()
        features_selected = []
        for ind in range(len(clf.feature_importances_)):
            if clf.feature_importances_[ind] > thresh:
                features_selected.append(features[ind])
        cont = 0
        with open(self.filename_out, "w") as output:
            for name in features_selected:
                writer = csv.writer(output, lineterminator='\n')
                writer.writerow([name])
                cont += 1
        print(">> %d features added." % cont)
        return features_selected
    
    def Select_by_name(self, features_wanted, filename1, filename2, filename3, filename4):
        # e.g. features_wanted = ['mirNA','TF','initial_features','GO']
        data = DataIntegration(filename1, filename2, filename3, filename4)
        X, Y = data.read_data()        
        features_selected = []
        GOs = [feat for feat in X.columns.to_list() if 'GO' in feat]
        miRNA = [feat for feat in X.columns.to_list() if 'miR' in feat or 'let' in feat]
        initial_features = X.columns.to_list()[0:5]
        TF = [feat for feat in X.columns.to_list() if feat not in GOs+miRNA+initial_features]
        if 'miRNA' in features_wanted:
            for feat in miRNA:
                features_selected.append(feat)
        if 'TF' in features_wanted:
            for feat in TF:
                features_selected.append(feat)
        if 'initial_features' in features_wanted:
            for feat in initial_features:
                features_selected.append(feat)
        if 'GOs' in features_wanted:
            for feat in GOs:
                features_selected.append(feat)
            
        cont = 0
        with open(self.filename_out, "w") as output:
            for name in features_selected:
                writer = csv.writer(output, lineterminator='\n')
                writer.writerow([name])
                cont += 1
        print(">> %d features added." % cont)
        
        return features_selected
     
    def Select_by_name_and_random_forest(self, features_wanted, filename1, filename2, filename3, filename4, thresh=0.001, n_est=50):
        # e.g. features_wanted = ['mirNA','TF','initial_features','GOs']
        data = DataIntegration(filename1, filename2, filename3, filename4)
        X, Y = data.read_data()        
        features_selected = []
        GOs = [feat for feat in X.columns.to_list() if 'GO' in feat]
        miRNA = [feat for feat in X.columns.to_list() if 'miR' in feat or 'let' in feat]
        initial_features = X.columns.to_list()[0:5]
        TF = [feat for feat in X.columns.to_list() if feat not in GOs+miRNA+initial_features]
        if 'miRNA' in features_wanted:
            for feat in miRNA:
                features_selected.append(feat)
        if 'TF' in features_wanted:
            for feat in TF:
                features_selected.append(feat)
        if 'initial_features' in features_wanted:
            for feat in initial_features:
                features_selected.append(feat)
        if 'GOs' in features_wanted:
            for feat in GOs:
                features_selected.append(feat)     
        # obtain an X set with only the subset of features wanted and then reduce with random forest
        X = X[features_selected]
        clf = ExtraTreesClassifier(n_estimators=n_est)
        clf = clf.fit(X, Y)
        SelectFromModel(clf, prefit=True)
        # X_new = model.transform(X)
        # find the features with tree value > thresh
        features = X.columns.to_list()
        features_selected = []
        for ind in range(len(clf.feature_importances_)):
            if clf.feature_importances_[ind] > thresh:
                features_selected.append(features[ind])
        cont = 0
        with open(self.filename_out, "w") as output:
            for name in features_selected:
                writer = csv.writer(output, lineterminator='\n')
                writer.writerow([name])
                cont += 1
        print(">> %d features added." % cont)
        return features_selected
       
        
     
class initial_features:
    def __init__(self, filename):
        self.filename = filename
        self.matrix = pd.read_csv("gene_attribute_matrix.csv", index_col=0)
        # use cancermine to assign roles
        self.cancer = pd.read_csv("cancermine.csv", index_col = 0)
        self.dictionary = {'Other':0,'Tumor_Suppressor':1,'Driver':2,'Oncogene':3}  

    @staticmethod
    def grchChr(version,nchr, grch38_r,grch37_r):
        if version == 'grch38':
            if nchr == 'X' or nchr == 'Y':
                grch_i = grch38_r[grch38_r.chr == nchr]
            else:
                grch_i = grch38_r[grch38_r.chr == int(nchr)]
        elif version == 'grch37':
            if nchr == 'X' or nchr == 'Y':
                grch_i = grch37_r[grch37_r.chr == nchr]
            else:
                grch_i = grch37_r[grch37_r.chr == int(nchr)]
        else:
            print('error in reading genome reference')
            grch_i = ''
        return grch_i
    
    def Role(self,gene_name):
        ind = None
        cont = 0    
        for index, gene in self.cancer.iterrows():
            if gene.gene_normalized == gene_name:
                ind = cont    
                break
            cont += 1
        if ind == None:
            role = 'Other'
        else:
            role = self.cancer.iloc[ind]['role']
        return role  
    
    def RetainedPerc(self, version, gene_name, nchr):
        grch38 = pd.read_csv("grch38.csv", usecols = ['chr','type','start','end','strand','gene_name'], chunksize=10000)
        grch37 = pd.read_csv("grch37.csv", usecols = ['chr','type','start','end','strand','gene_name'], chunksize=10000)
        grch38_r = pd.DataFrame()
        grch37_r = pd.DataFrame()
        for chunck in grch38:
            chunck = chunck[(chunck['type'] == 'gene')]
            grch38_r = grch38_r.append(chunck)  
        for chunck in grch37:
            chunck = chunck[(chunck['type'] == 'gene')]
            grch37_r = grch37_r.append(chunck) 
        grch_i = self.grchChr(version, nchr, grch38_r, grch37_r)           
        cont, ind = 0, None
        for index, gene in grch_i.iterrows():
            if gene.gene_name == gene_name:
                ind = cont
                break
            cont += 1
        if ind == None:
            print(f"gene {gene_name} was not found in geneome reference database")
            gene_start, gene_end, gene_strand = None, None, None
        else:
            gene_start = grch_i.iloc[ind]['start']
            gene_end = grch_i.iloc[ind]['end']
            gene_strand = grch_i.iloc[ind]['strand']

        return gene_start, gene_end, gene_strand
    
    def obtaining_partial_features(self, label, new_filename):
        # There are 3 possible datasets: DeePrior, Pegasus, Oncofuse
        error_message = ('The input file format did not fall into any of the usual cases.\n'
            'Please enter a file in a .csv format (or fusion_candidates.xlsx) with either one of the following sets of minimum required information separated by a tab:\n\n'
            '\t1. gene names, labels, breakpoints, chromosomes, genome version in 7 columns defined as '
            '"5pCommonName", "3pCommonName", "Coord5p", "Coord3p", "Chr5p", "Chr3p", "Version"\n\n'
            '\t2. gene names, labels, breakpoints, chromosomes, gene start, gene end in 10 columns defined as '
            '"Gene_Name1", "Gene_Name2", "Gene_Breakpoint1", "Gene_Breakpoint2", "Chr1", "Chr2", "Gene_Start1", "Gene_End1", "Gene_Start2", "Gene_End2"\n\n'
                 )
        min_pegasus_info = ["Gene_Name1", "Gene_Name2", "Gene_Breakpoint1", "Gene_Breakpoint2", "Chr1", "Chr2", "Gene_Start1", "Gene_End1", "Gene_Start2", "Gene_End2", "Strand1","Strand2"]
        min_deeprior_info = ["5pCommonName", "3pCommonName", "Coord5p", "Coord3p", "Chr5p", "Chr3p", "Version"] #"Label"
        if 'oncofuse' in self.filename:
            # Do Oncofuse stuff to parse the file
            myset = pd.read_csv(self.filename, sep='\t')
            gene_name5p_column_name, gene_name3p_column_name = "FPG5_NAME", "FPG3_NAME"
            label_column_name = "Label"   
        else:
            try:
                myset = pd.read_csv(self.filename, sep='\t')
                if set(min_pegasus_info).issubset(set(myset.columns.to_list())):
                    gene_name5p_column_name, gene_name3p_column_name = "Gene_Name1", "Gene_Name2"
                    breakpoint5p_column_name, breakpoint3p_column_name = "Gene_Breakpoint1", "Gene_Breakpoint2"
                    chr5p_column_name, chr3p_column_name = "Chr1", "Chr2"
                    gene_start5p_column_name, gene_start3p_column_name, gene_end5p_column_name, gene_end3p_column_name = "Gene_Start1", "Gene_Start2", "Gene_End1", "Gene_End2"
                    strand5p_column_name, strand3p_column_name = "Strand1", "Strand2"
                    label_column_name = ''
                    version_column_name = ''
                elif set(min_deeprior_info).issubset(set(myset.columns.to_list())):
                    gene_name5p_column_name, gene_name3p_column_name = "5pCommonName", "3pCommonName"
                    breakpoint5p_column_name, breakpoint3p_column_name = "Coord5p", "Coord3p"
                    chr5p_column_name, chr3p_column_name = "Chr5p", "Chr3p"
                    version_column_name = "Version"
                    if "Label" in myset.columns.to_list():
                        label_column_name = "Label"
                    else:
                        label_column_name = "NotFound"
                else:
                    print('unable to read set'+error_message)
            except:
                print(error_message)
        
        Fusion_pairs, Labels, strands5p, strands3p, same_strand, names5p, names3p, role5p, role3p, perc5p, perc3p = [], [], [], [], [], [], [], [], [], [], []
    
        for i in range(len(myset)):
            if i%20==0:
                print("iter %d/%d" %(i,len(myset)))
            # the fusion pair in the form Gene1_Gene2
            Fusion_pairs.append(myset.iloc[i][gene_name5p_column_name]+'_'+myset.iloc[i][gene_name3p_column_name])
            # the class/label, either 1 or 0
            if len(label_column_name)==0 or label==1 or label==0: # in Pegasus the information is given by the filename
                Labels.append(label)
            elif label_column_name=="NotFound":
                Labels.append(1) # if the dataset was in the DEEPrior format but "Label" column was not found assign 1
            else: # in DEEPrior the information is given as a column
                Labels.append(myset.iloc[i][label_column_name])

            # common names and ensg names of the gene pair
            names5p.append(myset.iloc[i][gene_name5p_column_name])
            names3p.append(myset.iloc[i][gene_name3p_column_name])

            # the role assigned by cancermine table
            role5p.append(self.dictionary[self.Role(myset.iloc[i][gene_name5p_column_name])])
            role3p.append(self.dictionary[self.Role(myset.iloc[i][gene_name3p_column_name])])
            
            try:
                gene_start5p, gene_end5p, strand_5p = self.RetainedPerc(myset.iloc[i][version_column_name], myset.iloc[i][gene_name5p_column_name], myset.iloc[i][chr5p_column_name])
                gene_start3p, gene_end3p, strand_3p = self.RetainedPerc(myset.iloc[i][version_column_name], myset.iloc[i][gene_name3p_column_name], myset.iloc[i][chr3p_column_name])
            except:
                try:
                    gene_start5p = myset.iloc[i][gene_start5p_column_name]
                    gene_start3p = myset.iloc[i][gene_start3p_column_name]
                    gene_end5p = myset.iloc[i][gene_end5p_column_name]
                    gene_end3p = myset.iloc[i][gene_end3p_column_name]
                    strand_5p = myset.iloc[i][strand5p_column_name]
                    strand_3p =  myset.iloc[i][strand3p_column_name]
                except NameError:
                    gene_start5p, gene_end5p, strand_5p,gene_start3p, gene_end3p, strand_3p='','','','','',''
                    
            try:
                if strand_5p == '+':
                    retained_gene5p = (myset.iloc[i][breakpoint5p_column_name] - gene_start5p)/(gene_end5p-gene_start5p)*100
                elif strand_5p == '-':
                    retained_gene5p = (gene_end5p - myset.iloc[i][breakpoint5p_column_name])/(gene_end5p-gene_start5p)*100
                if strand_3p == '+':
                    retained_gene3p = (gene_end3p - myset.iloc[i][breakpoint3p_column_name])/(gene_end3p-gene_start3p)*100
                elif strand_3p == '-':
                    retained_gene3p = (myset.iloc[i][breakpoint3p_column_name] - gene_start3p)/(gene_end3p-gene_start3p)*100
                perc3p.append(retained_gene3p)
                perc5p.append(retained_gene5p)
                strands5p.append(strand_5p)
                strands3p.append(strand_3p)
                if strand_5p == strand_3p:
                    same_strand.append(1)
                else:
                    same_strand.append(0)
                my_set = pd.DataFrame(data={'FusionPair':Fusion_pairs,
                                                  'Label':Labels,
                                                  'strand5p':strands5p,
                                                  'strand3p':strands3p,                                                                      
                                                  'same_strand':same_strand,
                                                  'name5p':names5p,
                                                  'name3p':names3p,    
                                                  'role5p':role5p,
                                                  'role3p':role3p,                                        
                                                  'perc5p':perc5p,
                                                  'perc3p':perc3p                              
                                                  })
            except NameError: # oncofuse case we have none of these info
                my_set = pd.DataFrame(data={'FusionPair':Fusion_pairs,
                                                  'Label':Labels,
                                                  'name5p':names5p,
                                                  'name3p':names3p,    
                                                  'role5p':role5p,
                                                  'role3p':role3p,                                                                     
                                                  })
             
            my_set.to_csv(new_filename, sep='\t', index=False)
        

class Transcription_factors_features:
    def __init__(self, filename):
        self.set_to_process = pd.read_csv(filename, sep='\t')
        while 'Unnamed' in self.set_to_process.columns[0]:
            self.set_to_process = self.set_to_process.drop([self.set_to_process.columns[0]],axis=1)
        self.matrix = pd.read_csv("gene_attribute_matrix.csv", index_col=0)

    @staticmethod
    def get_genes(myset):
        pairs = myset.FusionPair.tolist()
        genesA = [pair.split('_')[0] for pair in pairs]
        genesB = [pair.split('_')[1] for pair in pairs]
        genes = list(set(genesA + genesB))
        return genes

    @staticmethod
    def paired_matr(red_matrix, myset):
        # a matr with gene pairs on the rows e TF on the columns + a column with the labels
        cols = red_matrix.columns
        factor_matr = pd.DataFrame(columns=cols, index=myset['FusionPair'])
        Label = pd.DataFrame(index=myset['FusionPair'], columns=['Label'])
        for i in range(len(myset['FusionPair'])):
            if i % 20 == 0:
                print('iter %d/%d' % (i, len(myset['FusionPair'])))
            gene1 = myset.iloc[i]['FusionPair'].split('_')[0] # uses only the gene at the 5' position
            gene1_transcr = red_matrix.loc[gene1].to_list()
            factor_matr.iloc[i] = gene1_transcr
            Label.iloc[i]['Label'] = myset.iloc[i]['Label']
        matr = pd.concat([factor_matr, Label], axis=1)
        modified = matr.reset_index()
        modified = modified.set_index('FusionPair')
        return modified

    def build_TF_features(self, filename_out):
        gene_myset = self.get_genes(self.set_to_process)
        # obtaining a reduced matrix to speed up computational time that includes only the genes of the datasets
        new_index = gene_myset
        new_index = list(set(new_index))  # there were some duplicated genes
        red_matrix = self.matrix.reindex(new_index)
        red_matrix = red_matrix.fillna(0)  # where the TF association was unknown it is assumed 0
        factor_matr_final = self.paired_matr(red_matrix, self.set_to_process)
        factor_matr_final.to_csv(filename_out,sep='\t')

class Gene_ontologies_features:
    # generates the gene ontologies matrices for both X_train and x_test using the GO found in the training set
    def __init__(self, filename, whichset):
        self.set_to_process = pd.read_csv(filename,sep='\t')
        while 'Unnamed' in self.set_to_process.columns[0]:
            self.set_to_process = self.set_to_process.drop([self.set_to_process.columns[0]],axis=1)
        self.biomart = 'gene_attribute_matrix_info_biomart.txt'
        self.whichset = whichset

    def GO_from_mart_to_list(self, set_genes):  #trainset, testset1, testset2
        GO_set = []
        with open(self.biomart,'r') as f:
            cols = f.readline().split('\n')[0].split(',')
            for line in f:
                if line.split('\n')[0].split(',')[1] in set_genes:
                    if len(line.split('\n')[0].split(',')[-1])>1: # se non c'è nessuna GO va avanti
                        GO_set.append(line.split('\n')[0].split(','))
        GO_s = pd.DataFrame(data=GO_set,columns=cols)
        GO_s = GO_s.drop_duplicates(subset=['GO term accession','Gene name'])
        return GO_s

    @staticmethod
    def GO_matr(myset, GO_set, GO_ids):  # X_train/x_test
        fusion_pairs = myset['FusionPair'].to_list()
        # GO dataframe made out of the GOs found in the training set
        X_GO = pd.DataFrame(index=myset['FusionPair'], columns=GO_ids)
        for i in range(len(fusion_pairs)):
            # per ogni gene della coppia
            gene_name1 = fusion_pairs[i].split('_')[0]
            gene_name2 = fusion_pairs[i].split('_')[1]
            # seleziono le GO di quel gene
            GOs = GO_set[GO_set['Gene name'] == gene_name1]['GO term accession'].to_list() + \
                  GO_set[GO_set['Gene name'] == gene_name2]['GO term accession'].to_list()
            for go in GOs:
                X_GO.iloc[i][go] = 1
        X_GO = X_GO.fillna(0)
        X_GO['Label'] = myset['Label'].to_list()
        return X_GO

    def build_GO_features(self, GO_filename_out):
        folder_name = GO_filename_out.split('/')[0]
        # GO_tr/ts containing in each row gene names - gene term - ensg (one row for each GO, gene names repeted)
        GO_set = self.GO_from_mart_to_list(self.set_to_process['name3p'].to_list()+self.set_to_process['name5p'].to_list())
        if self.whichset == 'train':
            GO_ids = list(set(GO_set['GO term accession'].to_list()))  # using only the GOs from the training set
            with open(folder_name+'/GO_ids.txt', 'w') as f:
                for item in GO_ids:
                    f.write("%s\n" % item)
        elif self.whichset == 'train_again':  # In pegasus i file su cui trainare sono 2, allora bisogna tenere conto di tutti gli ID ottenuti
            GO_ids = list(set(GO_set['GO term accession'].to_list()))  # using only the GOs from the training set
            with open(folder_name+'/GO_ids.txt', 'a') as f:
                for item in GO_ids:
                    f.write("%s\n" % item)                    
        else:
            GO_ids = []
            with open(folder_name+'/GO_ids.txt', 'r') as f:
                for row in f:
                    line = row.strip('\n')
                    GO_ids.append(line)

        # X_train/test have Fusion pair on the rows and GOs on the columns with an additional Label column
        # the cell corresponding to FusionPairN - GO:M is marked 1 if any of the gene in the pair had GO:M among its ontologies
        X_set_GO = self.GO_matr(self.set_to_process, GO_set, GO_ids)
        X_set_GO.to_csv(GO_filename_out,sep='\t')

        return (len(X_set_GO.columns) - 1)


class miRNA_features(Transcription_factors_features):
    def __init__(self, filename):
        self.set_to_process = pd.read_csv(filename, sep='\t')  # obtain X_train
        while 'Unnamed' in self.set_to_process.columns[0]:
            self.set_to_process = self.set_to_process.drop([self.set_to_process.columns[0]],axis=1)
        self.matrix = pd.read_csv("miRNA_gene_matrix.csv", sep='\t', index_col=0)
        self.matrix = self.matrix.fillna(
            0)  # si assume che la prob sia nulla per quelle coppie per le quali non è definita

    @staticmethod
    def paired_matr(red_matrix,
                    myset):  # coppie di geni del set sulle righe e miRNA sulle colonne + una colonna con i label
        cols = red_matrix.columns
        factor_matr = pd.DataFrame(columns=cols, index=myset.FusionPair)
        Label = pd.DataFrame(index=myset.FusionPair, columns=['Label'])
        for i in range(len(myset.FusionPair)):
            if i % 20 == 0:
                print('iter %d/%d' % (i, len(myset.FusionPair)))
            gene1 = myset.iloc[i]['FusionPair'].split('_')[0]
            gene1_transcr = red_matrix.loc[gene1].to_list()
            gene2 = myset.iloc[i]['FusionPair'].split('_')[1]
            gene2_transcr = red_matrix.loc[gene2].to_list()
            genepair_transcr = []
            for ind in range(len(gene1_transcr)):  # assegnare la prob non nulla tra le due oppure quella massima
                if gene1_transcr[ind] == 0 and gene2_transcr[ind] > 0:
                    genepair_transcr.append(gene2_transcr[ind])
                elif gene1_transcr[ind] > 0 and gene2_transcr[ind] == 0:
                    genepair_transcr.append(gene1_transcr[ind])
                elif gene1_transcr[ind] > 0 and gene2_transcr[ind] > 0:
                    genepair_transcr.append(max(gene1_transcr[ind], gene2_transcr[ind]))
                else:
                    genepair_transcr.append(0.0)
            factor_matr.iloc[i] = genepair_transcr
            Label.iloc[i]['Label'] = myset.iloc[i]['Label']
        matr = pd.concat([factor_matr, Label], axis=1)
        modified = matr.reset_index()
        modified = modified.set_index('FusionPair')
        return modified

    def build_miRNA_features(self, filename_out):
        gene_set = self.get_genes(self.set_to_process)
        new_index = list(set(gene_set))  # c'erano geni duplicati
        red_matrix = self.matrix.reindex(new_index).fillna(0)
        factor_matr_set = self.paired_matr(red_matrix, self.set_to_process)
        factor_matr_set.to_csv(filename_out,sep='\t')
        return len(factor_matr_set.columns) - 1




class MLP(DataIntegration):
    def __init__(self, folder_name, train_filename_base='.', test_filename_base='.', val_filename_base='.', use_validation_set='.', feat_sel='all', learning_rate=0.001, training_epochs=300, batch_size=32
                 ):
        self.use_validation_set = use_validation_set
#        if use_validation_set=='True':
#             self.use_validation_set = True
#        else:
#            self.use_validation_set = False
#        self.use_validation_set = use_validation_set
        self.selected_features = feat_sel
        self.folder_name = folder_name
        self.trset_name = self.folder_name +'/'+ train_filename_base + 'structfeat.csv' #dataframes[0]

        self.trset_GO = self.folder_name +'/'+ train_filename_base + 'GO.csv' #dataframes[1]
        self.trset_TF = self.folder_name +'/'+ train_filename_base + 'TF.csv'#dataframes[2]
        self.trset_miRNA = self.folder_name +'/'+ train_filename_base + 'miRNA.csv'#dataframes[3]
        self.tsset_name = self.folder_name +'/'+ test_filename_base + 'structfeat.csv' #dataframes[4]
        self.tsset_GO = self.folder_name +'/'+ test_filename_base + 'GO.csv' #dataframes[5]
        self.tsset_TF = self.folder_name +'/'+ test_filename_base + 'TF.csv'#dataframes[6]
        self.tsset_miRNA = self.folder_name +'/'+ test_filename_base + 'miRNA.csv' #dataframes[7]
        self.valset_name = self.folder_name +'/'+ val_filename_base + 'structfeat.csv' #dataframes[8]
        self.valset_GO = self.folder_name +'/'+ val_filename_base + 'GO.csv' #dataframes[9]
        self.valset_TF = self.folder_name +'/'+ val_filename_base + 'TF.csv' #dataframes[10]
        self.valset_miRNA = self.folder_name +'/'+ val_filename_base + 'miRNA.csv' #dataframes[11]
                
        self.x_val, self.y_val = None, None
        self.learning_rate = learning_rate
        self.training_epochs = training_epochs
        self.batch_size = batch_size
        self.Y_pred, self.y_pred, self.Y_train,self.Y_test, self.results = '','','','',''
        self.model_filename=self.folder_name + '/trained_model_lr'+str(self.learning_rate)+'_nodes'


    def train_test_set(self):
        trdata = DataIntegration(self.trset_name, self.trset_GO, self.trset_TF, self.trset_miRNA)
        testdata = DataIntegration(self.tsset_name, self.tsset_GO, self.tsset_TF, self.tsset_miRNA)
        valdata = DataIntegration(self.valset_name, self.valset_GO, self.valset_TF, self.valset_miRNA)
        self.X_train, self.Y_train = trdata.read_data()
        self.X_test, self.Y_test = testdata.read_data()
        self.X_train, self.X_test = self.X_train[self.selected_features], self.X_test[self.selected_features]
        if self.use_validation_set=='True':
            random.seed(16)
            self.x_val, self.y_val = valdata.read_data()
            self.x_val = self.x_val[self.selected_features]
            ids_val = random.sample(range(len(self.X_train)), 200)
            ids_tr = [i for i in range(len(self.X_train)) if i not in ids_val]
            # validate on a combination of test set 1 and training set
            self.x_val = self.x_val.append(self.X_train.iloc[ids_val])
            self.y_val = self.y_val.append(self.Y_train.iloc[ids_val])
            self.X_train = self.X_train.iloc[ids_tr]
            self.Y_train = self.Y_train.iloc[ids_tr]
        if 'reactiveLN' in self.trset_name:  #balancing Pegasus training_set
            random.seed(16)
            print('reactiveLN file found in train set, now balancing Pegasus training_set...')
            pos_tr_samples = len(self.Y_train[self.Y_train==1])
            neg_tr_samples = len(self.Y_train[self.Y_train==0])
            neg_tr_samples_ind =list(self.Y_train[self.Y_train==0].reset_index().index)
            pos_tr_samples_ind =list(self.Y_train[self.Y_train==1].reset_index().index)
            print('Pegaus training set contained %d positive samples and %d negative samples' %(pos_tr_samples, neg_tr_samples))
            ids_bal = random.sample(neg_tr_samples_ind, pos_tr_samples)
            print('Randomly picked %d samples from the negative samples' % pos_tr_samples)
            ids = ids_bal + pos_tr_samples_ind
            self.X_train = self.X_train.reset_index().iloc[ids]
            self.Y_train = self.Y_train.reset_index().iloc[ids]  
            self.X_train = self.X_train.set_index(self.X_train.columns[0])
            self.Y_train = self.Y_train.set_index(self.Y_train.columns[0])
           
    def kfold_cross_validation(self, n_nodes, act_functs, k=10, drop = 0.2):
        self.train_test_set()
        l1,l2,l3,l4 = n_nodes
        a1,a2,a3,a4 = act_functs
        accuracies = []
        best_roc_auc = 0
        train_samples = self.X_train.reset_index().index.to_list()
        random.shuffle(train_samples)
        samples_per_fold = round(len(train_samples)/k)
        folds = [train_samples[x:x+samples_per_fold] for x in range(0, len(train_samples), samples_per_fold)]
        if len(folds)>10: # if an 11th fold was created with a few samples just add it to the second to last one
            last=folds[-1]
            second_to_last=folds[-2]
            remaining_folds = folds[:-2]
            new_last = second_to_last + last
            folds = remaining_folds + [new_last]
        cross_val_results = pd.DataFrame(
            columns=['lr', 'epochs', 'layer1', 'layer2', 'layer3', 'layer4', 'act_fun1', 'act_fun2',
            'act_fun3', 'act_fun4', 'confmat', 'accuracy', 'precision', 'recall','f1','AUC']) 
        for fold in folds:
            print('performing cross validation fold %d/%d' %(folds.index(fold)+1,k))
            kvalidation_set = self.X_train.iloc[fold]
            kvalidation_target = self.Y_train.iloc[fold]
            ktraining_set = self.X_train.iloc[[sample for sample in train_samples if sample not in fold]]
            ktraining_target = self.Y_train.iloc[[sample for sample in train_samples if sample not in fold]]
            n_input = len(ktraining_set.columns)
            model = Sequential()
            model.add(Dense(l1, input_dim=n_input, activation=a1))
            model.add(Dropout(drop))
            model.add(Dense(l2, activation=a2))
            model.add(Dropout(drop))
            model.add(Dense(l3, activation=a3))
            model.add(Dropout(drop))
            model.add(Dense(l4, activation=a4))
            model.add(Dropout(drop))
            model.add(Dense(1, activation='sigmoid'))
            opt = optimizers.Adam(lr=self.learning_rate, decay=1e-6)
            model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
            model.fit(ktraining_set, ktraining_target, epochs=self.training_epochs, batch_size=self.batch_size,
                      verbose=0)
            y_kpred = model.predict_classes(kvalidation_set)
            accuracies.append(metrics.accuracy_score(kvalidation_target, y_kpred))
            # save best cross validated model according to roc auc
            if metrics.roc_auc_score(kvalidation_target, y_kpred) > best_roc_auc:
                best_roc_auc = metrics.roc_auc_score(kvalidation_target, y_kpred)
                model.save(self.folder_name + '/cross_validated_model_best_rocauc.h5')
                
            print("acc for fold %d = %.3f" %(folds.index(fold)+1,metrics.accuracy_score(kvalidation_target, y_kpred)))
            cross_val_results.loc[folds.index(fold)] = [self.learning_rate, self.training_epochs, l1, l2, l3, l4, a1, a2, a3,
                          a4, metrics.confusion_matrix(kvalidation_target, y_kpred),
                          metrics.accuracy_score(kvalidation_target, y_kpred),
                          metrics.precision_score(kvalidation_target, y_kpred),
                          metrics.recall_score(kvalidation_target, y_kpred),
                          metrics.f1_score(kvalidation_target, y_kpred),
                          metrics.roc_auc_score(kvalidation_target, y_kpred)]

        return cross_val_results, np.mean(accuracies), np.std(accuracies)
       

    def train_model(self, n_nodes, drop=0.2):
        self.train_test_set()
        best_roc_auc = 0        
        layers = [n_nodes] #,[n_nodes[1], n_nodes[0], n_nodes[0], n_nodes[1]], [n_nodes[0]] * 4 , [n_nodes[0], n_nodes[0], n_nodes[1], n_nodes[2]]]
        #layers = [n_nodes]
        n_input = len(self.X_train.columns)
        self.results = pd.DataFrame(
            columns=['train/test', 'lr', 'dropout','epochs','epoch_stop', 'layer1', 'layer2', 'layer3', 'layer4', 'act_fun1', 'act_fun2',
                     'act_fun3', 'act_fun4', 'confmat', 'accuracy', 'precision', 'recall','f1','AUC'])
        cont = 0
        act_functs = [['sigmoid'] * 4, ['tanh'] * 4, ['relu']*4,['tanh', 'sigmoid', 'tanh', 'sigmoid'],['relu','sigmoid','relu','sigmoid'],['sigmoid','relu','sigmoid','relu'],['sigmoid','tanh','sigmoid','tanh'],['relu','tanh','relu','tanh'],['tanh','relu','tanh','relu']]
        #act_functs = [['sigmoid']*4]
        pat = 50
        for l1,l2,l3,l4 in layers:
        #for learn_rate in [0.002,0.004,0.006,0.008]:
            for a1, a2, a3, a4 in act_functs:
            #for dr in [0.3,0.4]:
                #l1,l2,l3,l4 = 512,512,512,512
                #a1,a2,a3,a4 ='relu','sigmoid','relu','sigmoid'
                model = Sequential()
                model.add(Dense(l1, input_dim=n_input, activation=a1))
                model.add(Dropout(drop))
                model.add(Dense(l2, activation=a2))
                model.add(Dropout(drop))
                model.add(Dense(l3, activation=a3))
                model.add(Dropout(drop))
                model.add(Dense(l4, activation=a4))
                model.add(Dropout(drop))
                model.add(Dense(1, activation='sigmoid'))
                opt = optimizers.Adam(lr=self.learning_rate, decay=1e-6)
                #opt = optimizers.Adam(lr=learn_rate, decay=1e-6)
                model.compile(loss='binary_crossentropy', optimizer=opt, metrics=['accuracy'])
                if self.use_validation_set=='True':
                    earlystop = EarlyStopping(monitor='val_accuracy', mode='auto', verbose=1, patience=pat, min_delta=0)
                    history = model.fit(self.X_train, self.Y_train, epochs=self.training_epochs, batch_size=self.batch_size,verbose=0, validation_data=(self.x_val, self.y_val),callbacks=[earlystop,ModelCheckpoint(filepath='best_model.h5',monitor='val_accuracy',save_best_only=True, verbose=0)])
                    model = load_model('best_model.h5')
                    epoch_stop = len(history.history['val_accuracy']) - pat
                else:
                    print('now fitting')
                    model.fit(self.X_train, self.Y_train, epochs=self.training_epochs, batch_size=self.batch_size,
                              verbose=0)
                    epoch_stop = self.training_epochs
                
                model.save(self.model_filename+'_dropout'+str(drop)+str(l1)+'-'+str(l2)+'-'+str(l3)+'-'+str(l4)+'_actfunc'+str(a1)+'-'+str(a2)+'-'+str(a3)+'-'+str(a4)+'.h5')
                self.Y_pred = model.predict_classes(self.X_train)
                self.results.loc[cont] = ['train', self.learning_rate, drop, self.training_epochs, epoch_stop, l1, l2, l3, l4, a1, a2, a3,
                                          a4, metrics.confusion_matrix(self.Y_train, self.Y_pred),
                                          metrics.accuracy_score(self.Y_train, self.Y_pred),
                                          metrics.precision_score(self.Y_train, self.Y_pred),
                                          metrics.recall_score(self.Y_train, self.Y_pred),
                                          metrics.f1_score(self.Y_train, self.Y_pred),                                           
                                          metrics.roc_auc_score(self.Y_train, self.Y_pred)]
                
                self.y_pred = model.predict_classes(self.X_test)            
                self.y_pred_proba = pd.DataFrame(model.predict_proba(self.X_test),columns=['probonco']) 
                try:
                    AUC_ = metrics.roc_auc_score(self.Y_test, self.y_pred)
                except:
                    AUC_ = ''
                self.results.loc[cont+1] = ['test', self.learning_rate, drop, self.training_epochs, epoch_stop, l1, l2, l3, l4, a1, a2,
                                          a3, a4, metrics.confusion_matrix(self.Y_test, self.y_pred),
                                          metrics.accuracy_score(self.Y_test, self.y_pred),
                                          metrics.precision_score(self.Y_test, self.y_pred),
                                          metrics.recall_score(self.Y_test, self.y_pred),
                                          metrics.f1_score(self.Y_test, self.y_pred), 
                                          AUC_]

               # save best cross validated model according to roc auc
                if metrics.roc_auc_score(self.Y_train, self.Y_pred) > best_roc_auc:
                    best_roc_auc = metrics.roc_auc_score(self.Y_train, self.Y_pred)
                    model.save(self.folder_name + '/trained_model_best_rocauc.h5') 
                cont += 2
                print('training stopped at epoch %d' %epoch_stop)
                print(self.results.accuracy[-2:])
                #print('Y_test:')
                #print(self.Y_test)
                #print('Y_pred:')
                #print(self.y_pred)
                pd.DataFrame(data=self.results).to_csv(self.folder_name + '/testing_MLP_Results_use_val_'+str(self.use_validation_set)+'_lr'+str(self.learning_rate)+'_drop'+str(drop)+str(l1)+str(l2)+str(l3)+str(l4)+str(a1)+str(a2)+str(a3)+str(a4)+'.csv', sep='\t')
                pd.DataFrame(self.X_test,columns=['X_test']).to_csv(self.folder_name + '/X_test_MLP__use_val_'+str(self.use_validation_set)+'_lr'+str(self.learning_rate)+'_drop'+str(drop)+str(l1)+str(l2)+str(l3)+str(l4)+str(a1)+str(a2)+str(a3)+str(a4)+'.csv', sep='\t')
                pd.DataFrame(self.Y_test.to_numpy(),columns=['Y_test']).to_csv(self.folder_name + '/Y_test_MLP__use_val_'+str(self.use_validation_set)+'_lr'+str(self.learning_rate)+'_drop'+str(drop)+str(l1)+str(l2)+str(l3)+str(l4)+str(a1)+str(a2)+str(a3)+str(a4)+'.csv', sep='\t')
                pd.DataFrame(self.y_pred,columns=['y_pred']).to_csv(self.folder_name + '/Y_pred_MLP__use_val_'+str(self.use_validation_set)+'_lr'+str(self.learning_rate)+'_drop'+str(drop)+str(l1)+str(l2)+str(l3)+str(l4)+str(a1)+str(a2)+str(a3)+str(a4)+'.csv', sep='\t')
                self.y_pred_proba.to_csv(self.folder_name + '/Y_prob_MLP__use_val_'+str(self.use_validation_set)+'_lr'+str(self.learning_rate)+'_drop'+str(drop)+str(l1)+str(l2)+str(l3)+str(l4)+str(a1)+str(a2)+str(a3)+str(a4)+'.csv', sep='\t')
        return self.results, self.X_train, self.X_test, self.x_val, self.Y_test, self.y_pred, self.y_pred_proba
        
    def test_model(self, model_to_load_filename): #n_nodes, drop, act_functs):
        #l1,l2,l3,l4 = n_nodes
        #a1,a2,a3,a4 = act_functs        
        #model = load_model(self.model_filename+'_dropout'+str(drop)+str(l1)+'-'+str(l2)+'-'+str(l3)+'-'+str(l4)+'_actfunc'+str(a1)+'-'+str(a2)+'-'+str(a3)+'-'+str(a4)+'.h5')
        
        self.train_test_set()
        model = load_model(model_to_load_filename)
        
        testdata = DataIntegration(self.tsset_name, self.tsset_GO, self.tsset_TF, self.tsset_miRNA)        
        self.X_test, self.Y_test = testdata.read_data()
        self.X_test = self.X_test[self.selected_features]
        
        self.y_pred = model.predict_classes(self.X_test)
        self.y_pred_proba = model.predict_proba(self.X_test)
        for i in range(len(X_test)):
            string_to_print = self.X_test['FusionPair'].iloc[i] + ' probability of oncogenity: '+ self.y_pred_proba['probonco'].iloc[i]
            print(string_to_print)
        try:
            AUC_score =  [metrics.roc_auc_score(self.Y_test, self.y_pred)]
        except:
            AUC_score = []
        results = {'learning rate': [self.learning_rate], 
                   'model name': [model_to_load_filename],
                   'confusion matrix': [metrics.confusion_matrix(self.Y_test, self.y_pred)],
                   'accuracy': [metrics.accuracy_score(self.Y_test, self.y_pred)],
                   'precision': [metrics.precision_score(self.Y_test, self.y_pred)],
                   'recall': [metrics.recall_score(self.Y_test, self.y_pred)],
                   'f1 score': [metrics.f1_score(self.Y_test, self.y_pred)], 
                   'AUC score': AUC_score}
        #print('testing results on model ' + model_to_load_filename +' :\n')
        #print(results)
        pd.DataFrame(data=self.results).to_csv(self.folder_name + '/testing_MLP_Results_use_val_'+str(self.use_validation_set)+'_lr'+str(self.learning_rate)+'_drop'+str(drop)+str(l1)+str(l2)+str(l3)+str(l4)+str(a1)+str(a2)+str(a3)+str(a4)+'.csv', sep='\t')
        pd.DataFrame(self.X_test,columns=['X_test']).to_csv(self.folder_name + '/X_test_MLP_'+model_to_load_filename.split('/')[-1].split('.')[0]+'.csv', sep='\t')
        pd.DataFrame(self.Y_test.to_numpy(),columns=['Y_test']).to_csv(self.folder_name + '/Y_test_MLP_'+model_to_load_filename.split('/')[-1].split('.')[0]+'.csv', sep='\t')
        pd.DataFrame(self.y_pred,columns=['y_pred']).to_csv(self.folder_name + '/Y_pred_MLP_'+model_to_load_filename.split('/')[-1].split('.')[0]+'.csv', sep='\t')
        self.y_pred_proba.to_csv(self.folder_name + '/Y_prob_MLP_'+model_to_load_filename.split('/')[-1].split('.')[0]+'.csv', sep='\t')
        
        return results, self.Y_test, self.y_pred, self.y_pred_proba

    def inference(self, model_to_load_filename):  # n_nodes, drop, act_functs):
        model = load_model(model_to_load_filename)

        testdata = DataIntegration(self.trset_name, self.trset_GO, self.trset_TF, self.trset_miRNA)
        self.X_test, self.Y_test = testdata.read_data()
        self.X_test = self.X_test[self.selected_features]

        self.y_pred = model.predict_classes(self.X_test)
        self.y_pred_proba = model.predict_proba(self.X_test)

        try:
            AUC_score = [metrics.roc_auc_score(self.Y_test, self.y_pred)]
        except:
            AUC_score = []
        results = {'learning rate': [self.learning_rate],
                   'model name': [model_to_load_filename],
                   'confusion matrix': [metrics.confusion_matrix(self.Y_test, self.y_pred)],
                   'accuracy': [metrics.accuracy_score(self.Y_test, self.y_pred)],
                   'precision': [metrics.precision_score(self.Y_test, self.y_pred)],
                   'recall': [metrics.recall_score(self.Y_test, self.y_pred)],
                   'f1 score': [metrics.f1_score(self.Y_test, self.y_pred)],
                   'AUC score': AUC_score}
        # print('testing results on model ' + model_to_load_filename +' :\n')
        # print(results)



        # pd.DataFrame(data=self.results).to_csv(
        #     self.folder_name + '/testing_MLP_Results_use_val_' + str(self.use_validation_set) + '_lr' + str(
        #         self.learning_rate) + '_drop' + str(drop) + str(l1) + str(l2) + str(l3) + str(l4) + str(a1) + str(
        #         a2) + str(a3) + str(a4) + '.csv', sep='\t')
        # pd.DataFrame(self.X_test, columns=['X_test']).to_csv(
        #     self.folder_name + '/X_test_MLP_' + model_to_load_filename.split('/')[-1].split('.')[0] + '.csv', sep='\t')
        # pd.DataFrame(self.Y_test.to_numpy(), columns=['Y_test']).to_csv(
        #     self.folder_name + '/Y_test_MLP_' + model_to_load_filename.split('/')[-1].split('.')[0] + '.csv', sep='\t')
        # pd.DataFrame(self.y_pred, columns=['y_pred']).to_csv(
        #     self.folder_name + '/Y_pred_MLP_' + model_to_load_filename.split('/')[-1].split('.')[0] + '.csv', sep='\t')
        # self.y_pred_proba.to_csv(
        #     self.folder_name + '/Y_prob_MLP_' + model_to_load_filename.split('/')[-1].split('.')[0] + '.csv', sep='\t')

        return results, self.X_test, self.y_pred, self.y_pred_proba

    def plot_results(self):
  
        print('AUC values:')
        for item in self.results.AUC:
            print('%.2f' % float(item))
        plt.figure()
        plt.plot([item[-1] for item in self.results.values if item[0] == 'train'],'bo', label='train_AUC')
        plt.plot([item[-1] for item in self.results.values if item[0] == 'test'],'ro', label='test_AUC')
        plt.legend()
        plt.ylim(0.5,1)
        plt.xticks(np.arange(len(self.results) / 2),
                   [','.join([str(item[4]), str(item[5]), str(item[6]), str(item[7])]) for item in self.results.values
                    if item[0] == 'train'], rotation=-30)
        plt.savefig(f'MLP_performances on {self.model_filename}_dropout{str(drop)}{l1}-{l2}-{l3}-{l4}_actfunc{a1}-{a2}-{a3}-{a4}.png')

