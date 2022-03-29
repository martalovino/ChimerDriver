# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:08:19 2020

@author: vener

The code runs in three mode: build, train or test
For "build" the feature dataset are built from scratch starting from the filename.
Build mode has been built to support data in the form of Pegasus, DEEPrior or Oncofuse.
Anyway it supports any dataset with the required information specified below
(the column names must match exactly one of the three cases and the data has to be separated by tabs)
One must build the training set first and then any other testing/validation set because
the gene ontologies features will be defined only once that is when the training set is built.
After building a piece of word will be added to each dataset provided for each of the 4 families of features:
X_benign.csv -> X_benign_structfeat.csv, X_benign_TF.csv, X_benign_GO.csv, X_benign_miRNA.csv


Once you have built each of the dataset you intend to use as a training set and testing set
you can go on and train the model. If the training dataset is divided into multipe files
the algorithm will concatenate each of the sets, provided that you insert the list of filenames
If the filename is single, no concatenation will happen.
If the filename is in the form: X_benign.csv+X_oncogenic.csv the algorithm will merge
the two datasets into a single training or testing dataset according to the command.
For example a command line like the following:
>> python FeatureFusion.py train Pegasus_data Pegasus_data/X_chimerDB.txt+Pegasus_data/X_chimerDB_frameshifter.txt
+Pegasus_data/X_reactiveLN.txt Pegasus_data/X_benign.txt+Pegasus_data/X_oncogenic.txt . False forest all 0.0005 0.01 0.2
will use both X_benign and X_oncogenic as a test set, this will be possible by merging

the two datasets before proceeding onto the training.

For the testing is the following:
>> python FeatureFusion.py test DEEPrior_data DEEPrior_data/training_set.csv
 DEEPrior_data/test_set.csv DEEPrior_data/test_set_1.csv True subset_forest 5_GO_TF 0.0005 0.01 0.2
"""
import os.path

import pandas as pd

from ChimerDriver_tools import *
from pathlib import Path
import sys

if __name__ == '__main__':
  

    action = sys.argv[1]

    #  CREATING ONCOFUSE FILES
    try:
        dummy = pd.read_csv('Oncofuse/oncofuse_TICDB.csv',sep='\t')
    except FileNotFoundError:
        Path("Oncofuse/").mkdir(parents=True, exist_ok=True)    
        fusion_candidates = pd.read_excel('fusion_candidates.xlsx', sheet_name='fusion_candidates', engine='openpyxl')
        fusion_candidates = fusion_candidates.drop(columns=[n for n in list(fusion_candidates.columns) if n not in ['DB_NAME',"FPG5_NAME", "FPG3_NAME"]])
        fusion_candidates['Label']=fusion_candidates['DB_NAME']
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('TICDB', 1)
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('NORM', 0)
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('RTH', 0)
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('CHIMERDB2A', 1)
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('CHIMERDB2B', 1)
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('CHIMERDB2C', 1)            
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('NGS1', 1)
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('NGS2', 1)
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('REFSEQ', 0)
        fusion_candidates['Label'] = fusion_candidates['Label'].replace('CGC', 0)
        for db_name in list(set(fusion_candidates['DB_NAME'].to_list())):
            fusion_candidates[fusion_candidates['DB_NAME']==db_name].drop(columns=['DB_NAME']).to_csv('Oncofuse/oncofuse_'+db_name+'.csv',sep='\t')
    
    if action == 'build':

        #action,folder_name,filename_dataset,istrainset,label_dataset = 'build','dummy','dummy_data_like_DEEPrior.csv','True','N'
        #action,folder_name,filename_dataset,istrainset,label_dataset = 'build','dummy','dummy_data_like_Pegasus.csv','False','0'    
        folder_name,filename_dataset,whichset,label_dataset = sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]
        filename_out = filename_dataset[:-4].split('.')[-1] # preso tutto il nome esclusa l'estensione e selezionato l'ultima stringa dopo il punto
        filename_base = filename_out.split('/')[-1] #prendo solo il nome base senza la cartella di provenienza
        
        # create folder if it doesn't exist
        Path(folder_name+"/").mkdir(parents=True, exist_ok=True)
        
        
        ########## CREATING STRUCTURAL FEATURES ###########
        print('Creating structural features...')
        IF = initial_features(filename_dataset)
        IF.obtaining_partial_features(label_dataset,folder_name+"/"+filename_base+'_structfeat.csv')

        ######### CREATING FEATURES TRANSCRIPTION FACTORS ###########
        print('Creating transcription factor features...')
        TF = Transcription_factors_features(folder_name+"/"+filename_base+'_structfeat.csv')
        TF.build_TF_features(folder_name+"/"+filename_base+'_TF.csv')
        #######################################   
       
        ######### CREATING FEATURES GENE ONTOLOGIES ########### already saved 'GO_oncofuse_testset_usando_solo_train.csv' already saved
        # the GOs found in the training set will be the ones used for the rest of the sets
        print('Creating gene ontology features...')
        GO = Gene_ontologies_features(folder_name+"/"+filename_base+'_structfeat.csv', whichset)
        num_GOs = GO.build_GO_features(folder_name+"/"+filename_base+'_GO.csv')
        #######################################


        ######### CREATING FEATURES miRNA ########### already saved "gene_miRNA_matrix_oncofuse_trset"
        print('Creating microRNA features...')
        miRNA = miRNA_features(folder_name+"/"+filename_base+'_structfeat.csv')
        num_miRNA = miRNA.build_miRNA_features(folder_name+"/"+filename_base+'_miRNA.csv')
        #######################################

    
    elif action == 'rename_cols': #rename columns and build dataset to fit the algorithm requests
        case, filename_in, filename_out, name_col_label, name_col_version = sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6]
        df = pd.read_csv(filename_in, sep='\t')
        if name_col_label=='0' or name_col_label=='1':
            df['Label']=int(name_col_label)
            labels = df['Label'].to_list()
        else:
            labels = int(df[name_col_label]).to_list()   

        if name_col_version=='grch37' or name_col_version=='grch38': # you can simply write grch37 to apply this version to the whole dataset
            df['Version'] = name_col_version
            version = df['Version'].to_list()
        else:
            version = df[name_col_version].to_list()

        if case == 'Starfusion': # nomi geni, colonna label/valore label, version, breakpoints
            
            df[['GeneName5p','GeneEnsg5p']] = df['LeftGene'].str.split('^',expand=True,)
            df[['GeneName3p','GeneEnsg3p']] = df['RightGene'].str.split('^',expand=True,)
            df[['Chr5p','Breakpoint5p','Strand5p']] = df['LeftBreakpoint'].str.split(':',expand=True,)
            df[['Chr3p','Breakpoint3p','Strand3p']] = df['RightBreakpoint'].str.split(':',expand=True,)
            df['Chr5p'],df['Chr3p']=df['Chr5p'].str[3:], df['Chr3p'].str[3:]

            df_out = pd.DataFrame({ 'FusionPair':df['GeneName5p'].str.cat(df['GeneName3p'],sep="_"),
                                   'Label': labels,
                                   'Version':version,
                                   '5pCommonName':df['GeneName5p'],
                                   'Chr5p': df['Chr5p'],
                                   'Coord5p': df['Breakpoint5p'],
                                   '5pStrand': df['Strand5p'],
                                   '3pCommonName':df['GeneName3p'],
                                   'Chr3p': df['Chr3p'],
                                   'Coord3p': df['Breakpoint3p'],
                                   '3pStrand': df['Strand3p']
                           })
            df_out.to_csv(filename_out, sep='\t')
        elif case == 'DEEPrior': # names, bps, version, label, strands (come DEEPrior, strands potrebbe essere omesso)
            name_col_gene5p, name_col_gene3p, name_col_label, name_col_cr5p, name_col_cr3p, name_col_bp5p, name_col_bp3p, strand5p, strand3p = sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11],sys.argv[12],sys.argv[13],sys.argv[14],sys.argv[15]
            
            df_out = pd.DataFrame({ 'FusionPair':df[name_col_gene5p].str.cat(df[name_col_gene3p],sep="_"),
                           'Label': labels,
                           'Version':version,
                           'Chr5p': df[name_col_cr5p],
                           'Coord5p': df[name_col_bp5p],
                           '5pStrand': df[name_col_strand5p],
                           'Chr3p': df[name_col_cr3p],
                           'Coord3p': df[name_col_bp3p],
                           '3pStrand': df[name_col_strand3p]                          
            })

            df_out.to_csv(filename_out, sep='\t')
        elif case == 'Pegasus': # names, bp2, start, end, label, strands (come Pegasus)
            name_col_gene5p, name_col_gene3p, name_col_cr5p, name_col_cr3p, name_col_bp5p, name_col_bp3p, strand5p, strand3p = sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10],sys.argv[11],sys.argv[12],sys.argv[13],sys.argv[14]
            df_out = pd.DataFrame({ 'FusionPair':df[name_col_gene5p].str.cat(df[name_col_gene3p],sep="_"),
                           'Label': labels,
                           'Version':version,
                           'Chr1': df[name_col_cr5p],
                           'Gene_Breakpoint1': df[name_col_bp5p],
                           'Strand1': df[name_col_strand5p],
                           'Gene_Start1': df[name_col_start5p],
                           'Gene_End1': df[name_col_end5p],
                           'Chr2': df[name_col_cr3p],
                           'Gene_Breakpoint2': df[name_col_bp3p],
                           'Strand2': df[name_col_strand3p],
                           'Gene_Start2': df[name_col_start3p],
                           'Gene_End2': df[name_col_end3p]         
                               })
            df_out.to_csv(filename_out, sep='\t')



    
    elif action == 'cross_val_model' or action == 'train_test_model':
    
        folder_name, train_filename,test_filename,val_filename, use_validation_set, user_feat_sel,feat_set,threshold,lr,drop = sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],float(sys.argv[9]),float(sys.argv[10]),float(sys.argv[11])


        new_filenames = []        
        for filename in [train_filename,test_filename,val_filename]:
            #filename = 'dummy_data_like_DEEPrior.csv+dummy_data_like_Pegasus.csv'
            if len(filename)>1:
                JD = JoinDatasets(filename)
                if JD.isToMerge():
                    new_filenames.append(JD.MergeDatasets(folder_name)) # filename_out_base unisce i due filename precedenti, non include il .csv che viene aggiunto dalla stessa funzione
                else:
                    new_filenames.append(filename.split('/')[-1].split('.')[0]+'_')
            else:
                new_filenames.append(filename.split('/')[-1].split('.')[0]+'_')     
        train_filename_base,test_filename_base,val_filename_base = new_filenames[0],new_filenames[1],new_filenames[2]

    
        feat_selections_dict = {'5':['initial_features'], 'TF':['TF'],  '5_TF':['initial_features','TF'],  'miR':['miRNA'],    '5_miR':['initial_features', 'miRNA'],
                               'TF_miR':['TF','miRNA'],   '5_TF_miR':['initial_features', 'miRNA','TF'],  'GO':['GOs'], '5_GO':['initial_features','GOs'],  'GO_TF':['GOs','TF'],  'GO_miR':['GOs','miRNA'],
                               '5_GO_miR':['initial_features','GOs','miRNA'], '5_GO_TF':['initial_features', 'GOs', 'TF'],'5_GO_TF_miR':['initial_features', 'GOs', 'TF', 'miRNA']}
    
        ######### FEATURE SELECTION ########### on the train set
        if user_feat_sel == 'forest':
            # select now :
            FeatSel = FeatureSelection(filename_out=folder_name + '/' + 'feat_sel'+feat_set+'.txt')
            features_selected = FeatSel.Select_with_trees(filename1=folder_name + '/' + train_filename_base+'structfeat.csv',
                                                          filename2=folder_name + '/' + train_filename_base+'GO.csv',
                                                          filename3=folder_name + '/' + train_filename_base+'TF.csv',
                                                          filename4=folder_name + '/' + train_filename_base+'miRNA.csv',
                                                          thresh = threshold)
        elif user_feat_sel == 'subset':
            # select by subset of feature:
            FeatSel = FeatureSelection(filename_out=folder_name + '/' + 'feat_sel'+feat_set+'.txt')
            # possible features: 'initial_features', 'GOs', 'TF', 'miRNA'
            features_selected = FeatSel.Select_by_name(feat_selections_dict[feat_set],
                                                       filename1=folder_name + '/' + train_filename_base+'structfeat.csv',
                                                       filename2=folder_name + '/' + train_filename_base+'GO.csv',
                                                       filename3=folder_name + '/' + train_filename_base+'TF.csv',
                                                       filename4=folder_name + '/' + train_filename_base+'miRNA.csv')
     
        elif user_feat_sel == 'subset_forest':
            # select by subset of feature AND random forest
            FeatSel = FeatureSelection(filename_out=folder_name + '/' + 'feat_sel'+feat_set+'.txt')
            # possible features: 'initial_features', 'GOs', 'TF', 'miRNA'
    
            features_selected = FeatSel.Select_by_name_and_random_forest(feat_selections_dict[feat_set],
                                                          filename1=folder_name + '/' + train_filename_base+'structfeat.csv',
                                                          filename2=folder_name + '/' + train_filename_base+'GO.csv',
                                                          filename3=folder_name + '/' + train_filename_base+'TF.csv',
                                                          filename4=folder_name + '/' + train_filename_base+'miRNA.csv',
                                                          thresh = threshold) 
        elif user_feat_sel == 'load':
            # or if already saved just load :
            features_selected = []
            with open('feat_sel'+feat_set+'.txt',"r") as f:
                for row in f:
                    features_selected.append(row.split("\n")[0])
    
    
        #######################################

        ########### MLP ###############

        MultilayerPerceptron = MLP(folder_name, train_filename_base, test_filename_base, val_filename_base, use_validation_set,feat_sel=features_selected,learning_rate = lr,training_epochs = 5000,batch_size=32)

        if action == 'cross_val_model':
            print('now performing cross validation on 10 folds with nodes: 512,256,128,64 and activation functions: 4 sigmoids')
            cross_val_results, mean_acc, std_acc = MultilayerPerceptron.kfold_cross_validation(n_nodes=[512,256,128,64], act_functs=['sigmoid']*4, k=10, drop = drop)
            print('cross validation on 10 folds:\nmean of accuracy: %.3f,\nstd of accuracy = %.3f' %(mean_acc, std_acc))
            cross_val_results.to_csv(folder_name + '/cross_validation_MLP_'+user_feat_sel+'_'+feat_set+'_tresh'+str(threshold)+'_lr'+str(lr)+'_drop'+str(drop)+'.csv', sep='\t')

        elif action == 'train_test_model':
            # TRAINING & TESTING
            training_results, X_train, X_test, x_val, _, _, _ = MultilayerPerceptron.train_model(n_nodes=[512,256,128,64], drop= drop) 
            #results, Y_test, y_pred, y_pred_proba = MultilayerPerceptron.test_model([512,256,128,64], drop, ['tanh']*4)  # here you can change dropout if you want

            num_GO, num_miRNA, num_TF, num_init = 0, 0, 0, 0
            for feat in features_selected:
                if 'GO:' in feat:
                    num_GO+=1
                elif'miR-' in feat:
                    num_miRNA+=1
                elif 'let-' in feat:
                    num_miRNA+=1

                elif 'perc' in feat or 'same_strand' in feat or 'role' in feat:
                    num_init+=1
                else:
                    num_TF+=1
            print('%d feat tot\n%d GO\n%d miRNA\n%d TF\n%d initial features' %(len(features_selected),num_GO,num_miRNA,num_TF,num_init))

    elif action == 'load_test_model':
        folder_name,filename_dataset = sys.argv[2],sys.argv[3]
        whichset = 'test'
        label_dataset = 'N'
        filename_out = filename_dataset[:-4].split('.')[-1] # preso tutto il nome esclusa l'estensione e selezionato l'ultima stringa dopo il punto
        filename_base = filename_out.split('/')[-1] #prendo solo il nome base senza la cartella di provenienza
        
        # create folder if it doesn't exist
        Path(folder_name+"/").mkdir(parents=True, exist_ok=True)
        
        
        ########## CREATING STRUCTURAL FEATURES ###########
        print('Creating structural features...')
        IF = initial_features(filename_dataset)
        IF.obtaining_partial_features(label_dataset,folder_name+"/"+filename_base+'_structfeat.csv')

        ######### CREATING FEATURES TRANSCRIPTION FACTORS ###########
        print('Creating transcription factor features...')
        TF = Transcription_factors_features(folder_name+"/"+filename_base+'_structfeat.csv')
        TF.build_TF_features(folder_name+"/"+filename_base+'_TF.csv')
        #######################################   
       
        ######### CREATING FEATURES GENE ONTOLOGIES ########### already saved 'GO_oncofuse_testset_usando_solo_train.csv' already saved
        # the GOs found in the training set will be the ones used for the rest of the sets
        print('Creating gene ontology features...')
        GO = Gene_ontologies_features(folder_name+"/"+filename_base+'_structfeat.csv', whichset)
        num_GOs = GO.build_GO_features(folder_name+"/"+filename_base+'_GO.csv')
        #######################################


        ######### CREATING FEATURES miRNA ########### already saved "gene_miRNA_matrix_oncofuse_trset"
        print('Creating microRNA features...')
        miRNA = miRNA_features(folder_name+"/"+filename_base+'_structfeat.csv')
        num_miRNA = miRNA.build_miRNA_features(folder_name+"/"+filename_base+'_miRNA.csv')
        #######################################
    
    
    
    

        # LOADING & TESTING
        test_filename = sys.argv[3]
        feature_selected_to_load_filename = os.path.join(folder_name, 'feat_selall.txt')
        model_to_load_filename = 'best_model.h5'
        features_selected = []
        with open(feature_selected_to_load_filename,"r") as f:
            for row in f:
                features_selected.append(row.split("\n")[0])
        MultilayerPerceptron_load = MLP(folder_name, test_filename.split('/')[-1].split('.')[0] + '_', feat_sel=features_selected)

        print('Prediction in progress, please wait...')
        results, X_test, y_pred, y_pred_proba = MultilayerPerceptron_load.inference(model_to_load_filename)  # here you can change dropout if you want

        # adding predictions to the original file
        file_to_be_read = os.path.join(folder_name, test_filename.split('/')[-1].split('.')[0]) + '_structfeat.csv'
        final_output = pd.read_csv(file_to_be_read, sep='\t')
        # col = ['FusionPair','Version','Chr5p','Coord5p','5pStrand','5pCommonName','Chr3p','Coord3p','3pStrand','3pCommonName']
        # final_output = original[col]
        final_output['Class Prediction'] = [int(x) for x in list(y_pred)]
        final_output['Probability Prediction'] = [float(x) for x in list(y_pred_proba)]

        # save results to file!
        final_output.to_csv(os.path.join(folder_name, 'ChimerDriverPredictions.csv'), sep='\t')
        print('\n\nResults successfully saved at ' + os.path.join(folder_name, 'ChimerDriverPredictions.csv'))

