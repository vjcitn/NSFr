
# coding: utf-8

# # Demonstration of NSforest with 10X genomics PBMC 10k antibody capture data
# 
# This is a slight modification to the notebook published at https://github.com/JCVenterInstitute/NSForest to work with a specific dataset.  The file Ab10k.tsv is derived from the 10x [PBMC 10K expression plus protein dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3). 
# 
# The 17 Ab features were used to produce a kmeans clustering with 9 clusters.  The file Ab10k.tsv was constructed from the SingleCellExperiment instance Ab10k.rda that is part of github.com/vjcitn/pcmp.

# In[1]:


#Necessary libraries

import numpy as np
import pandas as pd
import graphviz
import numexpr
import itertools
from subprocess import call

from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.metrics import fbeta_score
from sklearn.metrics import accuracy_score


# In[2]:


#Set Parameters

##Input data

dataFull = pd.read_table("%%TSVNAME%%", index_col = 0)


# In[3]:


dataFull.keys()


# In[4]:


#Creates dummy columns for one vs all Random Forest modeling
dataDummy = pd.get_dummies(dataFull, columns=["%%CLUSTVBLNAME%%"], prefix = "", prefix_sep = "")

#Creates matrix of cluster median expression values
medianValues = dataFull.groupby(by="%%CLUSTVBLNAME%%").median()
medianValues.to_csv('Function_medianValues.csv')     

#Finding the number of clusters and printing that to screen (sanity check)
PrecolNum = len(dataFull.columns)
PostcolNum = len(dataDummy.columns)
adjustedColumns = PrecolNum-1
clusters2Loop=PostcolNum-PrecolNum
#print clusters2Loop


####Random Forest parameters
rfTrees=%%NUMTREES%% #Number of trees
threads=%%NUMTHREADS%%     #Number of threads to use, -1 is the greedy option where it will take all available CPUs/RAM

####Filtering and ranking of genes from random forest parameters

Median_Expression_Level = 0
InformativeGenes = %%NUM_INF_GENES%% #How many top genes from the Random Forest ranked features will be evaluated for binariness 
Genes_to_testing = %%NUM_TOP_RANKED%%    #How many top genes ranked by binary score will be evaluated in permutations by fbeta-score (as the number increases the number of permutation rises exponentially!)

#### fbeta-score parameters                   

betaValue = 0.5 #Set values for fbeta weighting. 1 is default f-measure. close to zero is Precision, greater than 1 weights toward Recall


# In[5]:


#Function list


def randomForest(column, dataDummy,PrecolNum,rfTrees,threads):
    x_train = dataDummy[list(dataDummy.columns[0:PrecolNum-1])]
    names = dataDummy.columns[0:PrecolNum-1]
    y_train =dataDummy[column]
    rf = RandomForestClassifier(n_estimators = rfTrees, n_jobs = threads, random_state=123456)
    rf.fit(x_train, y_train)
    Ranked_Features=sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), names),reverse=True)
    return Ranked_Features

def rankInformative(Ranked_Features):
    RankedList = []
    midcounter = 0
    for x in Ranked_Features:
        midcounter +=1
        RankedList.append(x[1])
        rankedDict[column] = RankedList
        if midcounter==30:
            break
    return RankedList             

def negativeOut(x, column, medianValues,Median_Expression_Level):
    Positive_RankedList_Complete = []
    for i in x:
        if medianValues.loc[column, i] > Median_Expression_Level:
            print i
            print medianValues.loc[column, i]
            Positive_RankedList_Complete.append(i)
        else:
            print i
            print medianValues.loc[column, i]
            print "Is Right Out!"
    return Positive_RankedList_Complete

def binaryScore(Positive_RankedList_Complete, informativeGenes, medianValues, column):
    Positive_RankedList=list(Positive_RankedList_Complete[0:InformativeGenes])
    Median_RF_Subset=medianValues.loc[:, Positive_RankedList]
    Rescaled_Matrix=pd.DataFrame()
    
    for i in Positive_RankedList:
        Target_value=medianValues.loc[column, i]
        Rescaled_values=Median_RF_Subset[[i]].divide(Target_value)
        Rescaled_Matrix=pd.concat([Rescaled_Matrix,Rescaled_values],axis=1)
    difference_matrix=Rescaled_Matrix.apply(lambda x: 1-x, axis=1)
    difference_matrix_clean = difference_matrix.where(difference_matrix > 0, 0)
    ColumnSums=difference_matrix_clean.sum(0)
    rescaled = ColumnSums/clusters2Loop

    # Double sort so that for ties, the RF ranking prevails!     
    Ranked_Features_df=pd.DataFrame(Ranked_Features)
    Ranked_Features_df.rename(columns={1: 'Symbol'}, inplace=True)
    Ranked_Features_df_indexed=Ranked_Features_df.set_index("Symbol")
    rescaled_df=pd.DataFrame(rescaled)
    binaryAndinformation_Ranks=rescaled_df.join(Ranked_Features_df_indexed,lsuffix='_scaled', rsuffix='_informationGain')
    binaryAndinformation_Ranks.sort_values(by=['0_scaled','0_informationGain'],ascending= [False, False], inplace = True)
    Binary_ranked_Genes=binaryAndinformation_Ranks.index.tolist()
    Binary_RankedList=list(Binary_ranked_Genes[0:Genes_to_testing])
    Binary_scores=rescaled.to_dict()
    global Binary_store_DF 
    Binary_store_DF = Binary_store_DF.append(binaryAndinformation_Ranks)
    return Binary_RankedList

def DT_cutOffs(x, column):
    cut_dict = {}
    for i in x:
        filename=str(i)
        y_train =dataDummy[column]
        x_train = dataDummy[i]
        X = x_train[:, None]
        clf = tree.DecisionTreeClassifier(max_leaf_nodes=2)
        clf = clf.fit(X, y_train)
        threshold = clf.tree_.threshold
        cut_dict[i]=threshold[0]
    return cut_dict

def queryGenerator(x, cut_dict):
    queryList = []
    for i in x:
        str1 = i
        current_value = cut_dict.get(str1)
        queryString1 = str(str1)+'>='+ str(current_value)
        queryList.append(queryString1)
    return queryList

def permutor(x):
    binarylist2 = x
    combs = []
    for i in xrange(1, len(x)+1):
        els = [list(x) for x in itertools.permutations(binarylist2, i)]
        combs.extend(els)
    return combs

def fbetaTest(x, column,testArray, betaValue):
    fbeta_dict = {}
    for list in x:
        testArray['y_pred']= 0
        betaQuery = '&'.join(list)
        Ineq1=dataFull.query(betaQuery)
        testList=Ineq1.index.tolist()
        testArray.loc[testList, 'y_pred'] = 1
        f1 = fbeta_score(testArray['y_true'], testArray['y_pred'], average= 'binary', beta=betaValue)        
        dictName = column+"&"+betaQuery
        fbeta_dict[dictName] = f1 
    return fbeta_dict


# In[10]:


#Core analysis 
rankedDict =  {}  ###gives us the top ten features from RF
f1_store_1D = {}
Binary_score_store_DF=pd.DataFrame()
DT_cutoffs_store={}

for column in dataDummy.columns[PrecolNum-1:PostcolNum]:
       
        ## Run Random Forest and get a ranked list 
        Ranked_Features= randomForest(column, dataDummy, PrecolNum, rfTrees, threads)
        RankedList = rankInformative(Ranked_Features)
        
        ## Setup testArray for f-beta evaluation
        testArray = dataDummy[[column]]
        testArray.columns = ['y_true']
        
        #Rerank according to expression level and binary score
        Positive_RankedList_Complete = negativeOut(RankedList, column, medianValues, Median_Expression_Level) 
        Binary_store_DF = pd.DataFrame()
        Binary_RankedList = binaryScore(Positive_RankedList_Complete, InformativeGenes, medianValues, column)
        
        Binary_score_store_DF_extra = Binary_store_DF.assign(clusterName = column)
        #print Binary_score_store_DF_extra
        Binary_score_store_DF = Binary_score_store_DF.append(Binary_score_store_DF_extra)

        
        #Get expression cutoffs for f-beta testing
        cut_dict = DT_cutOffs(Binary_RankedList,column)
        DT_cutoffs_store[column]=cut_dict
        
        #Generate expression queries and run those queries using fbetaTest() function
        queryInequalities=queryGenerator(Binary_RankedList, cut_dict)
        FullpermutationList = permutor(queryInequalities)
        #print len(FullpermutationList)
        f1_store = fbetaTest(FullpermutationList, column, testArray, betaValue)
        f1_store_1D.update(f1_store)
       
    
    
#Report generation and cleanup         
        
f1_store_1D_df = pd.DataFrame() #F1 store gives all results.
f1_store_1D_df = pd.DataFrame.from_dict(f1_store_1D, orient='index', dtype=str)
f1_store_1D_df.columns = ["f-measure"]
f1_store_1D_df['markerCount'] = f1_store_1D_df.index.str.count('&')
f1_store_1D_df.reset_index(level=f1_store_1D_df.index.names, inplace=True)

f1_store_1D_df_done= f1_store_1D_df['index'].apply(lambda x: pd.Series(x.split('&')))

NSForest_Results_Table=f1_store_1D_df.join(f1_store_1D_df_done)

NSForest_Results_Table_Fin = pd.DataFrame()
NSForest_Results_Table_Fin = NSForest_Results_Table[NSForest_Results_Table.columns[0:4]]

for i, col in enumerate(NSForest_Results_Table.columns[4:11]):
    splitResults= NSForest_Results_Table[col].astype(str).apply(lambda x: pd.Series(x.split('>='))) 
    firstOnly = splitResults[0]
    Ascolumn = firstOnly.to_frame()
    Ascolumn.columns = [col]
    NSForest_Results_Table_Fin = NSForest_Results_Table_Fin.join(Ascolumn)
    
    
NSForest_Results_Table_Fin.rename(columns={0:'clusterName'},inplace=True) #rename columns by position
NSForest_Results_Table_Fin.sort_values(by=['clusterName','f-measure','markerCount'],ascending= [True, False, True], inplace = True)

#Write outs
Binary_score_store_DF.to_csv('Binary_scores_Supplmental_results.csv')
NSForest_Results_Table_Fin.to_csv('NS-Forest_v2_results.csv')       


#Subsets of full results 
max_grouped = NSForest_Results_Table_Fin.groupby(by="clusterName")["f-measure"].max()
max_grouped.df=pd.DataFrame(max_grouped)
max_grouped.df.to_csv('NSForest_v2_maxF-scores.csv')

NSForest_Results_Table_Fin["f-measureRank"] = NSForest_Results_Table_Fin.groupby(by="clusterName")["f-measure"].rank(ascending=False)
topResults = NSForest_Results_Table_Fin["f-measureRank"] < 50
NSForest_Results_Table_top = NSForest_Results_Table_Fin[topResults]
NSForest_Results_Table_top.to_csv('NSForest_v2_topResults.csv') 


    
        

