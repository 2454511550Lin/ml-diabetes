# load the rna expression data from the csv file
import pandas as pd
import pickle
import numpy as np
import scipy

# load the genes id
def get_data(tissue,sign,module):
    
    path = 'cleaned/module_dict.pickle'
    with open(path, 'rb') as f:
        module_dict = pickle.load(f)

    path = 'cleaned/{}_rna.csv'.format(tissue)
    # load the csv file
    df = pd.read_csv(path, index_col=0)

    # gat the data with genes as the columns
    
    if module == 'all':
        X = df
    else:
        # get the target module genes ids
        genes = module_dict[(tissue,sign)][module]
        # check that all the genes exist in df column index
        assert all([gene in df.columns for gene in genes])
        X = df[genes]

    # load the glucose data
    traits = ['glucose','weight', 'insulin', 'triglyceride']

    Y = pd.DataFrame(columns = traits)
    for i in range(len(traits)):
        trait = traits[i]
        path = 'cleaned/clinic_traits/{}.csv'.format(trait)
        temp = pd.read_csv(path, index_col=0)['wk10']
        Y[trait] = temp
    
    # only keep the mouse with id appers in X
    Y = Y.loc[X.index]

    # load the train/test split 
    test_index = pickle.load(open('cleaned/train_test_split.pickle', 'rb'))
    
    # get the train and test data
    
    X_test = X.loc[test_index]
    Y_test = Y.loc[test_index]

    # compute the remaining index as training index
    train_index = X.index.difference(test_index)
    X_train = X.loc[train_index]
    Y_train = Y.loc[train_index]

    return X_train, X_test, Y_train, Y_test

# compute the column-wise spearman correlation 
def spearman(y_true, y_pred):
    corr = []
    for i in range(4):
        spearman_corr = scipy.stats.spearmanr(y_true.iloc[:,i], y_pred[:,i])
        corr.append(spearman_corr[0])
    return corr

def fill_miss_rna(tissue):

    path = 'raw_data\F2expr.ns.500mice.{}.WGCNA.csv'.format(tissue)
    # read the csv file
    rna = pd.read_csv(path, index_col=0)

    # flip column and row
    rna = rna.transpose()

    # get rid of the first row 
    rna = rna.iloc[1:]

    # make the first colume as the index
    # rna.index = rna.iloc[:,0]

    # change the dtype to float
    rna = rna.astype(float)

    rna_corr = rna.transpose().corr()

    rank = 2

    while rna.isnull().sum().sum() != 0:
        print("Rank: {}, {} rows still have missing value".format(rank, rna.isnull().sum().sum()))
        missing_row_indices = rna[rna.isnull().any(axis=1)].index
        # for each row, find the index of the row that has the highest correlation, excluding itself
        most_corr_row = rna_corr.apply(lambda x: x.nlargest(rank).index.tolist()[-1], axis=1)
        for i in missing_row_indices:
        # fill in just the missing value
            rna.loc[i] = rna.loc[i].fillna(rna.loc[most_corr_row[i]])
        rank += 1

    # save the rna to csv file
    rna.to_csv("./cleaned/{}_rna.csv".format(tissue))


import pandas as pd
import numpy as np
def fill_missing_glucose():
    # load in the clinic phenotypes
    clinic = pd.read_csv('raw_data/(F2)B6_BTBR OB mice_clinic_traits.csv', index_col=0)
    
    # load in rna expression data
    tissues = ['islet', 'liver', 'adipose', 'kidney', 'gastroc']
    indices = []
    for tissue in tissues:
        rna = pd.read_csv('cleaned/{}_rna.csv'.format(tissue), index_col=0)
        indices.append(rna.index)
    
    # get the union of all the indices
    rna_index = indices[0]
    for i in range(1, len(indices)):
        rna_index = rna_index.union(indices[i]) 
    
    glucose_name_tags = ["Unnamed: 7","Unnamed: 13","Unnamed: 19","Unnamed: 25"]
    
    # get the mouse id, conver it to string
    mouse_id = rna_index.to_list()
    mouse_id_clinic = clinic.index.to_list()
    
    # get rid of 'Mouse' at the beginning of each element of mouse id
    mouse_id = [i[5:] for i in mouse_id] 
    
    glucose = np.zeros((len(mouse_id), 4))
    # use mouse id and name tag to get the glucose level
    for i in range(len(mouse_id)):
        if mouse_id[i] not in mouse_id_clinic:
            print("mouse id {} not in clinic data".format(mouse_id[i]))
            continue
        for j in range(4):
            glucose[i,j] = clinic.loc[mouse_id[i], glucose_name_tags[j]]
    
    # check if there is missing value
    print("There are {} missing values in glucose".format(np.isnan(glucose).sum()))
    
    # for each missing value, fill it in using linear regression
    for i in range(len(glucose)):
        for j in range(4):
            if np.isnan(glucose[i,j]):
                # get the index of the mice that has no missing value
                index = np.where(~np.isnan(glucose[:,j]))[0]
                # get the glucose level of mice that has no missing value
                glucose_level = glucose[index,j]
                # use linear regression to predict the missing value
                glucose[i,j] = np.poly1d(np.polyfit(index, glucose_level, 1))(i)
    
    print("After regression, there are {} missing values in glucose".format(np.isnan(glucose).sum()))
    
    # save the glucose level to csv file
    glucose = pd.DataFrame(glucose, index=rna_index, columns=['wk4', 'wk6', 'wk8', 'wk10'])
    glucose.to_csv("./cleaned/glucose.csv")


def fill_missing_clinic_traits():
    # load in the clinic phenotypes
    clinic = pd.read_csv('raw_data/(F2)B6_BTBR OB mice_clinic_traits.csv', index_col=0)
    
    # load in rna expression data
    tissues = ['islet', 'liver', 'adipose', 'kidney', 'gastroc']
    indices = []
    for tissue in tissues:
        rna = pd.read_csv('cleaned/{}_rna.csv'.format(tissue), index_col=0)
        indices.append(rna.index)
    
    # get the union of all the indices
    rna_index = indices[0]
    for i in range(1, len(indices)):
        rna_index = rna_index.union(indices[i]) 
    
    # stands for weight, glucose, insulin, TRIGLYCERIDE
    name_tags = ["Unnamed: 23","Unnamed: 25","Unnamed: 26","Unnamed: 27"] 
    
    # get the mouse id, conver it to string
    mouse_id = rna_index.to_list()
    mouse_id_clinic = clinic.index.to_list()
    
    # get rid of 'Mouse' at the beginning of each element of mouse id
    mouse_id = [i[5:] for i in mouse_id] 
    
    glucose = np.zeros((len(mouse_id), 4))
    # use mouse id and name tag to get the traits
    for i in range(len(mouse_id)):
        if mouse_id[i] not in mouse_id_clinic:
            print("mouse id {} not in clinic data".format(mouse_id[i]))
            continue
        for j in range(4):
            glucose[i,j] = clinic.loc[mouse_id[i], name_tags[j]]
    
    # check if there is missing value
    print("There are {} missing values in glucose".format(np.isnan(glucose).sum()))
    
    # for each missing value, fill it in using linear regression
    for i in range(len(glucose)):
        for j in range(4):
            if np.isnan(glucose[i,j]):
                # get the index of the mice that has no missing value
                index = np.where(~np.isnan(glucose[:,j]))[0]
                # get the glucose level of mice that has no missing value
                glucose_level = glucose[index,j]
                # use linear regression to predict the missing value
                glucose[i,j] = np.poly1d(np.polyfit(index, glucose_level, 1))(i)
    
    print("After regression, there are {} missing values in glucose".format(np.isnan(glucose).sum()))
    
    # save the glucose level to csv file
    glucose = pd.DataFrame(glucose, index=rna_index, columns=['wk4', 'wk6', 'wk8', 'wk10'])
    glucose.to_csv("./cleaned/glucose.csv")
    