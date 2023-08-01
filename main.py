import pickle
import pandas as pd
import sklearn
import numpy as np
import scipy

from multiprocessing import Manager
from multiprocessing import Process

from sklearn.model_selection import GridSearchCV
# load the dic file

from utils import *

# knn regression model
# use sklearn knn regression method to train a model then predict
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import Lasso

    
def run(tissue, sign ,module , method, return_dict):
    # get data
    X_train, X_test, Y_train, Y_test = get_data(tissue, sign,module)

    # to control the randomness
    rds = np.random.RandomState(0)

    if method == 'knn':
    # train a model
        model = KNeighborsRegressor()
        n_neighbors = [5,10,20,30,40,50]
        metric = ['correlation','minkowski']
        param_grid = {"n_neighbors": n_neighbors,'metric':metric}

    if method == 'lasso':
        alpha= [0.1,1,10,100,1000,5000,10000]
        param_grid = {"alpha": alpha}
        model = Lasso(random_state=rds,selection='random',max_iter=3000,tol=0.0005)
    
    else:
        raise ValueError('method "{}" not supported'.format(method))
    
    gsh = GridSearchCV(estimator=model, param_grid=param_grid)
    gsh.fit(X_train, Y_train)

    # predict
    Y_pred = gsh.predict(X_test)

    # compute the average mean percentage error
    def mean_percentage_error(y_true, y_pred):
        return np.mean(np.abs((y_true - y_pred) / y_true)) * 100

    per_err = mean_percentage_error(Y_test, Y_pred)
    # compute the average percentage error of per_err
    
    # compute the mean percentage error
    per_err_mean = np.mean(per_err,axis=0)
    # merge per_err and per_err_mean as pandas dataframe
    per_err['mean'] = per_err_mean
    
    # compute the average mean squared error
    return_dict[(tissue, sign,module)] = per_err
    
    return

import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--method', type=str, help='ML method to use')
    method = parser.parse_args().method

    tissues = ['islet', 'liver', 'adipose', 'kidney', 'gastroc']
    signs = ['SIGNED', 'UNSIGNED']

    path = 'cleaned/module_dict.pickle'
    with open(path, 'rb') as f:
        module_dict = pickle.load(f)

    # run different module using python process
    process_lst = []
    return_dict = Manager().dict()

    for tissue in tissues:
        for sign in signs:
            for module in module_dict[(tissue,sign)].keys():
                p = Process(target=run, args=(tissue, sign, module, method, return_dict))
                process_lst.append(p)

    batch_size = 20
    # run the processes and get their return values
    for i in range(0, len(process_lst), batch_size):
        batch = process_lst[i:i+batch_size]
        for p in batch:
            p.start()
        for p in batch:
            p.join()
        print('finished batch {}'.format(i))

    # now put return_dict into a dataframe
    # now put dic into a dataframe
    columns = ['tissue', 'network', 'module', 'glucose','weight', 'insulin', 'triglyceride','mean']
    result = []
    for key in return_dict.keys():
        # parse the key, ('tissue', 'sign', 'module') is in key, 'glucose','weight', 'insulin', 'triglyceride','mean' is in value
        tissue, sign, module = key
        # get the value
        value = return_dict[key]
        # get the glucose, weight, insulin, triglyceride, mean
        glucose, weight, insulin, triglyceride, mean = value
        # put them into a list
        result.append([tissue, sign, module, glucose, weight, insulin, triglyceride, mean])
    df = pd.DataFrame(result, columns = columns)
    # save the dataframe
    df.to_csv('{}_result.csv'.format(method))