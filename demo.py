import pickle
import pandas as pd
import sklearn
import numpy as np
import scipy

from multiprocessing import Manager
from multiprocessing import Process
# load the dic file

from utils import *

# knn regression model
# use sklearn knn regression method to train a model then predict
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

    
def run(tissue, sign,module,return_dict):
    # get data
    X_train, X_test, Y_train, Y_test = get_data(tissue, sign,module)

    # use only the last week glucose
    Y_train = Y_train.iloc[:,-1]
    Y_test = Y_test.iloc[:,-1]

    # train a model
    knn = KNeighborsRegressor(n_neighbors=20)
    knn.fit(X_train, Y_train)

    # predict
    Y_pred = knn.predict(X_test)

    # compute the average mean percentage error
    def mean_percentage_error(y_true, y_pred):
        return np.mean(np.abs((y_true - y_pred) / y_true)) * 100

    per_err = mean_percentage_error(Y_test, Y_pred)
    # compute the average mean squared error
    return_dict[(tissue, sign,module)] = per_err
    return


if __name__ == '__main__':
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
                p = Process(target=run, args=(tissue, sign,module,return_dict))
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
    
    print(return_dict)

    # now put return_dict into a dataframe
    df = pd.DataFrame.from_dict(return_dict, orient='index')
    # save the dataframe
    df.to_csv('knn.csv')