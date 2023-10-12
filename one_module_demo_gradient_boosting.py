import pickle
import pandas as pd
import sklearn
import numpy as np
import scipy

# load the dic file
from utils import *

path = 'cleaned/module_dict.pickle'
with open(path, 'rb') as f:
    module_dict = pickle.load(f)

tissues = ['islet', 'liver', 'adipose', 'kidney', 'gastroc']
signs = ['SIGNED', 'UNSIGNED']
X_train, X_test, Y_train, Y_test = get_data('liver', 'SIGNED', 'yellow')

print(X_train.shape, X_test.shape, Y_train.shape, Y_test.shape)

# use sklearn gradient boosting regression method to train a model then predict
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GridSearchCV

# train a model
rds = np.random.RandomState(0)
model = GradientBoostingRegressor(random_state=0)
params = model.get_params()
print(params)

param_grid = {"max_depth": [3, 5, 7, 9],
              "subsample": [0.5, 0.8, 1],
              "n_estimators": [50, 100, 150],
              "learning_rate": [0.05, 0.1, 0.15, 0.2]
              }
Y_pred = {}
# gsh = GridSearchCV(estimator=model, param_grid=param_grid)
for name in Y_train.keys():
    print(name)
    gsh = GridSearchCV(estimator=model, param_grid=param_grid)
    gsh.fit(X_train, Y_train[name])
    # regressor = GradientBoostingRegressor(max_depth=3, n_estimators=100, random_state=0).fit(X_train, Y_train[name])
    # regressor.fit(X_train, Y_train[name])
    Y_pred[name] = gsh.predict(X_test)
    # # gsh.fit(X_train, Y_train[name])
    # # predict
    # print("s")
    # print(Y_train[name])
    # regressor.fit(X_train, Y_train[name])
    # print("done")
    # Y_pred[name] = regressor.predict(X_test)


# # visualize the gsh.cv_results_
# rank = gsh.cv_results_['rank_test_score']
# para = gsh.cv_results_['params']
# # print para based on rank
# for i in rank:
#     print(para[int(i) - 1])


# compute the average mean percentage error
def mean_percentage_error(y_true, y_pred):
    return np.mean(np.abs((y_true - y_pred) / y_true), axis=0) * 100


Y_pred = pd.DataFrame.from_dict(Y_pred)
per_err = mean_percentage_error(Y_test, Y_pred.to_numpy())
# compute the mean percentage error
per_err_mean = np.mean(per_err, axis=0)
print(per_err_mean)
# merge per_err and per_err_mean as pandas dataframe
per_err['mean'] = per_err_mean
print(per_err)
