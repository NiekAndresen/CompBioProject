import numpy as np
import scipy.spatial.distance
import itertools
import time

def zero_one_loss(y_true, y_pred):
    return (np.sign(y_true) != np.sign(y_pred)).sum() / float(len(y_true))
    
#assumes classes 0 and 1
def false_discovery_rate(y_true, y_pred):
    return np.logical_and(y_pred, ~y_true).sum() / float(y_pred.sum())

def cv(X, y, method, parameters, nfolds=10, nrepetitions=5, loss_function=zero_one_loss):
    #see if two classes are present
    different_labels = set()
    for ll in y:
        different_labels.add(ll)
        if len(different_labels) > 1:
            break
    if len(different_labels) < 2:
        print("(XVAL) ERROR: only one class given")
        return None
    N = len(X)
    set_size = N // nfolds
    nof_param_sets = len(list(itertools.product(*parameters.values())))#needed for remaining time estimation
    param_sets = (dict(zip(parameters, x)) for x in itertools.product(*parameters.values()))
    loss = []
    best_loss = np.inf
    avg_fold_time = 0. #average time of past iterations
    total_folds = nof_param_sets * nrepetitions * nfolds #how many iterations will be run in total
    fold_count = 0 #counting the elapsed folds
    for setI,param_set in enumerate(param_sets):
        print("(XVAL) Testing parameter set %d of %d."%(setI+1,nof_param_sets))
        losses = np.zeros([nrepetitions, nfolds])
        for rep in range(nrepetitions):
            time_remaining = (total_folds-fold_count)*avg_fold_time
            print("(XVAL)   Repetition %3d of %3d. Remaining: %.2fs"%(rep+1, nrepetitions, time_remaining))
            #determine training and test set indices for each fold
            valid_folding = False
            while not valid_folding:
                test_indices = []
                train_labels = []
                fold_indices = np.random.choice(np.arange(N), N, replace=False)
                for fold in range(nfolds):
                    #the last training set is bigger to include all points at least once
                    max_index = N if fold==nfolds-1 else (fold+1)*set_size
                    test_indices += [fold_indices[fold*set_size:max_index]]
                    train_labels += [np.delete(y, test_indices[-1])]
                    #test if all folds contain samples from both classes to be able to train
                    different_labels = set()
                    valid_folding = False
                    for ll in train_labels[-1]:
                        different_labels.add(ll)
                        if len(different_labels) > 1:
                            valid_folding = True
                            break
                    if not valid_folding:
                        break
            #build sets, train and determine loss for each fold
            for i,fold in enumerate(range(nfolds)):
                start_time = time.time()
                Xtest = X[test_indices[i]]
                Xtrain = np.delete(X, test_indices[i], axis=0)
                mObj = method(**param_set)
                mObj.fit(Xtrain, train_labels[i])
                prediction = mObj.predict(Xtest)
                losses[rep,fold] = loss_function(y[test_indices[i]],prediction)
                fold_count += 1
                avg_fold_time = (avg_fold_time * (fold_count-1) + (time.time() - start_time)) / fold_count
        mloss = losses.mean()
        loss += [mloss]
        if mloss < best_loss:
            best_loss = mloss
            best_params = param_set
    methodObject = method(**best_params)
    methodObject.fit(X,y)
    methodObject.cvloss = best_loss
    return methodObject
