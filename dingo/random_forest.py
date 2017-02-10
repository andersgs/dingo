'''
Some functions to fit a random forest
'''

import sklearn.ensemble
import pandas

def learn(X,y, n_trees = 10, criterion = 'entropy', max_features = "sqrt", max_depth = None, min_samples_split = 2, min_samples_leaf = 1, min_weight_fraction_leaf = 0, max_leaf_nodes = None, min_impurity_split = 1e-7, boostrapt = True, oob_score = False, n_jobs = 10, random_state = None, warm_start = False, class_weight = 'balanced_subsample'):
    rf = sklearn.ensemble.RandomForestClassifier(n_estimators = n_trees, \
                                                criterion = criterion, \
                                                max_features = max_features, \
                                                max_depth = max_depth, \
                                                min_samples_split = min_samples_split, \
                                                min_samples_leaf = min_samples_leaf, \
                                                min_weight_fraction_leaf = min_weight_fraction_leaf, \
                                                max_leaf_nodes = max_leaf_nodes, \
                                                min_impurity_split = min_impurity_split, \
                                                oob_score = oob_score, \
                                                n_jobs = n_jobs, \
                                                random_state = random_state, \
                                                warm_start = warm_start, \
                                                class_weight = class_weight \
                                                )
    rf.fit(X, y)
    return rf

def importance(rf, kmers):
    d = {"kmer": kmers,
        "importance": rf.feature_importances_}
    d = pandas.DataFrame(d)
    d.sort_values(by = "importance", ascending = 0)
    return d
