'''
Some functions to fit a random forest
'''

import sklearn.ensemble
import pandas
import progressbar

bar = progressbar.ProgressBar()

def test_max_features(max_features):
    if (max_features not in ['sqrt', 'auto', 'log2', None]):
        try:
            max_features = int(max_features)
        except ValueError:
            print("max_features has to be an integer or one of 'sqrt', 'auto', 'log2' or None.")
            raise
    return max_features

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
    importance = rf.estimators_[0].feature_importances_
    for est in bar(rf.estimators_[1:]):
        importance += est.feature_importances_
    importance = importance/rf.n_estimators
    d = {"kmer": kmers,
        "importance": importance}
    d = pandas.DataFrame(d)
    d = d.sort_values(by = "importance", ascending = 0)
    d = d.loc[d.importance > 0]
    return d
