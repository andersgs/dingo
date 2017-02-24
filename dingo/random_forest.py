'''
Some functions to fit a random forest
'''

import multiprocessing
import numpy
import sklearn.ensemble
import pandas
import progressbar
import sys

bar = progressbar.ProgressBar()

def test_max_features(max_features):
    if (max_features not in ['sqrt', 'auto', 'log2', None]):
        try:
            max_features = int(max_features)
        except ValueError:
            print("max_features has to be an integer or one of 'sqrt', 'auto', 'log2' or None.")
            raise
    return max_features

def learn(X,y, n_trees = 10, criterion = 'entropy', max_features = "sqrt", max_depth = None, min_samples_split = 2, min_samples_leaf = 1, min_weight_fraction_leaf = 0, max_leaf_nodes = None, min_impurity_split = 1e-7, bootstrap = False, oob_score = False, n_jobs = 10, random_state = None, warm_start = True, class_weight = None):
    rf = sklearn.ensemble.RandomForestClassifier(n_estimators = n_trees, \
                                                criterion = criterion, \
                                                max_features = max_features, \
                                                max_depth = max_depth, \
                                                min_samples_split = min_samples_split, \
                                                min_samples_leaf = min_samples_leaf, \
                                                min_weight_fraction_leaf = min_weight_fraction_leaf, \
                                                max_leaf_nodes = max_leaf_nodes, \
                                                min_impurity_split = min_impurity_split, \
                                                bootstrap = bootstrap, \
                                                oob_score = oob_score, \
                                                n_jobs = n_jobs, \
                                                random_state = random_state, \
                                                warm_start = warm_start, \
                                                class_weight = class_weight, \
                                                verbose = 1
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
    return d

def chunks(l, n):
  for i in iter(range(0, len(l), n)):
    yield l[i:i+n]


def map_importance(estimators):
    importance = estimators[0].feature_importances_
    for est in estimators[1:]:
        importance += est.feature_importances_
    return importance

def parallel_importance(rf, kmers, n_workers = 8, start = 0):
    pool = multiprocessing.Pool(processes=n_workers,)
    total_jobs = len(rf.estimators_[start:])
    n_chunks = int(total_jobs / n_workers)
    if n_chunks < 1:
        n_chunks = 1
    print("Number of chunks: {}".format(n_chunks), file = sys.stderr)
    partitioned_estimators = list(chunks(rf.estimators_[start:], n_chunks))
    importance_list = pool.map(map_importance, partitioned_estimators)
    importance = numpy.sum(importance_list, axis = 0)/rf.n_estimators
    d = {"kmer": kmers,
        "importance": importance}
    pool.close()
    return d

def incremental_learn(rf, old_importance, X, y, n_tries = 50, n_inc = 5, n_workers = 8):
    print("Incrementing model...", file = sys.stderr)
    for i in range(n_tries):
        print("Try number {}...".format(i+1), file  = sys.stderr)
        cur_estimators = rf.n_estimators
        rf.n_estimators = cur_estimators + n_inc
        rf.fit(X,y)
        print("Current number of trees: {}...".format(cur_estimators))
        print("New total estimators: {}...".format(rf.n_estimators))
        incremental_importance = parallel_importance(rf, kmers = old_importance['kmer'], n_workers = n_workers, start = cur_estimators)
        new_batch_importance = incremental_importance['importance'] + (cur_estimators/rf.n_estimators) * old_importance['importance']
        # new_importance_check = rf.feature_importances_
        # print("Testing batch is approximately the same expected: {}.".format(numpy.testing.assert_allclose(new_batch_importance, new_importance_check)))
        print("Correlation kmer importance between runs: {}".format(numpy.corrcoef(new_batch_importance, old_importance['importance'])[0,1]), file = sys.stderr)
        old_importance['importance'] = new_batch_importance
    return [rf,old_importance]

def rank_importance(importance):
    d = pandas.DataFrame(importance)
    d = d.sort_values(by = "importance", ascending = 0)
    d = d.loc[d.importance > 0]
    return d
