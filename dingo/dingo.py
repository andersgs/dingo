'''
The main bits of dingo
'''


# system imports
import click
import os
import sys

# local imports
from jellyfish import JellyFish
import input_reader
import random_forest

@click.command()
@click.option("-k", "--ksize", help = "Kmer size to search for.", default = 31, show_default = True)
@click.option("-a", "--hashsize", help = "Hash size", default = '5M', show_default = True)
@click.option("--min_number", help = "Minimum number of kmer observations to count.", default = 10, show_default = True)
@click.option("--sreads", help = "Number of files to read simultaneously by jellyfish", default = 2, show_default = True)
@click.option("--nbytes", help = "Size of the number used for counting kmers", default = 1, show_default = True)
@click.option("-s", "--single_end", is_flag = True, default = False, help = "Data is single end")
@click.option("-f", "--force", is_flag = True, default = False, help = "Write over previous analysis")
@click.option("-o", "--outdir", help = "Output folder")
@click.option("-i", "--input_file", help = "Input file.")
@click.option("--kmer_fa", help = "A FASTA of kmers to count", default = 'allcount.fa')
@click.option("-t", "--threads", help = "Number of threads to run Jellyfish", default = 16, show_default = True)
@click.option("-p","pickled_matrix", help = "Use a pickeled matrix", default = ".kmer_table.pickle", show_default = True)
### random forest options
@click.option("-n", "--n_trees", help = "Number of trees to grow", default = 10, show_default = True)
@click.option("-c", "--criterion", help = 'Criterion to decide on optimal split <entropy|gini>', default = "entropy", show_default = True)
@click.option("-m", "--max_features", help = "Maximum number of features to consider for each tree", default = "sqrt", show_default = True)
def main(input_file, ksize, hashsize, min_number, sreads, nbytes, single_end, force, outdir, threads, n_trees, criterion, max_features, kmer_fa, pickled_matrix):
    # check that necessary software exists, otherwise quit
    jf = JellyFish()
    jf.exists()

    # do some parameter checking
    max_features = random_forest.test_max_features(max_features)

    # check that outdir already exists, if NOT force, then quit

    # load input file --- check that paths exist otherwise quit
    data = input_reader.read(input_file)

    # create response variable
    y = [s[1] for s in data]

    # create output folder structure --- if can't write quit

    # run jellyfish to identify all the kmers
    # only run it if necessary
    if (os.path.exists(pickled_matrix) and not force):
        print("Found a pickled kmer matrix, going to use it...")
        X,kmers = jf.load_kmertable(pickle_file = pickled_matrix)
    else:
        if (kmer_fa == None or not os.path.isfile(kmer_fa)):
            jf.count_all_mers(data, ksize, hashsize, threads = threads, min_number = min_number, simult_read = sreads, n_bytes = nbytes)
        else:
            print("Found {}, so skipping counting kmers across all samples" .format(kmer_fa), file = sys.stderr)
        # run jellyfish to count kmers in individual isolates
        jf.count_ind_mers(data, ksize, hashsize, threads = threads, min_number = min_number, simult_read = sreads, n_bytes = nbytes)
        # merge individual jellyfish results to generate our input matrix
        print("Generating kmer table...", file = sys.stderr)
        X,kmers = jf.join_counts(data)
    print("Learning about the kmers...", file = sys.stderr)
    # run random forests to learn something
    learn = random_forest.learn(X = X, y = y, n_trees = n_trees, criterion = criterion, max_features = max_features)
    print("Computing importance of kmers...", file = sys.stderr)
    kmer_imp = random_forest.importance(learn, kmers)
    print("Making predictions...", file = sys.stderr)
    print(learn.predict(X), file = sys.stderr)
    #print(learn.predict_log_proba(X), file = sys.stderr)
    print(kmer_imp.head(), file = sys.stderr)
    kmer_imp.to_csv("junk.csv")
    pass

if __name__ == "__main__":
    main()
