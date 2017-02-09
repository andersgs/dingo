'''
The main bits of dingo
'''


# system imports
import click

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
### random forest options
@click.option("-n", "--n_trees", help = "Number of trees to grow", default = 10, show_default = True)
@click.option("-c", "--criterion", help = 'Criterion to decide on optimal split <entropy|gini>', default = "entropy", show_default = True)
@click.option("-m", "--max_features", help = "Maximum number of features to consider for each tree", default = "sqrt", show_default = True)
def main(input_file, ksize, hashsize, min_number, sreads, nbytes, single_end, force, outdir, n_trees, criterion, max_features):
    # check that necessary software exists, otherwise quit
    jf = JellyFish()
    jf.exists()

    # check that outdir already exists, if NOT force, then quit

    # load input file --- check that paths exist otherwise quit
    data = input_reader.read(input_file)

    # create response variable
    y = [s[1] for s in data]

    # create output folder structure --- if can't write quit

    # run jellyfish to identify all the kmers
    jf.count_all_mers(data, ksize, hashsize, min_number = min_number, simult_read = sreads, n_bytes = nbytes)
    # run jellyfish to count kmers in individual isolates
    jf.count_ind_mers(data, ksize, hashsize, min_number = min_number, simult_read = sreads, n_bytes = nbytes)
    # merge individual jellyfish results to generate our input matrix
    X,kmers = jf.join_counts(data)
    # run random forests to learn something
    learn = random_forest.learn(X = X, y = y, n_trees = n_trees, criterion = criterion, max_features = max_features)
    kmer_imp = random_forest.importance(learn, kmers)
    print(kmer_imp.head())
    pass

if __name__ == "__main__":
    main()
