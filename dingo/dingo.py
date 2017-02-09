'''
The main bits of dingo
'''


# system imports
import click

# local imports
from jellyfish import JellyFish
import input_reader

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
def main(input_file, ksize, hashsize, min_number, sreads, nbytes, single_end, force, outdir):
    # check that necessary software exists, otherwise quit
    jf = JellyFish()
    jf.exists()

    # check that outdir already exists, if NOT force, then quit

    # load input file --- check that paths exist otherwise quit
    data = input_reader.read(input_file)

    # create output folder structure --- if can't write quit

    # run jellyfish to identify all the kmers
    jf.count_all_mers(data, ksize, hashsize, min_number = min_number, simult_read = sreads, n_bytes = nbytes)
    # run jellyfish to count kmers in individual isolates
    jf.count_ind_mers(data, ksize, hashsize, min_number = min_number, simult_read = sreads, n_bytes = nbytes)
    # merge individual jellyfish results to generate our input matrix
    X = jf.join_counts(data)
    # run random forests to learn something
    pass

if __name__ == "__main__":
    main()
