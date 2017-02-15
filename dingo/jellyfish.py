'''
Some auxiliary functions to check on, and run jellyfish
'''

import distutils.version
import os
import shlex
import subprocess
import sys
import pandas
import numpy as np

jellyfish_min_version = "2.2.4"

def run(cmd):
    '''
    Process generator
    '''
    print("Running: {}".format(cmd), file = sys.stderr)
    #p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    #out, err = p.communicate()
    p = subprocess.check_call(cmd, shell = True)
    return 0

def read_table(sample_id, make_binary = True):
    '''
    Quickly read a table of kmer counts

    With make_binary, the output is always 0,1. But that can change later
    '''
    tab = pandas.read_csv("{}.txt".format(sample_id), delimiter = '\t', names = ['kmer', '{}'.format(sample_id)])
    if make_binary:
        tab[sample_id] = np.where(tab[sample_id] > 0, 1, 0)
    return tab

def join_tables(master, new_table, on):
    '''
    Quickly join two tables on column ON
    '''
    return pandas.merge(master, new_table, on = on, sort = False)

def find_var_rows(row):
    return np.count_nonzero(row) < len(row)

class JellyFish:
    def __init__(self, force = False):
        self.cmd = ''
        self.version = ''
        self.force = force
    def exists(self):
        try:
            p = subprocess.Popen(shlex.split('jellyfish --version'), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            out,err = p.communicate()
            cmd,version = out.split()
            self.cmd = cmd.decode()
            self.version = version.decode()
            if distutils.version.LooseVersion(self.version) >= distutils.version.LooseVersion(jellyfish_min_version):
                print("Found jellyfish version {}.... OK!".format(self.version))
            else:
                raise ValueError("Version of jellyfish found is less than {}, please update to run dingo".format(min_version))
        except ValueError:
            print("Did not find jellyfish on the path.")
    def __build_input(self, path, clear = False):
        if clear:
            gen = open("generator.txt", 'w')
        else:
            gen = open("generator.txt", 'a')
        path = os.path.abspath(path)
        if path.endswith('.gz'):
            gen.write("gunzip -c {}\n".format(path))
        else:
            gen.write("cat {}\n".format(path))
        gen.close()
        return
    def count_all_mers(self, tab, ksize, hash_size, threads = 16, output_file = 'allcount', min_number = 10, simult_read = 2, n_bytes = 1):
        cmd = self.cmd + ' count -s {} -m {} -G {} --out-counter-len {} -C -L {} -o {} -g {} -t {}'.format(hash_size, ksize, simult_read, n_bytes, min_number, output_file, 'generator.txt', threads)
        for s in tab:
            self.__build_input(s[3])
        p = run(cmd)
        cmd = self.cmd + ' dump -o {0}.fa {0}'.format(output_file)
        p = run(cmd)
    def count_ind_mers(self, tab, ksize, hash_size, threads = 16, infile = 'allcount.fa', min_number = 10, simult_read = 2, n_bytes = 1):
        for s in tab:
            output_file = s[0]
            if os.path.exists("{}.txt".format(output_file)) and not self.force:
                print("File {}.txt already exists... Skipping kmer counting!".format(output_file))
            else:
                self.__build_input(s[3], clear = True)
                cmd = self.cmd + ' count -s {} -m {} -G {} --out-counter-len {} -C -o {}.jf --if {} -g {} -t {}'.format(hash_size, ksize, simult_read, n_bytes, output_file, infile, "generator.txt", threads)
                p = run(cmd)
                cmd = self.cmd + ' dump -ct -o {0}.txt {0}.jf'.format(output_file)
                p = run(cmd)
    def join_counts(self, tab):
        master = read_table(tab[0][0])
        for s in tab[1:]:
            master = join_tables(master, read_table(s[0]), on = 'kmer')
        master = master.loc[master.apply(find_var_rows, axis = 1),] # remove kmers present in all samples
        master = master.transpose()
        print("Found {} variable kmers.".format(master.shape[1]), file = sys.stderr)
        master = master.values
        return [master[1:], master[0]]
