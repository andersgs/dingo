'''
Some auxiliary functions to check on, and run jellyfish
'''

import distutils.version
import os
import shlex
import subprocess
import sys
import pandas

jellyfish_min_version = "2.2.4"

def run(cmd):
    '''
    Process generator
    '''
    print("Running: {}".format(cmd), file = sys.stderr)
    p = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = p.communicate()
    return [out, err]

def read_table(sample_id):
    return pandas.read_csv("{}.txt".format(sample_id), delimiter = '\t', names = ['kmer', '{}'.format(sample_id)])

def join_tables(master, new_table, on):
    return pandas.merge(master, new_table, on = on, sort = False)

class JellyFish:
    def __init__(self):
        self.cmd = ''
        self.version = ''
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
    def __build_input(self, path):
        path = os.path.abspath(path)
        if path.endswith('.gz'):
            return " <(zcat {})".format(path)
        else:
            return " {}".format(path)
    def count_all_mers(self, tab, ksize, hash_size, output_file = 'allcount', min_number = 10, simult_read = 2, n_bytes = 1):
        cmd = self.cmd + ' count -s {} -m {} -F {} --out-counter-len {} -C -L {} -o {}'.format(hash_size, ksize, simult_read, n_bytes, min_number, output_file)
        for s in tab:
            cmd += self.__build_input(s[3])
        p = run(cmd)
        cmd = self.cmd + ' dump -o {0}.fa {0}'.format(output_file)
        p = run(cmd)
    def count_ind_mers(self, tab, ksize, hash_size, infile = 'allcount.fa', min_number = 10, simult_read = 2, n_bytes = 1):
        for s in tab:
            output_file = s[0]
            cmd = self.cmd + ' count -s {} -m {} -F {} --out-counter-len {} -C -o {}.jf --if {}'.format(hash_size, ksize, simult_read, n_bytes, output_file, infile) + self.__build_input(s[3])
            p = run(cmd)
            cmd = self.cmd + ' dump -ct -o {0}.txt {0}.jf'.format(output_file)
            p = run(cmd)
    def join_counts(self, tab):
        master = read_table(tab[0][0])
        for s in tab[1:]:
            master = join_tables(master, read_table(s[0]), on = 'kmer')
        return master
