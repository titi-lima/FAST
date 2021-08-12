'''
This file is part of an ICSE'18 submission that is currently under review.
For more information visit: https://github.com/icse18-FAST/FAST.

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This software is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this source.  If not, see <http://www.gnu.org/licenses/>.
'''

import os
import re
import sys
import time
from pathlib import Path
#from struct import pack, unpack

# TODO: relative import. fix this latter.
sys.path.append(".")
from py import lsh



def parse_tc(path):
    parsed_tc = ''
    with open(path, 'r') as fin:
        line = fin.readline()
        while line:
            if line != "\n":
                newline = line.strip()
                parsed_tc += newline + " "
            line = fin.readline()
    parsed_tc = re.sub(' +', ' ', parsed_tc.strip())
    return parsed_tc


# TODO: store a single pickle file for each test
def preprocess_test_cases(inputdir, outputdir, singlefile=True):
    if singlefile:
        with open(JTeC_preproc_single, 'w') as fout, open(JTeC_preproc_map, 'w') as mapfile:
            for path in Path(JTeC_dir).rglob('*.java'):  # java files
                tc = parse_tc(path)  # parses the test files to remove line breaks and empty lines
                fout.write(tc + '\n')
                mapfile.write(str(path.relative_to(JTeC_dir)) + '\n')
    else:
        tcID = 1
        with open(JTeC_preproc_map, 'w') as mapfile:
            for path in Path(JTeC_dir).rglob('*.java'):  # java files
                outfile = os.path.join(JTeC_preproc_dir, '{}.txt'.format(tcID))
                with open(outfile, 'w') as fout:
                    tc = parse_tc(path)  # parses the test files to remove line breaks and empty lines
                    fout.write(tc)
                    mapfile.write(str(path.relative_to(JTeC_dir)) + '\n')
                    tcID += 1
        

# TODO: skip preprocessing if already done
def storeSignatures():
    mh_t = time.perf_counter()    
    with open(mapfileloc, 'w') as mapfile, open(sigfileloc, 'w') as sigfile:
        for path in Path(JTeC_dir).rglob('*.java'):  # java files
            tc = parse_tc(path)  # parses the test files to remove line breaks and empty lines
            tc_shingles = set()
            for i in range(len(tc) - k + 1):
                tc_shingles.add(hash(tc[i:i + k]))            
            sig = lsh.tcMinhashing((None, set(tc_shingles)), hashes)

            #print(sig)

            for hash_ in sig:
                #sigfile.write(repr(unpack('>d', hash_)[0]))
                sigfile.write(hash_)
                sigfile.write(" ")
            sigfile.write("\n")

            mapfile.write(str(path.relative_to(JTeC_dir)) + '\n')
    mh_time = time.perf_counter() - mh_t
    with open(sigtimefileloc, "w") as fout:
            fout.write(repr(mh_time))    



if __name__ == '__main__':
    #JTeC_dir = '/home/breno/research/JTEC/JTeC-Bundle/JTeC/'
    JTeC_dir = 'scalability/input/JTeC'
    mapfileloc = 'scalability/input/JTeC_map.txt'
    sigfileloc = 'scalability/input/JTeC.sig'
    sigtimefileloc = 'scalability/input/JTeC_sigtime.txt'
    #-----
    JTeC_preproc_dir = 'scalability/input/JTeC_preproc'
    JTeC_preproc_single = 'scalability/input/JTeC_preproc_all.txt'
    JTeC_preproc_map = 'scalability/input/JTeC_preproc_all_map.txt'

    if not os.path.exists(JTeC_preproc_dir):
        os.makedirs(JTeC_preproc_dir)

    # TODO: add all FAST parameters in a config file
    # FAST parameters
    k, n, r, b = 5, 10, 1, 10

    hashes = [lsh.hashFamily(i) for i in range(n)]

    #storeSignatures()
    
    preprocess_test_cases(JTeC_dir, JTeC_preproc_dir, singlefile=True)
    #preprocess_test_cases(JTeC_dir, JTeC_preproc_dir, singlefile=False)