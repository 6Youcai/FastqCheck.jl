#!/usr/bin/env python3

import os
import sys
import gzip
import numpy as np

# fqcheck
def fqcheck(infile, read_length = 151):
    # open file
    if infile.endswith(".gz"):
        hand = gzip.open(infile)
    else:
        # hand = open(infile)
        print("Only support gzip file now")
        sys.argv[-1]
    # init array
    m = (read_length, 48)
    out_array = np.zeros(m, 'int64')
    order_dict = {"A":0, "T":1, "C":2, "G":3, "N":4}
    reads = 0
    for line1 in hand:
        line2 = hand.readline().decode()[:-1] # trim the \n
        _, line4 = hand.readline(), hand.readline().decode()
        reads += 1
        Length = len(line2)
        for i in range(Length):
            base = line2[i]
            base_in_array = order_dict[base]
            out_array[i, base_in_array] += 1
            qual = ord(line4[i]) - 33 # phred 33
            qual_in_array = qual + 5 # qual 0 is behind ATCGN
            out_array[i, qual_in_array] += 1
    hand.close()
    return(out_array, reads)

def ration(a, t, c, g, n):
    Sum = a+t+c+g+n
    ra = a/Sum
    rt = t/Sum
    rc = c/Sum
    rg = g/Sum
    rn = n/Sum
    return(ra, rt, rc, rg, rn)

def error_rate(in_array):
    all_base = np.sum(in_array)
    sum_qual = 0
    for i in range(43):
        tmp = np.sum(in_array[:, i])
        sum_qual += tmp * i
    qual = sum_qual/all_base
    print(qual)
    E = 10**(-1*qual/10)
    return(E)

def printArray(in_array, reads, outfile):
    f = open(outfile, "w")
    all_base = np.sum(in_array[:, :5])
    q20 = np.sum(in_array[:, 25:])/all_base*100
    q30 = np.sum(in_array[:, 35:])/all_base*100
    q40 = np.sum(in_array[:, 45:])/all_base*100
    gc = np.sum(in_array[:, 2:4])/all_base*100
    E = error_rate(in_array[:, 5:])
    # print info
    print("#reads:", reads, "base:", all_base, file = f)
    print("#Q20:", q20, "#Q30:", q30, "#Q40:", q40, file =f)
    print("#GC content:", gc, "error:", E, file = f)
    # print header
    print("pos A T C G N ", file = f, end = '') # trim the \n
    for q in range(43):
        print(q, file = f, end = ' ')
    print(file = f)
    # print body
    reads_len = in_array.shape[0]
    for pos in range(reads_len):
        A = in_array[pos, 0]
        T = in_array[pos, 1]
        C = in_array[pos, 2]
        G = in_array[pos, 3]
        N = in_array[pos, 4]
        if A+T+C+G+N == 0:
            return 0
        ra, rt, rc, rg, rn = ration(A, T, C, G, N)
        print(pos, ra, rt, rc, rg, rn, file = f, sep = ' ', end = '')
        for q_plus in range(5, 48):
            qual_count = in_array[pos, q_plus]
            print(" ", qual_count, end = '', file = f)
        print(file = f)
    f.close()
    return(1)

if __name__ == "__main__":
    if sys.version_info < (3,0):
        print("python2 is not supported, please use python3")
        sys.exit(-1)
    if len(sys.argv) != 2:
        print("\n\tUsage:\tpython3 fqcheck.py <reads.fastq.gz>")
        exit(-1)
    infile = os.path.abspath(sys.argv[1])
    file_name = os.path.split(infile)[1]
    out_name = file_name.split('.')[0] + ".fqcheck"
    # the logest read length at present
    resultArray, reads = fqcheck(infile)
    printArray(resultArray, reads, out_name)
