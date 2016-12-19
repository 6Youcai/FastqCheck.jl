#!/usr/bin/env python3

import os
import sys
import gzip
import numpy as np

# fqcheck
def fqcheck(infile, read_length):
    # open file
    if infile.endswith(".gz"):
        hand = gzip.open(infile)
    else:
        hand = open(infile)
    # init array
    m = (read_length, 48)
    out_array = np.zeros(m, int16)
    order_dict = {"A":1, "T":2, "C":3, "G":4, "N":5}
    reads = 0
    for line1 in hand:
        line2 = hand.readline()
        _, line4 = hand.readline(), hand.readline()
        reads += 1
        Length = len(line2)
        for i in range(Length):
            base = line2[i]
            base_in_array = order_dict[base]
            out_array[i, base_in_array] += 1
            qual = ord(line4[i]) - 33 # phred 33
            # the first 5 col is ATCGN, and qual start from 0
            qual_in_array = qual + 6
            out_array[i, qual_in_array] += 1
    hand.close()
    return(out_array)

def ration(a, t, c, g, n):
    Sum = a+t+c+g+n
    ra = a/Sum
    rt = t/Sum
    rc = c/Sum
    rg = g/Sum
    rn = n/Sum
    return(ra, rt, rc, rg, rn)

def error_rate(in_array)
    all_base = sum(in_array)
    sum_qual = 0
    for i in range(43):
        tmp = sum(in_array[:, i])
        sum_qual += tmp * i
    qual = sum_qual/all_base
    E = 10**(-1*qual/10)

def printArray(in_array, reads, outfile)
    f = open(outfile, "w")
    all_base = sum(in_array[:, 5:])
    q20 = sum(in_array[:, 25:])/all_base
    q30 = sum(in_array[:, 35:])/all_base
    q40 = sum(in_array[:, 45:])/all_base
    gc = sum(in_array[:, 2:3])/all_base
    E = error_rate(in_array[:, 5:])
    # print info
    print("#reads:", reads, "base:", all_base, file = f)
    print("#Q20:", q20, "#Q30:", q30, "#Q40:", q40, file =f)
    print("#GC content:", qc, "error%:", E, file = f)
    # print header
    print("pos A T C G N", file = f, end = '') # trim the \n
    for q in range(43):
        print(q, file = f, end = '') # trim \n
    print(file = f)
    # print body
    reads_len = in_array.shape[0]
    for pos in range(reads_len):
        A = in_array[pos, 0]
        T = in_array[pos, 1]
        C = in_array[pos, 2]
        G = in_array[pos, 3]
        N = in_array[pos, 4]
        if A+T+C+G+N == 0
            return 0
        ra, rt, rc, rg, rn = ration(A, T, C, G, N)
        print(pos, ra, rt, rc, rg, rn, file = f, sep = ' ', end = '')
        for q_plus in range(5, 49):
            qual_count = in_array[pos, q_plus]
            print(" ", qual_count, end = '', file = f)
        print(file = f)
    f.close()
    return(1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("\n\tUsage:\tpython3 <reads.fastq.gz>")
        exit(-1)
    infile = os.path.abspath(sys.argv[1])
    file_name = os.path.split(infile)[1]
    out_name = file_name.split('.')[0] + ".fqcheck"
    # the logest read length at present
    resultArray, reads = fqcheck(infile, reads_len = 151)
    printArray(resultArray, out_name)















    
