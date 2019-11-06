#!/usr/bin/env python
# coding: utf-8


from scipy.stats import binom, hypergeom
from tqdm import tqdm
import itertools

import matplotlib as mpl
mpl.use('Agg')

# take into account error correction that can correct up to te errors
def Binom_prob(prob, ne, te):
    res = binom.sf(te, ne, p=prob)
    return res

# calculate the theoretical predictions
def checktheory(thres, n, ne, p, te, s):
    # calculate failure probability for a certain number of ones in se
    pfaildict = {}

    for se in range(0, n + 1):
        tmp_0 = [binom.sf((thres + i + se) / 2, se, p=0.5) * s[i] for i in s]
        tmp_1_1 = [binom.sf((thres + i + se + 1) / 2, se, p=0.5) * s[i] for i in s]
        tmp_1_2 = [binom.sf((thres + i + se - 1) / 2, se, p=0.5) * s[i] for i in s]
        pfail = sum(tmp_0) + (sum(tmp_1_1) + sum(tmp_1_2)) * 0.5
        pfaildict[se] = pfail

    # set everything to zero
    fail = 0

    # loop over all norm values
    for l1, l2 in tqdm(itertools.combinations_with_replacement(range(0, n + 1), 2), leave=False, total=n * (n + 1) / 2):
        # probability of a certain norm
        # pl = P[||s||2] * P[||c||2]
        pl1 = binom.pmf(l1, n=n, p=p)
        pl2 = binom.pmf(l2, n=n, p=p)
        pl = pl1 * pl2
        if l1 != l2:
            pl *= 2

        # skip if probability is too small
        if pl < 2**-100: ## 200
            continue

        # calculate the probability of a failure
        failtmp = 0
        # loop over all possible number of nonzero elements in se
        for se1 in range(max(0, l1 + l2 - n), min(l1, l2) + 1):
            # probability of number of nonzero elements in se
            pse = hypergeom.pmf(k=se1, M=n, n=l1, N=l2)
            # probability of failure for a certain se1
            pfail = pfaildict[se1]
            # weighted average share
            failtmp += pse * pfail*0.5 # failtmp = pb
        
        # for new model, take error correction into account
        fail += pl * Binom_prob(failtmp, ne, te) # Binom_prob : 1 - Binom(d, lm, pb)

    return fail


def main():
    # maximum error correction to plot
    te = 15
    # scheme to use
    q = 256
    p = 64
    n = 512 #1024 for strong
    LizarMong = {}
    LizarMong['thres'] = q / 4 - q / (2 * p)
    LizarMong['s'] = {-1: 1. / 8, 0: 3. / 4, 1: 1. / 8}   #hs=128 with Comfort
    #LizarMong['s'] = {-1: 1. / 16, 0: 7. / 8, 1: 1. / 16}   #hs=128 with Strong
    LizarMong['e'] = {-1: 1. / 4, 0: 1. / 2, 1: 1. / 4}
    LizarMong['sprime'] = LizarMong['s']
    LizarMong['eprime'] = 0
    LizarMong['eprimeprime'] = 0
    LizarMong['n'] = n
    LizarMong['n2'] = 1  # ??
    LizarMong['name'] = 'LizarMong'
    scheme = LizarMong

    # load some parameters
    n = 2 * scheme['n']
    s = scheme['s']
    thres = scheme['thres']
    alg = scheme['name']
    ne = n / 2 - 1
    p = s[1] + s[-1]
    errors = list(range(0, te + 1))

    test = checktheory(thres, n, ne, p, errors, s)
    print(test) 
    f = open("LizarMong_result.txt", 'w')
    f.write(str(test))
    f.close()



if __name__ == '__main__':
    main()






