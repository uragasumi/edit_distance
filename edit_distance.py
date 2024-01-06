#  
#  This library provides the following functions:
# 
#  (1) edit_distance(r, h)
#  (2) calc_wer(r, h)
#  (3) align_hyp_to_ref(h, r)
#  (4) interval_hyp_to_ref(a, h_ivl)
#  (5) align_intervals_hyp_to_ref(h, r, h_ivls)
#
#  Shinsuke Sakai [ sakai Mon Dec 25 23:51:15 2023 ]
#                 [ sakai Sat Dec 30 23:58:00 2023 ]

import sys
import numpy as np


def str_to_int(tok, tok2i):
    '''
RETURNS
    integer representation of the token 'tok'
    '''

    if tok not in tok2i:
        tok2i[tok] = len(tok2i)
    
    return tok2i[tok]


def str_list_to_int_array(l, tok2i):
    '''
RETURNS
    a      a numpy array of length len(l) + 1 representing string tokens in
           the input list l.  It has an azdditional zero-th cell at the left
           end of the array to make 1-indexing possible for 'a'.

    '''
    a = np.zeros(len(l) + 1, dtype=int)
    a[0] = -1

    for i, tok in enumerate(l, 1):   # 1-indexing
        a[i] = str_to_int(tok, tok2i)
    
    return a


def edit_distance(r, h):
    '''
SYNOPOSIS
    edit_distance(r, h) computes edit distance aka Levenshtein distance using
    a Wagner and Fischer's algorithm. (They seem to have variants.)

PARAMETERS
    r    reference token sequence, a python list of strings.
    h    hypothesized token sequence, a python list of strings.

RETURNS
    m    number of tokens in reference.
    b    backtrace matrix
    '''

    tok2i = {}
    m = len(r)
    n = len(h)
    rix = str_list_to_int_array(r, tok2i)  # len: m + 1
    hix = str_list_to_int_array(h, tok2i)  # len: n + 1
    
    d = np.zeros((m+1, n+1), dtype=int)
    b = np.ones((m+1, n+1), dtype=int) * -1   # [1]
    
    # source prefixes can be transformed into empty string by
    # dropping all characters
    d[:, 0] = range(m+1)
    b[1:, 0] = 0                              # [1]
    # target prefixes can be reached from empty source prefix
    # by inserting every character
    d[0, :] = range(n+1)
    b[0, 1:] = 1                              # [1]

    n_del = n_ins = n_subst = 0

    for j in range(1, n+1):
        for i in range(1, m+1):
            if rix[i] == hix[j]:
                subst_cost = 0
            else:
                subst_cost = 1

            del_ins_subst = (d[i-1, j] + 1,             # deletion
                             d[i, j-1] + 1,             # insertion
                             d[i-1, j-1] + subst_cost)  # subst or right
            

            min_arg = np.argmin(del_ins_subst)
            
            d[i, j] = del_ins_subst[min_arg]
            b[i, j] = min_arg
            # backtrace + subst flag
            if min_arg == 2:
                b[i, j] = min_arg + subst_cost
            else:
                b[i, j] = min_arg

    # end for
    e = d[m, n]
    
    return e, b

#    n_ins, n_del, n_subst = count_ops(b, m, n)
#    return m, e, n_ins, n_del, n_subst

def count_ops(b, m, n):
    '''
SYNOPSIS    
    count the number of insert, delete, and substitute operations along the
    optimal path from (m, n) back to (0, 0).

PARAMETERS
    b      a backtrace matrix of the size (m+1) x (n+1).
    m      #  of ref tokens
    n      #  of hyp tokens
RETURNS
    n_del
    n_ins
    n_subst

    '''

    n_del = n_ins = n_subst = 0
    i = m
    j = n
    while not (i == 0 and j == 0):
        # print("[i, j] =", [i, j])
        # print("b[i, j] =", b[i, j])

        if b[i, j] == 0:
            n_del = n_del + 1
            i = i - 1
        elif b[i, j] == 1:
            n_ins = n_ins + 1
            j = j - 1
        else:
            if b[i, j] == 3:
                n_subst = n_subst + 1
            i = i - 1
            j = j - 1
    # end while
    
    return n_ins, n_del, n_subst


def calc_wer(r, h):
    '''
SYNOPOSIS
    calc_wer(r, h) computes edit distance aka Levenshtein distance using
    a Wagner and Fischer's algorithm. (They seem to have variants.)

PARAMETERS
    r    reference token sequence, a python list of strings.
    h    hypothesized token sequence, a python list of strings.

RETURNS
    m    number of tokens in reference.
    e    edit distance, an integer.
    i    number of insertions
    d    number of deletions
    s    number of substitutions
    '''

    m = len(r)  # number of tokens in reference
    n = len(h)  # number of tokens in hypothesis

    # e: edit distance
    # b: backtrace matrix
    e, b = edit_distance(r, h)

    n_ins, n_del, n_subst = count_ops(b, m, n)

    return m, e, n_ins, n_del, n_subst


def make_alignment_table(b, m, n):
    '''
SYNOPSIS
    make_alignment_table(b, m, n) makes a table of lengh n + 1 that maps a
    position j in the hypothesis to the matching positions (i_1, i_2, ..)  in
    the reference based on the optimal path from (m, n) back to (0, 0) in b.

PARAMETERS
    b      a backtrace matrix of the size (m+1) x (n+1).
    m      #  of ref tokens
    n      #  of hyp tokens
RETURNS
    a      an alignment table of length n + 1. An entry a[j] is a 1- or
           2-tuple representing an interval in the reference. 
    '''
    # n_del = n_ins = n_subst = 0
    a = [()] * (n + 1)
    i = m
    j = n

    # - meanings of b[i, j] values:
    #   b[i, j] == 0  del
    #              1  ins
    #              2  right
    #              3  subt

    while not (i == 0 and j == 0):
        # print("[i, j] =", [i, j])
        # print("b[i, j] =", b[i, j])
        # print('a[', j, ']:', a[j])
        a[j] = (i,) + a[j]

        if b[i, j] == 0:
            # n_del = n_del + 1
            i = i - 1
        elif b[i, j] == 1:
            # n_ins = n_ins + 1
            j = j - 1
        else:
            # if b[i, j] == 3:
            #     n_subst = n_subst + 1
            i = i - 1
            j = j - 1
    # end while

    # print('n_ins, n_del, n_subst =', n_ins, n_del, n_subst)
    # return n_ins, n_del, n_subst
    return a


def make_ref_to_hyp_alignment_table(b, m, n):
    '''
SYNOPSIS
    make_alignment_table(b, m, n) makes a table of lengh n + 1 that maps a
    position j in the hypothesis to the matching positions (i_1, i_2, ..)  in
    the reference based on the optimal path from (m, n) back to (0, 0) in b.

PARAMETERS
    b      a backtrace matrix of the size (m+1) x (n+1).
    m      #  of ref tokens
    n      #  of hyp tokens
RETURNS
    a      an alignment table of length n + 1. An entry a[j] is a 1- or
           2-tuple representing an interval in the reference. 
    '''
    # n_del = n_ins = n_subst = 0
    a = [()] * (m + 1)
    i = m
    j = n

    # - meanings of b[i, j] values:
    #   b[i, j] == 0  del
    #              1  ins
    #              2  right
    #              3  subt

    while not (i == 0 and j == 0):
        # print("[i, j] =", [i, j])
        # print("b[i, j] =", b[i, j])
        # print('a[', j, ']:', a[j])
        a[i] = (j,) + a[i]
        # a[j] = (i,) + a[j]

        if b[i, j] == 0:
            # n_del = n_del + 1
            i = i - 1
        elif b[i, j] == 1:
            # n_ins = n_ins + 1
            j = j - 1
        else:
            # if b[i, j] == 3:
            #     n_subst = n_subst + 1
            i = i - 1
            j = j - 1
    # end while

    # print('n_ins, n_del, n_subst =', n_ins, n_del, n_subst)
    # return n_ins, n_del, n_subst
    return a


def align_hyp_to_ref(h, r):
    '''
SYNOPSIS
    align_hyp_to_ref(h, r) generates an alignment table that helps mapping 
    an interval in h to an interval in r.

PARAMETERS
    r    reference token sequence, a python list of strings.
    h    hypothesized token sequence, a python list of strings.
RETURNS
    a    alignment table.
'''

    m = len(r)  # number of tokens in reference
    n = len(h)  # number of tokens in hypothesis

    e, b = edit_distance(r, h)
    
    a = make_alignment_table(b, m, n)

    return a


def align_ref_to_hyp(r, h):
    '''
SYNOPSIS
    align_hyp_to_ref(h, r) generates an alignment table that helps mapping 
    an interval in h to an interval in r.

PARAMETERS
    r    reference token sequence, a python list of strings.
    h    hypothesized token sequence, a python list of strings.
RETURNS
    a    alignment table.
'''

    m = len(r)  # number of tokens in reference
    n = len(h)  # number of tokens in hypothesis

    e, b = edit_distance(r, h)
    
    a = make_ref_to_hyp_alignment_table(b, m, n)

    return a


def interval_hyp_to_ref(a, h_ivl):
    # def interval_hyp_to_ref(a, h, r, h_ivl):
    '''
SYNPOSIS
    interval_hyp_to_ref(a, j1, j2) gets interval [i1, i2] in the reference r
    corresponding to the interval [j1, j2] in hypothesis h.

PARAMETERS
    a         alignment table generated by align_hyp_to_ref(h, r).
    h_ivl = (j1, j2)  a pair of 0-based indices representing an interval in h


RETURNS
    (i1, i2)   2-tuple representing an interval in r.

EXAMPLE
    given h = [a, b, a, c, a],
          r = [a, b, b, a, c, c, a],  and
          a = align_hyp_to_ref(h, r)

    interval_hyp_to_ref(a, h, r, (0, 2)) returns (0, 3)
    '''
    
    j1, j2 = h_ivl
    # we access alignment table with indices modified to 1-based
    i1 = a[j1+1][0]   # first tuple element
    i2 = a[j2+1][-1]  # last tuple element

    # print('h:', h[j1:j2+1])
    # print('r:', r[i1-1:i2-1+1])

    # convert from 1-based to 0-based
    return i1-1, i2-1


def interval_ref_to_hyp(a, r_ivl):
    '''
SYNPOSIS
    interval_hyp_to_ref(a, j1, j2) gets interval [i1, i2] in the reference r
    corresponding to the interval [j1, j2] in hypothesis h.

PARAMETERS
    a         alignment table generated by align_hyp_to_ref(h, r).
    h_ivl = (j1, j2)  a pair of 0-based indices representing an interval in h


RETURNS
    (i1, i2)   2-tuple representing an interval in r.

EXAMPLE
    given r = [a, b, b, a, c, c, a],
          h = [a, b, a, c, a],        and
          a = align_ref_to_hyp(r, h)

    interval_ref_to_hyp(a, (0, 2)) returns (0, 1)
    '''
    
    i1, i2 = r_ivl
    # we access alignment table with indices modified to 1-based
    j1 = a[i1+1][0]   # first tuple element
    j2 = a[i2+1][-1]  # last tuple element

    # print('r:', r[i1:i2+1])
    # print('h:', h[j1-1:j2-1+1])

    # convert from 1-based to 0-based
    return j1-1, j2-1

    

def insert_parens(tokens, ivls, lpar, rpar):
    '''
PARAMETERS
    tokens       a list of strings.
    ivls         a list of 2-tuples representing intervals in tokens.
    lpar
    rpar

RETURNS
    tokens with parens inserted to highlight the designated intervals.

EXAMPLE
    tokens = ['a', 'a', 'b', 'b', 'a', 'a', 'c', 'c', 'a']
    ivls = [(2,3), (6,7)]
    lpar = '('
    rpar = ')'

    returns ['a', 'a', '(', 'b', 'b', ')', 'a', 'a', '(', 'c', 'c', ')', 'a']
    '''

    tokens_p = tokens.copy()

    gap = 0
    for ivl in ivls:
        left, right = ivl
        tokens_p.insert(left + gap, lpar)
        gap = gap + 1
        tokens_p.insert(right + 1 + gap, rpar)
        gap = gap + 1
    
    return tokens_p


def align_intervals_hyp_to_ref(h, r, h_ivls):
    '''
PARAMETERS
    h             hypoethesis, a list of strings.
    r             reference, a list of strings.
    h_ivls        a list of intervals in h represented as 2-tuples.

RETURNS
    r_ivls        a list of intervals in r represented as 2-tuples.
    hp            hypoethesis with parens inserted.
    rp            reference with parens inserted.
    '''

    # these may be customized by optional parameters
    lpar = '('
    rpar = ')'

    a = align_hyp_to_ref(h, r)
    r_ivls = []
    for h_ivl in h_ivls:
        r_ivl = interval_hyp_to_ref(a, h_ivl)
        r_ivls.append(r_ivl)
        
    hp = insert_parens(h, h_ivls, lpar, rpar)
    rp = insert_parens(r, r_ivls, lpar, rpar)

    return r_ivls, hp, rp


def align_intervals_ref_to_hyp(r, h, r_ivls):
    '''
PARAMETERS
    r             reference, a list of strings.
    h             hypoethesis, a list of strings.
    r_ivls        a list of intervals in r represented as 2-tuples.

RETURNS
    h_ivls        a list of intervals in h represented as 2-tuples.
    rp            reference with parens inserted.
    hp            hypoethesis with parens inserted.
    '''

    # these may be customized by optional parameters
    lpar = '('
    rpar = ')'

    a = align_ref_to_hyp(r, h)
    h_ivls = []

    for r_ivl in r_ivls:
        h_ivl = interval_ref_to_hyp(a, r_ivl)
        h_ivls.append(h_ivl)
        
    rp = insert_parens(r, r_ivls, lpar, rpar)
    hp = insert_parens(h, h_ivls, lpar, rpar)

    return h_ivls, rp, hp


