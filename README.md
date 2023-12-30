# edit_distance
Edit distance, word error rate, hypothesis-reference alignment


## This library provides the following functions:

- edit_distance(r, h)
- calc_wer(r, h)
- align_hyp_to_ref(h, r)
- interval_hyp_to_ref(a, h_ivl)
- align_intervals_hyp_to_ref(h, r, h_ivls)

  Shinsuke Sakai [ sakai Mon Dec 25 23:51:15 2023 ]
                 [ sakai Sat Dec 30 23:58:00 2023 ]
---
```
align_intervals_hyp_to_ref(h, r, h_ivls):
    
PARAMETERS
    h             hypoethesis, a list of strings.
    r             reference, a list of strings.
    h_ivls        a list of intervals in h represented as 2-tuples.

RETURNS
    r_ivls        a list of intervals in r represented as 2-tuples.
    hp            hypoethesis with parens inserted.
    rp            reference with parens inserted.
```
 
