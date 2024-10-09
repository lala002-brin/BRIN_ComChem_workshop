
# Battery DFTB parameter set

This set of parameters is aiming for battery applications.

The electronic parameters for C, N, O, H, P, S, F, Cl, Br, I, Na, K, Mg, Ca, Zn
were regenerated based on the documentation of the 3ob-3-1 set[1-4].
The elec. part of Cu are regenerated from 3ob-Cu[5].
However, the some of the electronic parameters are not consistent in the 3ob set.
The following table shows those pairs having large deviation (in Hartree) of the SK integral between the renegerated and original electronic parts.
One should be aware that the result may be different from the original parameters when these pairs are involved in the system.

Pair    |   Max.E.   |    RMSE
--------|------------|-------------
F-C     | 0.00209297 | 0.00040939
F-N     | 0.00271642 | 0.00051622
F-O     | 0.00332760 | 0.00063843
F-H     | 0.00159405 | 0.00054389
F-P     | 0.00233669 | 0.00055443
F-S     | 0.00207399 | 0.00044730
F-Ca    | 0.01083224 | 0.00235338
F-Zn    | 0.03141643 | 0.00369441  
Br-H    | 0.54246325 | 0.05957708
I-Ca    | 0.43732046 | 0.09510129
K-F     | 0.09962418 | 0.01716912
K-Cl    | 0.10937548 | 0.02094110
K-Br    | 0.62274550 | 0.07870410
K-I     | 0.10837735 | 0.02186691
K-Na    | 0.04304122 | 0.00844340
K-Mg    | 0.04452440 | 0.00827042
K-Ca    | 0.09857526 | 0.02593881
K-Zn    | 0.06009092 | 0.00703405

It should be safe to use the F-X(except K-F) repulsive potentials directly.
The long range part of K-F was modified using the IBI for reproducing 
the K-F RDF in KIB system. However, the short-range repulsive potential might also 
need to benchmarked.


## Modifications

* S-N, S-F repulsive potentials for FSA, TFSA[6]
* B-{B, C, N, O, H, F}
* Li-{C, N, O, H, F, S}
* K-F long range repulsive potential from IBI for K-F RDF

## Missing Repulsive Potentials

* Li-{P, Cl, Br, I, Na, K, Mg, Ca, Zn, Cu, B}
* B-{P, S, Cl, Br, I, Na, K, Mg, Ca, Zn, Cu, Li}
* Cu-{F, Cl, Br, I, Na, K, Mg, Ca, Zn, B, Li}

## Known issues

* K-F dimer energy is shifted up due to the difference of the electronic part.
* Li-N repulsive potential overestimate the desolvation energy ~20 kcal/mol

## Max L

    *from 3ob set
    C   p
    N   p
    O   p
    H   s
    P   d
    S   d
    F   p
    Cl  d
    Br  d
    I   d
    Na  p
    K   p
    Mg  p
    Ca  p
    Zn  d
    Cu  d
    -------
    B   p
    Li  p

## Hubbard derivatives (Original)

    C  -0.1492
    N  -0.1535
    O  -0.1575
    H  -0.1857
    P  -0.14
    S  -0.11
    F  -0.1623
    Cl -0.0697
    Br -0.0573
    I  -0.0433
    Na -0.0454
    K  -0.0339
    Mg -0.02
    Ca -0.0340
    Zn -0.03
    ------------------------
    Cu  -0.20 -0.0575 -0.0575 (d,p,s-Hubbard deriv)
    Li -0.04548
    B  -0.145

## Spin constants

    C  = { -0.03061752 -0.02504967 -0.02504967 -0.02265064 }
    O  = { -0.03523542 -0.02956329 -0.02956329 -0.02784892 }
    N  = { -0.03317737 -0.02754849 -0.02754849 -0.02545455 }
    H  = { -0.07173556 }
    P  = { -0.02062161 -0.01608727  0.00020336
           -0.01608727 -0.01490371 -0.00017681
            0.00020336 -0.00017681 -2.63949078
        }
    S  = { -0.02137207 -0.01699037  0.00026614
           -0.01699037 -0.01548853  0.00004212
            0.00026614  0.00004212 -3.18002259
        }
    F  = { -0.03697653 -0.03124096 -0.03124096 -0.02990088 }
    Cl = { -0.02191485 -0.01774113  0.00009189
           -0.01774113 -0.01605922 -0.00002294
            0.00009189 -0.00002294 -3.48878109
        }
    Br = { -0.01849037 -0.01439000  0.00008071
           -0.01439000 -0.01378058 -0.00014687
            0.00008071 -0.00014687 -1.58834980
        }
    I  = { -0.01444503 -0.01129546  0.00005365
           -0.01129546 -0.01143000 -0.00041630
            0.00005365 -0.00041630 -1.56920339
        }
    Na = { -0.01523980 -0.01346111 -0.01346111 -0.02499995 }
    K  = { -0.01079106 -0.01079534 -0.01079534 -0.01785043 }
    Mg = { -0.01666926 -0.01303186 -0.01303186 -0.01895754 }
    Ca = { -0.01195978 -0.01052233 -0.01052233 -0.01364146 }
    Zn = { -0.01680321 -0.01206216 -0.00276473
           -0.01206216 -0.02801945 -0.00070146
           -0.00276473 -0.00070146 -0.01925079
        }
    Cu = { -0.01702603 -0.01239199 -0.00280982
           -0.01239199 -0.05275212 -0.00052475
           -0.00280982 -0.00052475 -0.01737911
        }
    B  = { -0.02727469 -0.02171668 -0.02171668 -0.01915499 }
    Li = { -0.01966335 -0.01680606 -0.01680606 -0.01845789 }

## References

[1] J. Chem. Theory Comput., 2013, 9, 338-354.
[2] J. Chem. Theory Comput., 2014, 10, 1518-1537.
[3] J. Phys. Chem. B, 2015, 119, 1062–1082.
[4] J. Chem. Theory Comput., 2015, 11, 332–342.
[5] J. Chem. Theory Comput., 2015, 11, 4205.

# Appendix

## Hubbard Parameter (Re-Calculated)

    Sym U(d, p, s)
    H   0.41962078
    C   0.36466849 0.39904098
    O   0.49540525 0.52859485
    N   0.43088883 0.46411371
    P   0.28935528 0.28935528 0.33758121
    S   0.32877539 0.32877539 0.37485085
    F   0.55864068 0.59247710
    Cl  0.36688325 0.36688325 0.41145798
    Br  0.32765529 0.32765529 0.38037390
    I   0.28422165 0.28422165 0.33125639
    Na  0.09056012 0.16525644
    K   0.08172847 0.13627568
    Mg  0.15096373 0.22463958
    Ca  0.12862766 0.17652060
    Zn  0.52318882 0.15666630 0.26628325
    Cu  0.42472067 0.09943886 0.23826767
    B   0.29613013 0.33375716
    Li  0.13184802 0.17406794

## Hubbard derivatives (Re-Calculated)

    Sym U(d, p, s)
    C  -0.17924029 -0.15883485
    N  -0.15423157 -0.13474689
    O  -0.15815389 -0.13933670
    H  -0.18784076
    P  -0.04312579 -0.04312579 -0.04469802
    S  -0.04818050 -0.04818050 -0.05314414
    F  -0.16162133 -0.14318475
    Cl -0.02275617 -0.02275617 -0.02266281
    Br -0.05097151 -0.05097151 -0.02918050
    I  -0.07883505 -0.07883505 -0.07957930
    Na  0.09474448 -0.05074831
    K   0.02325357 -0.03575783
    Mg -0.00920686 -0.06928404
    Ca -0.02565264 -0.04222498
    Zn -0.15424958 -0.22249612 -0.05700515
    ------------------------
    Cu -0.20124542 -0.35477667 -0.07747923
    Li -0.61347016 -0.06994545
    B  -0.14500867 -0.12068370




