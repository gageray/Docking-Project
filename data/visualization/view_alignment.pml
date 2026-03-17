cmd.read_pdbstr('''
HETATM    1  C01 FYP Z 402       1.708   5.441   0.261  1.00 12.47           C  
HETATM    2  C02 FYP Z 402       0.204   5.209   0.383  1.00 12.47           C  
HETATM    3  O03 FYP Z 402      -0.087   3.802   0.205  1.00 12.47           O  
HETATM    4  C04 FYP Z 402      -0.312   3.074   1.314  1.00 12.47           C  
HETATM    5  O05 FYP Z 402      -0.172   3.579   2.407  1.00 12.47           O  
HETATM    6  C06 FYP Z 402      -0.730   1.670   1.205  1.00 12.47           C  
HETATM    7  C07 FYP Z 402      -0.055   0.657   0.602  1.00 12.47           C  
HETATM    8  C08 FYP Z 402       1.290   0.758  -0.053  1.00 12.47           C  
HETATM    9  N09 FYP Z 402       1.247   0.236  -1.418  1.00 12.47           N  
HETATM   10  C10 FYP Z 402       0.759  -0.961  -1.748  1.00 12.47           C  
HETATM   11  O11 FYP Z 402       0.800  -1.254  -2.924  1.00 12.47           O  
HETATM   12  C12 FYP Z 402       0.183  -1.990  -0.885  1.00 12.47           C  
HETATM   13  C13 FYP Z 402      -0.530  -1.737   0.280  1.00 12.47           C  
HETATM   14  N14 FYP Z 402      -0.806  -0.456   0.738  1.00 12.47           N  
HETATM   15  C15 FYP Z 402      -1.933  -0.110   1.401  1.00 12.47           C  
HETATM   16  N16 FYP Z 402      -1.888   1.153   1.675  1.00 12.47           N  
HETATM   17  C17 FYP Z 402      -1.037  -2.785   1.032  1.00 12.47           C  
HETATM   18  C18 FYP Z 402      -0.823  -4.085   0.643  1.00 12.47           C  
HETATM   19  C19 FYP Z 402      -0.105  -4.350  -0.507  1.00 12.47           C  
HETATM   20  C20 FYP Z 402       0.396  -3.321  -1.267  1.00 12.47           C  
HETATM   21  C22 FYP Z 402       1.789   1.101  -2.460  1.00 12.47           C  
HETATM   22  F21 FYP Z 402       0.101  -5.630  -0.882  1.00 12.47           F  
CONECT    1    2
CONECT    2    3
CONECT    3    4
CONECT    4    5    5    6
CONECT    6    7    7   16
CONECT    7    8   14
CONECT    8    9
CONECT    9   10   21
CONECT   10   11   11   12
CONECT   12   13   13   20
CONECT   13   14   17
CONECT   14   15
CONECT   15   16   16
CONECT   17   18   18
CONECT   18   19
CONECT   19   20   20   22
END
''', 'crystal_pose')

cmd.read_pdbstr('''
ATOM      1  C   UNL     1      -0.426   3.113   1.259  1.00  0.00           C  
ATOM      2  O   UNL     1      -0.541   3.874   2.197  1.00  0.00           O  
ATOM      3  O   UNL     1      -0.227   3.523   0.002  1.00  0.00           O  
ATOM      4  C   UNL     1      -0.128   4.939  -0.231  1.00  0.00           C  
ATOM      5  C   UNL     1       1.300   5.430  -0.046  1.00  0.00           C  
ATOM      6  C   UNL     1      -0.487   1.654   1.353  1.00  0.00           C  
ATOM      7  N   UNL     1      -0.584   1.020   2.565  1.00  0.00           N  
ATOM      8  C   UNL     1      -0.455   0.709   0.344  1.00  0.00           C  
ATOM      9  C   UNL     1      -0.610  -0.250   2.306  1.00  0.00           C  
ATOM     10  N   UNL     1      -0.527  -0.503   0.968  1.00  0.00           N  
ATOM     11  C   UNL     1      -0.327   0.747  -1.137  1.00  0.00           C  
ATOM     12  C   UNL     1      -0.520  -1.745   0.316  1.00  0.00           C  
ATOM     13  N   UNL     1       0.995   0.259  -1.498  1.00  0.00           N  
ATOM     14  C   UNL     1      -1.350  -2.756   0.794  1.00  0.00           C  
ATOM     15  C   UNL     1       0.309  -1.972  -0.791  1.00  0.00           C  
ATOM     16  C   UNL     1       1.351  -1.017  -1.297  1.00  0.00           C  
ATOM     17  C   UNL     1      -1.379  -3.992   0.180  1.00  0.00           C  
ATOM     18  C   UNL     1       0.288  -3.225  -1.395  1.00  0.00           C  
ATOM     19  O   UNL     1       2.481  -1.441  -1.542  1.00  0.00           O  
ATOM     20  C   UNL     1      -0.558  -4.202  -0.914  1.00  0.00           C  
ATOM     21  C   UNL     1       1.971   1.247  -1.900  1.00  0.00           C  
ATOM     22  F   UNL     1      -0.579  -5.411  -1.528  1.00  0.00           F  
CONECT    1    2    2    3    6
CONECT    3    4
CONECT    4    5
CONECT    6    7    8    8
CONECT    7    9    9
CONECT    8   10   11
CONECT    9   10
CONECT   10   12
CONECT   11   13
CONECT   12   14   14   15
CONECT   13   16   21
CONECT   14   17
CONECT   15   16   18   18
CONECT   16   19   19
CONECT   17   20   20
CONECT   18   20
CONECT   20   22
END
''', 'qm_conformer')

hide everything
show sticks, crystal_pose
show sticks, qm_conformer
color cyan, crystal_pose
color magenta, qm_conformer
zoom
