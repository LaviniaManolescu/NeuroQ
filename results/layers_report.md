# QAOA pe p = 1, 2, 3 (INTERP)

Seed 71, 2048 măsurători, COBYLA din 2 porniri (prima INTERP din p−1) pe statevector; codificarea mică a fiecărei instanțe. Adâncimea transpilată: rz, sx, x, cx, nivel 3, fără hartă de qubiți.

| instanța | p | optim (freq) | oglinda (freq) | optim+oglindă | cel mai probabil | freq | cost | adâncime | adâncime transpilată |
|---|---|---|---|---|---|---|---|---|---|
| notebook | 1 | 01101 20.70% | 10010 22.12% | 42.82% (2 stări) | 10010 | 22.12% | 1 | 10 | 30 |
| notebook | 2 | 01101 38.43% | 10010 41.46% | 79.88% (2 stări) | 10010 | 41.46% | 1 | 16 | 50 |
| notebook | 3 | 01101 41.60% | 10010 43.70% | 85.30% (2 stări) | 10010 | 43.70% | 1 | 22 | 70 |
| illumina | 1 | 0011001010 0.15% | 1100110101 0.39% | 0.54% (2 stări) | 0010110111 | 1.37% | 14 | 18 | 54 |
| illumina | 2 | 0011001010 9.18% | 1100110101 8.01% | 17.19% (2 stări) | 0011001010 | 9.18% | 0 | 24 | 74 |
| illumina | 3 | 0011001010 10.25% | 1100110101 9.13% | 19.38% (2 stări) | 0011001010 | 10.25% | 0 | 30 | 94 |
| nanopore | 1 | 010001001 0.98% | 101110110 1.17% | 2.15% (2 stări) | 101110010 | 4.30% | 25 | 15 | 45 |
| nanopore | 2 | 010001001 6.45% | 101110110 6.64% | 13.09% (2 stări) | 101110110 | 6.64% | 7 | 24 | 74 |
| nanopore | 3 | 010001001 8.84% | 101110110 10.94% | 19.78% (2 stări) | 101110110 | 10.94% | 7 | 33 | 103 |
