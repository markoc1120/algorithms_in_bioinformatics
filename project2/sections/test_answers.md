# Sequence Alignment Analysiss
## Question 1

Linear gap cost alignment (g(k)=5*k):

Optimal score: 226.0

Optimal alignment:

```
Aligned sequence 1 (105 bp):
tatgga-gagaataaaagaactgagagatct-aatgtcgcagtcccgcac
-tcgcgagatact-cactaagac-cactgtggaccatatggccataatca
aaaag
```
```
Aligned sequence 2 (105 bp):
-atggatgtcaatccga-ctctacttttcctaaaaattccagcgcaaaat
gccataag-caccacattcccttatactggagatcct-cca-tacagcca
tggaa
```
## Question 2

Affine gap cost alignment (g(k)=5+5k):

Optimal score: 266.0

Optimal alignment:

```
Aligned sequence 1 (103 bp):
tatggagagaataaaagaactgagagatct-aatgtcgcagtcccgcac-
tcgcgagatactcactaagac-cactgtggaccatatggccataatcaaa
aag
```
```
Aligned sequence 2 (103 bp):
-atggatgtcaatccgactctacttttcctaaaaattccagcgcaaaatg
ccataagcaccacattcccttatactggagatcctcca--tacagccatg
gaa
```
## Question 3

Score matrix for linear gap cost (g(k)=5k):

```
      seq1  seq2  seq3  seq4  seq5
seq1     0   226   206   202   209
seq2   226     0   239   223   220
seq3   206   239     0   219   205
seq4   202   223   219     0   210
seq5   209   220   205   210     0
```

## Question 4

Score matrix for affine gap cost (g(k)=5+5k):

```
      seq1  seq2  seq3  seq4  seq5
seq1     0   266   242   243   256
seq2   266     0   283   259   254
seq3   242   283     0   269   243
seq4   243   259   269     0   247
seq5   256   254   243   247     0
```

