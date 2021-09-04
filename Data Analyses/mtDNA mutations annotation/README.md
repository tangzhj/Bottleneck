# Calculate the number of all the possible synonymous and non-synonymous SNV in mitochondrial genome with PAML

## Step 1 Generate PAML input

We obtained CDS of single gene in mitochondrial genome and  generate human_same_condon.txt as input to PAML. 

## Step 2 PAML calculation

- PAML can calculate the calculate the number of all the possible synonymous (defined as S value) and non-synonymous (defined as N value) SNV sites on the CDS region. The command line is as follows:

```shell
codeml human_same.ctl
```

- We performed the same step to all 13 genes on mitochondrial genome, the output files can be viewed in Result directory.

## Step 3 Results integration

- We assumed S value of the complete mitochondrial genome was equal to the sum of S values of every single gene in mitochondria. We expected the same from N value.
- The result table can be viewed in Result directory.  

## anovar.sh
- We use anovar software to calculate variant type.