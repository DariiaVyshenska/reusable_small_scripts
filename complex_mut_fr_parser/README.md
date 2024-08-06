# Notes on usage (draft)

## IMPORTANT

Input bam file MUST be indexed!

## running tests

Note: test files were too large to push to GitHub. Please, contact the owner if  you want to use them.
    Single test file:
        python -m unittest tests.test_utils

    All test files:
        python -m unittest discover -s tests -v

## running on a test file from command line

./main.py ./test_data/coreF44.fastq.gz.bam 6747 .

Expected full output:
    CODON,FREQUENCY
    GCT,0.4972
    ACT,0.3437
    GTT,0.1075
    CTG,0.0283
    TGA,0.0085
    GCG,0.0021
    GNT,0.0014
    TTG,0.0014
    ACC,0.0014
    GTC,0.0007
    CTC,0.0007
    ACN,0.0007
    ANT,0.0007
    CCT,0.0007
    GCN,0.0007
    GAA,0.0007
    GCC,0.0007
    ACG,0.0007
    GAT,0.0007
    GCA,0.0007
    TNA,0.0007

