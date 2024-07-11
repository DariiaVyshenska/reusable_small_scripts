This program requires Python >= 3.11.8

Installing requirements:
    pip install -r requirements.txt

Getting docs:
    python main.py -h

Running with test files:
    python main.py ./test_data/test_full.fastq.vcf ./test_data/test.gb ./test_out

Running tests.
    Single test file:
        python -m unittest tests.test_utils

    All test files:
        python -m unittest discover -s tests -v
