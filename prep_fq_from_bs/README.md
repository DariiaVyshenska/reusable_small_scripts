# Prep_fq_from_bs.sh

## Overview
This bash script facilitates consistently copying, concatenating, and renaming BaseSpace sequencing files for streamlined storage and analysis preparation. Ensure execution privileges are set for your user by running the command `chmod +x prep_fq_from_bs.sh`.

## Usage
`./prep_fq_from_bs.sh <lane> <input_path> <output_path>`

**Where**:
- `<lane>`: The lane number, either 1 or 2.
- `<input_path>`: Directory path containing subdirectories with .fastq.gz files from a BaseSpace project.
- `<output_path>`: Path to an empty directory for storing processed files. Existing data may be overwritten.

> [!IMPORTANT]
> Input files must adhere to the naming pattern: `<uniq-sample-id>_L00<lane-number>_R<read-number>_001.fastq.gz`. For example, R1 reads generated on lane #1 for sample id `my_awesom_sample` will have a file name `my_awesom_sample_L001_R1_001.fastq.gz`


Example of calling the script:
```bash
./prep_fq_from_bs.sh 2 test_input_dir/two_lanes test_out_dir
```
#### Input

Example of the input directory structure and file naming:

*one lane*
```
test_input/one_lane
├── 240215_beads_1_L1-ds.537d05de13b544d98bbb288568246a83
│   ├── 240215_beads_1_S5_L001_R1_001.fastq.gz
│   └── 240215_beads_1_S5_L001_R2_001.fastq.gz
├── 240215_beads_2_L1-ds.40d36f77fe0f4d4ebbfc8db175241ebf
│   ├── 240215_beads_2_S6_L001_R1_001.fastq.gz
│   └── 240215_beads_2_S6_L001_R2_001.fastq.gz
├── 240215_pos_1_L1-ds.c2c587d4529843b9a936c9507974b9a8
│   ├── 240215_pos_1_S1_L001_R1_001.fastq.gz
│   └── 240215_pos_1_S1_L001_R2_001.fastq.gz
└── 240215_pos_2_L1-ds.dc4740b1736f44bba59e07da0596f76b
    ├── 240215_pos_2_S2_L001_R1_001.fastq.gz
    └── 240215_pos_2_S2_L001_R2_001.fastq.gz
```
*two lanes*
```
test_input/two_lanes
├── 240315_A4_L001-ds.db61922576934404aeb6fef9f4c7d7b3
│   ├── 240315-A4_S31_L001_R1_001.fastq.gz
│   └── 240315-A4_S31_L001_R2_001.fastq.gz
└── 240315_A4_L002-ds.08b406f6ed7a4a78ac2994c387907e9b
    ├── 240315-A4_S31_L002_R1_001.fastq.gz
    └── 240315-A4_S31_L002_R2_001.fastq.gz
```
#### Output

> [!CAUTION]
> The specified output directory must be empty or not exist to prevent data loss from overwriting.

For **one-lane** input, the script will:
- in the specified by the user `<output_path>`, it will create a new directory named `raw_reads` and copy each file with the updated file name: `<uniq-sample-id>_R<read-number>.fastq.gz`.

For **two-lane** input, the script will:
- in the specified by the user `<output_path>`, it will create a new directory named `raw_reads` and copy each file with the updated file name: `<uniq-sample-id>_L00<lane-number>_R<read-number>.fastq.gz`.
- in the specified by the user `<output_path>`, it will create a new directory named `concat_raw_reads` and concatinate files from two lanes into one file for each read # with the updated name: `<uniq-sample-id>_R<read-number>.fastq.gz`


#### Best Practices
- Verify the input path for correct file naming and organization before execution.
- Ensure the output directory is empty or specify a new directory to avoid unintentional data loss.

By adhering to these guidelines and instructions, users can efficiently prepare sequencing data from BaseSpace for further analysis or storage.
