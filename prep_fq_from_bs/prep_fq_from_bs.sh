#!/bin/bash

if [[ "$#" -ne 3 ]]; then
	echo "Usage: $0 <lane> <input_path> <output_path>"
	echo "<lane>: The lane number, either 1 or 2."
	echo "<input_path>: Directory path containing subdirectories with .fastq.gz files from a BaseSpace project."
	echo "<output_path>: Path to an empty directory for storing processed files. Existing data may be overwritten."
	echo 
	echo "Example of usage:"
	echo "./prep_fq_from_bs.sh 2 test_input_dir/two_lanes test_out_dir"
	
	exit 1
fi

# Assigning arguments to variables
lane="$1"
input_path="$2"
output_path="$3"

raw_reads_dir="${output_path}/raw_reads"
concat_dir="${output_path}/concat_raw_reads"

# Check if .fastq.gz files exist in the input directory
if ! find "$input_path" -mindepth 2 -maxdepth 2 -type f -name '*.fastq.gz' | read -r _; then
    echo "Error: No .fastq.gz files found in the input directory."
    exit 1
fi

# creating output dir for raw reads
mkdir -p "$raw_reads_dir" || { echo "Failed to create 'raw_reads' directory in specified output path."; exit 1; }
echo "Directory 'raw_reads' created successfully in ${output_path}."
echo

# reusable function to copy files to raw reads directory with new suffix
copy_and_rename() {
	local curr_suffix="$1"
	local repl_suffix="$2"

	for curr_file_path in "${input_path}"/*/*"${curr_suffix}".fastq.gz; do
		local file_basename=$(basename "$curr_file_path" "$curr_suffix".fastq.gz)
		local new_file_path="$raw_reads_dir/${file_basename}${repl_suffix}.fastq.gz"
		if cp "$curr_file_path" "$new_file_path"; then
			echo "File copied: $new_file_path"
		else
			echo "Error copying $curr_file_path to $new_file_path"
			exit 1
		fi
	done
}

# if working with one lane files: copy to raw_reads with cleaned file name
if [[ $lane == 1 ]]; then
	copy_and_rename "_L001_R1_001" "_R1"
	copy_and_rename "_L001_R2_001" "_R2"

	echo 
	echo "All done!"
# else if working with two lane files: copy to raw_reads with cleaned file name and then concatinate and store
# concatinated files with clean name in a separate directory - concat_dir.
elif [[ $lane == 2 ]]; then
	copy_and_rename	"_R1_001" "_R1"
	copy_and_rename "_R2_001" "_R2"
	echo

	mkdir -p "$concat_dir" || { echo "Failed to create 'concat_dir' directory in specified output path."; exit 1; }
	echo "Directory 'concat_dir' created successfully in ${output_path}."
	echo

	for read_type in R1 R2; do
		for file_path in "${raw_reads_dir}"/*_L001_${read_type}.fastq.gz; do
			file_basename=$(basename "$file_path" "_L001_${read_type}.fastq.gz")
			second_lane_file="${raw_reads_dir}/${file_basename}_L002_${read_type}.fastq.gz"
			concat_file_path="${concat_dir}/${file_basename}_${read_type}.fastq.gz"
			
			if [[ -f "$second_lane_file" ]]; then
				echo "Concatenating file: $concat_file_path ..."
				cat "$file_path" "$second_lane_file" > "$concat_file_path" && echo "...done"
				echo
			else
				echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    				echo
    				echo "Error: Missing matching file for concatenation: $second_lane_file. No concatenation done for $file_basename"
				echo
				echo '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
				echo
			fi
		done
	done

	echo
	echo "All done!"
else
    echo "The second argument identifies numbers of lanes used, must be of value either 1 or 2"
fi
