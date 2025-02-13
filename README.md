# Bioinformatics_249
This repository is a collection of scripts and submissions for the course CS 249 Algorithms in Bioinformatics.

## Project 0: Warmup

The objective was to implement a basic pattern matching algorithm to identify Alu sequences in different human genome assemblies.

### Implementation Description

The ```AluYFinder``` class is the core component that handles pattern matching. It contains the KMP (Knuth-Morris-Pratt) pattern matching algorithm implementation through its exact_match method, which finds exact occurrences of the AluY sequence in the genome. The class also implements inexact matching through the inexact_match method, which uses Levenshtein distance to find sequences that differ from the AluY pattern by at most one mutation (insertion, deletion, or substitution). However, the approach to Levenshtein distance is not true to the dynamic programming paradigm. This was done to enable local execution of the script.

The ```read_fasta_chunks``` function handles file I/O operations, specifically reading from zipped FASTA files containing genome sequences. It's designed to work with NCBI's directory structure and handles both GRCh38 (hg38) and T2T CHM13v2.0 assemblies. The function reads the file in chunks to manage memory efficiently. 

The ```process_chunk``` function is the workhorse that processes individual chunks of sequence data. It takes a chunk of genome sequence and finds both exact and inexact matches of the AluY pattern within that chunk. This function is designed to work independently, making it suitable for parallel processing. 

The ```process_assembly``` function orchestrates the parallel processing of the genome assembly [^1]. It creates a pool of worker processes, distributes chunks of the genome to these processes, and collects the results. This function manages the multiprocessing aspect of the program, ensuring efficient CPU utilization while keeping memory usage under control. 

Finally, the ```main``` function ties everything together. It initializes the AluY sequence, processes each assembly file, and reports the results including the number of matches found, processing time, and memory usage. It also handles any errors that might occur during processing.




[^1]: Claude ver3.5 Sonnet was used to assist in making the algorithm capable of running across multiple cores. (Prompt_Given: Suggest methods to implement multi-core processing for achieving the objective keeping in mind the limitations of main memory.) Based on its suggestions, the ```process_assembly``` function was designed to keep CPU usage in the range 80 - 93% with one core left for system processes.



