## Project 0: Warmup

The objective was to implement a basic pattern matching algorithm to identify Alu sequences in different human genome assemblies.

### Implementation Description ([Source_Code](script4.py))

The ```AluYFinder``` class is the core component that handles pattern matching. It contains the KMP (Knuth-Morris-Pratt) pattern matching algorithm implementation through its exact_match method, which finds exact occurrences of the AluY sequence in the genome. The class also implements inexact matching through the inexact_match method, which uses Levenshtein distance to find sequences that differ from the AluY pattern by at most one mutation (insertion, deletion, or substitution). However, the approach to Levenshtein distance is not true to the dynamic programming paradigm. This was done to enable local execution of the script. The ```read_fasta_chunks``` function handles file I/O operations, specifically reading from zipped FASTA files containing genome sequences. It's designed to work with NCBI's directory structure and handles both GRCh38 (hg38) and T2T CHM13v2.0 assemblies. The function reads the file in chunks to manage memory efficiently.  The ```process_chunk``` function is the workhorse that processes individual chunks of sequence data. It takes a chunk of genome sequence and finds both exact and inexact matches of the AluY pattern within that chunk. This function is designed to work independently, making it suitable for parallel processing. The ```write_matches_to_file``` function handles the output of match locations to tab-separated files. It creates separate files for exact and inexact matches, with each entry containing the chromosome name, start position, and end position. For inexact matches, it also includes the number of mismatches observed. The ```process_assembly``` function orchestrates the parallel processing of the genome assembly [^1]. It creates a pool of worker processes, distributes chunks of the genome to these processes, and collects the results. After processing, it writes the match locations to files in a 'results' directory, with separate files for exact and inexact matches from each assembly.  This function manages the multiprocessing aspect of the program, ensuring efficient CPU utilization while keeping memory usage under control. Finally, the ```main``` function ties everything together. It initializes the AluY sequence, processes each assembly file, and reports the results including the number of matches found, processing time, and memory usage. It also handles any errors that might occur during processing.

#### Requirements
- python 3.11 

### Results and Discussion

The script produced the following output

```
Starting AluY finder...
Results will be saved in the 'results' directory

Processing hg38.zip:
Exact matches found: 3
Inexact matches found: 50
Results saved to:
  - results/hg38_exact_matches.tsv
  - results/hg38_inexact_matches.tsv
Processing time: 1543.95 seconds
Peak memory usage: 1003168.00 KB

Processing T2T.zip:
Exact matches found: 2
Inexact matches found: 62
Results saved to:
  - results/T2T_exact_matches.tsv
  - results/T2T_inexact_matches.tsv
Processing time: 1386.78 seconds
Peak memory usage: 1456656.00 KB
```

After discussion with colleagues, it was ascertained that the exact_match method shows results that are expected and accurate. However, inexact matches are slightly inflated. This can be attributed to the approach for Levenshtein distance, which has resulted in overcounting as the mutations in current implementation are treated as independant. Furthermore, overlaps resulted in a sliding window effect where a single variant AluY sequence could be counted multiple times.

#### After processing the data, we found that the GRCh38 assembly had 13 inexact matches with a maximum of 1 mismatch, while the T2T CHM13v2.0 assembly had 16 inexact matches under the same condition.

The GRCh38 assembly took approximately 25 minutes to process with a peak memory usage of 1.003 GB, while the T2T CHM13v2.0 assembly processed faster at 23 minutes but required nearly double the memory at 1.456 GB. This higher memory requirement for T2T likely reflects its more complete coverage of complex genomic regions.
From a biological perspective, the T2T assembly yielded more total matches (64) compared to GRCh38 (53), suggesting better capture of Alu sequences in previously unresolved regions. Both assemblies showed a strong preference for inexact matches, with T2T having a slightly higher proportion (97% vs 94% in GRCh38). This aligns with the biological understanding that Alu elements accumulate mutations over time.

### Matches found in GRCh38 assembly ([Exact_Matches](hg38_exact_matches.tsv))([Inexact_Matches](hg38_inexact_matches.tsv))
### Matches found in T2T CHM13v2.0 assembly ([Exact_Matches](T2T_exact_matches.tsv))([Inexact_Matches](T2T_inexact_matches.tsv))

[^1]: Claude ver3.5 Sonnet was used to assist in making the algorithm capable of running across multiple cores. (Prompt_Given: Suggest methods to implement multi-core processing for achieving the objective keeping in mind the limitations of main memory.) Based on its suggestions, the ```process_assembly``` function was designed to keep CPU usage in the range 80 - 93% with one core left for system processes.
