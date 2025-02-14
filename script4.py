import zipfile
import time
import resource
from typing import List, Tuple, Generator
import io
import multiprocessing as mp
from multiprocessing import Manager
import os

class AluYFinder:
    def __init__(self, pattern: str):
        # Convert pattern to uppercase for consistent comparison
        self.pattern = pattern.upper()
        self.lps = self._compute_lps()
    
    def _compute_lps(self) -> List[int]:
        length = 0
        lps = [0] * len(self.pattern)
        i = 1
        while i < len(self.pattern):
            if self.pattern[i] == self.pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    def exact_match(self, text: str) -> List[Tuple[int, int]]:
        matches = []
        i = 0
        j = 0
        text = text.upper()
        
        while i < len(text):
            if self.pattern[j] == text[i]:
                i += 1
                j += 1
            
            if j == len(self.pattern):
                matches.append((i - j, i))
                j = self.lps[j - 1]
            elif i < len(text) and self.pattern[j] != text[i]:
                if j != 0:
                    j = self.lps[j - 1]
                else:
                    i += 1
        
        return matches

    def inexact_match(self, text: str, max_mismatches: int = 1) -> List[Tuple[int, int, int]]:
        matches = []
        text = text.upper()
        n = len(text)
        m = len(self.pattern)
        
        for i in range(n - m + 1):
            # Check for substitutions
            mismatches = 0
            match_found = True
            
            for j in range(m):
                if text[i + j] != self.pattern[j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        match_found = False
                        break
            
            if match_found:
                matches.append((i, i + m, mismatches))
            
            # Check for insertions
            if i < n - m:
                mismatches = 0
                match_found = True
                for j in range(m):
                    if j < m - 1 and text[i + j] != self.pattern[j]:
                        mismatches += 1
                        if mismatches > max_mismatches:
                            match_found = False
                            break
                if match_found:
                    matches.append((i, i + m - 1, mismatches))
            
            # Check for deletions
            if i < n - m - 1:
                mismatches = 0
                match_found = True
                for j in range(m):
                    if j < m and text[i + j + 1] != self.pattern[j]:
                        mismatches += 1
                        if mismatches > max_mismatches:
                            match_found = False
                            break
                if match_found:
                    matches.append((i, i + m + 1, mismatches))
        
        return matches

def process_chunk(args):
    """Process a single chunk of sequence data."""
    chr_name, sequence, chunk_start, pattern = args
    finder = AluYFinder(pattern)
    
    exact_matches = [
        (chr_name, start + chunk_start, end + chunk_start) 
        for start, end in finder.exact_match(sequence)
    ]
    
    inexact_matches = [
        (chr_name, start + chunk_start, end + chunk_start, mismatches) 
        for start, end, mismatches in finder.inexact_match(sequence)
    ]
    
    return exact_matches, inexact_matches

def read_fasta_chunks(filename: str, chunk_size: int = 1024*1024) -> Generator[Tuple[str, str, int], None, None]:
    """Read FASTA file in chunks and yield (chromosome, sequence_chunk, chunk_start)."""
    if filename.endswith('.zip'):
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            if filename == 'hg38.zip':
                fasta_name = 'ncbi_dataset/data/GCA_000001405.29/GCA_000001405.29_GRCh38.p14_genomic.fna'
            elif filename == 'T2T.zip':
                fasta_name = 'ncbi_dataset/data/GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna'
            
            with zip_ref.open(fasta_name) as f:
                text_f = io.TextIOWrapper(f, encoding='utf-8')
                yield from process_fasta_file(text_f, chunk_size)

def process_fasta_file(file_handle, chunk_size: int) -> Generator[Tuple[str, str, int], None, None]:
    current_chr = ""
    current_seq = []
    chunk_start = 0
    
    for line in file_handle:
        if line.startswith('>'):
            if current_chr and current_seq:
                sequence = ''.join(current_seq)
                for i in range(0, len(sequence), chunk_size):
                    chunk = sequence[i:i + chunk_size]
                    yield current_chr, chunk, i
            current_chr = line.strip()[1:].split()[0]
            current_seq = []
            chunk_start = 0
        else:
            current_seq.append(line.strip())
    
    if current_chr and current_seq:
        sequence = ''.join(current_seq)
        for i in range(0, len(sequence), chunk_size):
            chunk = sequence[i:i + chunk_size]
            yield current_chr, chunk, i

def write_matches_to_file(matches: List[Tuple], filename: str, match_type: str):
    """Write matches to a tab-delimited file."""
    os.makedirs('results', exist_ok=True)
    output_file = os.path.join('results', filename)
    
    with open(output_file, 'w') as f:
        # Write header
        if match_type == 'exact':
            f.write("Chromosome\tStart\tEnd\n")
        else:  # inexact
            f.write("Chromosome\tStart\tEnd\tMismatches\n")
            
        # Write matches
        for match in matches:
            if match_type == 'exact':
                chr_name, start, end = match
                f.write(f"{chr_name}\t{start}\t{end}\n")
            else:  # inexact
                chr_name, start, end, mismatches = match
                f.write(f"{chr_name}\t{start}\t{end}\t{mismatches}\n")

def process_assembly(filename: str, aluy_sequence: str) -> Tuple[List, List, float]:
    """Process a single assembly file using multiple processes."""
    num_processes = max(1, os.cpu_count() - 1)
    start_time = time.time()
    
    chunks = [(chr_name, sequence, chunk_start, aluy_sequence) 
             for chr_name, sequence, chunk_start in read_fasta_chunks(filename)]
    
    exact_matches = []
    inexact_matches = []
    
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_chunk, chunks)
        
        for exact, inexact in results:
            exact_matches.extend(exact)
            inexact_matches.extend(inexact)
    
    # Write results to files
    assembly_name = filename.split('.')[0]
    write_matches_to_file(exact_matches, f'{assembly_name}_exact_matches.tsv', 'exact')
    write_matches_to_file(inexact_matches, f'{assembly_name}_inexact_matches.tsv', 'inexact')
    
    processing_time = time.time() - start_time
    return exact_matches, inexact_matches, processing_time

def main():
    # AluY consensus sequence
    aluy_sequence = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    
    # Process each assembly
    assemblies = ['hg38.zip', 'T2T.zip']
    
    print("Starting AluY finder...")
    print(f"Results will be saved in the 'results' directory")
    
    for assembly in assemblies:
        print(f"\nProcessing {assembly}:")
        try:
            exact_matches, inexact_matches, processing_time = process_assembly(
                assembly, aluy_sequence
            )
            
            # Report results
            print(f"Exact matches found: {len(exact_matches)}")
            print(f"Inexact matches found: {len(inexact_matches)}")
            print(f"Results saved to:")
            print(f"  - results/{assembly.split('.')[0]}_exact_matches.tsv")
            print(f"  - results/{assembly.split('.')[0]}_inexact_matches.tsv")
            print(f"Processing time: {processing_time:.2f} seconds")
            
            # Report memory usage
            memory_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
            print(f"Peak memory usage: {memory_usage:.2f} KB")
            
        except Exception as e:
            print(f"Error processing {assembly}: {str(e)}")

if __name__ == "__main__":
    main()
