#!/usr/bin/env python3
"""
OLC Assembler Script
This script implements an Overlap-Layout-Consensus (OLC) assembly algorithm for assembling DNA sequences 
from FASTQ files. It uses MinHash sketches for efficient overlap detection and supports multi-threaded 
execution for performance optimization.
Usage:
    python olc_assembler_v2.py --fastq <input_fastq_files> --output <output_prefix> [options]
Arguments:
    --fastq            List of input FASTQ files to be assembled (required).
    --output           Prefix for output files (required).
    --min-overlap      Minimum overlap length between reads (default: 20).
    --min-identity     Minimum overlap identity between reads (default: 0.9).
    --kmer-size        K-mer size for MinHash sketching (default: 15).
    --sketch-size      Number of k-mers to retain in each MinHash sketch (default: 100).
    --threads          Number of CPU cores to use for parallel processing (default: 1).
                       Use 0 to utilize all available cores.
Outputs:
    - A FASTA file containing the assembled contigs (<output_prefix>.fasta).
    - A text file with assembly metrics (<output_prefix>_metrics.txt).
Assembly Metrics:
    - Total assembly length: Total length of all assembled contigs.
    - Number of contigs: Total number of assembled contigs.
    - GC content: Percentage of GC bases in the assembly.
    - Largest contig: Length of the largest contig.
    - N50: Contig length at which 50% of the assembly is contained in contigs of this length or longer.
    - N90: Contig length at which 90% of the assembly is contained in contigs of this length or longer.
    - L50: Number of contigs that make up 50% of the assembly.
Example:
    python olc_assembler_v2.py --fastq sample1.fastq sample2.fastq --output assembly_result \
        --min-overlap 30 --min-identity 0.95 --threads 4

"""



import os
import argparse
import multiprocessing
from functools import partial
from collections import defaultdict, deque, Counter
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class OLCAssembler:
    def __init__(self, min_overlap_length=20, min_overlap_identity=0.9, k=15, sketch_size=100, threads=1):
        self.min_overlap_length = min_overlap_length
        self.min_overlap_identity = min_overlap_identity
        self.k = k  # k-mer size for MinHash sketching
        self.sketch_size = sketch_size  # number of k-mers to keep in each sketch
        self.threads = threads
        self.reads = []
        self.read_ids = {}
        self.read_qualities = []
        self.overlap_graph = defaultdict(list)
        self.contained_reads = set()
        self.contigs = []
        self.assembly_metrics = {}
    
    def load_fastq(self, fastq_files):
        """Load reads and quality scores from FASTQ files"""
        read_id = 0
        for fastq_file in fastq_files:
            for record in SeqIO.parse(fastq_file, "fastq"):
                self.reads.append(str(record.seq))
                self.read_ids[read_id] = record.id
                
                # Convert quality scores to probabilities
                quality_scores = record.letter_annotations["phred_quality"]
                error_probs = [10**(-q/10) for q in quality_scores]
                self.read_qualities.append(error_probs)
                
                read_id += 1
        
        print(f"Loaded {len(self.reads)} reads from {len(fastq_files)} files")
    
    def _compute_minhash_sketch(self, sequence):
        """Compute a MinHash sketch of a sequence"""
        # Generate all k-mers from the sequence
        kmers = set()
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i+self.k]
            kmers.add(kmer)
        
        # Hash each k-mer and keep the smallest hashes
        hashes = [(hash(kmer) & 0xFFFFFFFF) for kmer in kmers]  # 32-bit hash values
        hashes.sort()
        
        # Return the sketch (smallest hash values)
        return hashes[:self.sketch_size] if len(hashes) >= self.sketch_size else hashes
    
    def _estimate_jaccard_similarity(self, sketch1, sketch2):
        """Estimate Jaccard similarity between two sequences using their MinHash sketches"""
        # Count number of common hash values
        common = len(set(sketch1) & set(sketch2))
        return common / len(set(sketch1) | set(sketch2)) if (set(sketch1) | set(sketch2)) else 0
    
    def _find_containment(self):
        """Identify reads that are contained within other reads"""
        print("Identifying contained reads...")
        for i in range(len(self.reads)):
            if i in self.contained_reads:
                continue
                
            read_i = self.reads[i]
            for j in range(len(self.reads)):
                if i == j or j in self.contained_reads:
                    continue
                    
                read_j = self.reads[j]
                
                # Check if read_i is contained in read_j
                if len(read_i) < len(read_j):
                    if read_i in read_j:
                        self.contained_reads.add(i)
                        break
        
        print(f"Identified {len(self.contained_reads)} contained reads")
    
    def _find_suffix_prefix_overlap(self, seq1, seq2, min_overlap):
        """Find the best suffix-prefix overlap between seq1 and seq2"""
        best_overlap = 0
        best_identity = 0
        
        for overlap_len in range(min_overlap, min(len(seq1), len(seq2)) + 1):
            suffix = seq1[-overlap_len:]
            prefix = seq2[:overlap_len]
            
            matches = sum(s == p for s, p in zip(suffix, prefix))
            identity = matches / overlap_len
            
            if identity > best_identity:
                best_identity = identity
                best_overlap = overlap_len
        
        return best_overlap, best_identity
    
    def _compute_sketches_batch(self, indices):
        """Compute MinHash sketches for a batch of reads"""
        return [(i, self._compute_minhash_sketch(self.reads[i])) for i in indices]
    
    def _compute_overlaps_batch(self, work_batch):
        """Process a batch of read pairs to find overlaps using MinHash"""
        jaccard_threshold = 0.1  # Adjust based on your data
        results = []
        
        for i, j in work_batch:
            if i == j:
                continue
                
            if i in self.contained_reads or j in self.contained_reads:
                continue
            
            # Quick filter using MinHash Jaccard similarity
            sketch_i = self.sketches.get(i, [])
            sketch_j = self.sketches.get(j, [])
            
            if not sketch_i or not sketch_j:
                continue
                
            similarity = self._estimate_jaccard_similarity(sketch_i, sketch_j)
            
            if similarity >= jaccard_threshold:
                # Compute actual overlap
                overlap_len, identity = self._find_suffix_prefix_overlap(
                    self.reads[i], self.reads[j], self.min_overlap_length
                )
                
                if overlap_len >= self.min_overlap_length and identity >= self.min_overlap_identity:
                    # Add edge with weight based on overlap length and identity
                    weight = overlap_len * identity
                    results.append((i, j, overlap_len, weight))
        
        return results
    
    def _generate_valid_read_pairs(self):
        """Generate all valid pairs of read indices for overlap computation"""
        n = len(self.reads)
        pairs = []
        for i in range(n):
            if i not in self.contained_reads:
                for j in range(n):
                    if j not in self.contained_reads and i != j:
                        pairs.append((i, j))
        return pairs
    
    def _split_work(self, items, num_workers):
        """Split a list of items into approximately equal chunks"""
        chunk_size = max(1, len(items) // num_workers)
        return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]
    
    def _check_containment(self, read_pairs):
        """Helper method to check for read containment in a batch of read pairs"""
        contained = set()
        for i, j in read_pairs:
            if len(self.reads[i]) < len(self.reads[j]):
                if self.reads[i] in self.reads[j]:
                    contained.add(i)
        return contained
    
    def _find_containment_parallel(self):
        """Identify reads that are contained within other reads in parallel"""
        print("Identifying contained reads...")
        
        # Generate all read pairs for comparison
        n = len(self.reads)
        all_pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
        
        # Split work across processes
        if self.threads > 1:
            batches = self._split_work(all_pairs, self.threads)
            
            with multiprocessing.Pool(processes=self.threads) as pool:
                results = pool.map(self._check_containment, batches)
            
            # Combine results
            for contained_set in results:
                self.contained_reads.update(contained_set)
        else:
            self.contained_reads.update(self._check_containment(all_pairs))
        
        print(f"Identified {len(self.contained_reads)} contained reads")
    
    def compute_overlaps(self):
        """Compute overlaps using MinHash sketches for initial filtering, in parallel"""
        # Step 1: Compute MinHash sketches for all reads in parallel
        print(f"Computing sketches for all reads using {self.threads} threads...")
        
        read_indices = list(range(len(self.reads)))
        batches = self._split_work(read_indices, self.threads)
        
        if self.threads > 1:
            with multiprocessing.Pool(processes=self.threads) as pool:
                sketch_results = pool.map(self._compute_sketches_batch, batches)
            
            # Combine sketch results
            self.sketches = {}
            for batch in sketch_results:
                for i, sketch in batch:
                    self.sketches[i] = sketch
        else:
            # Single-threaded fallback
            self.sketches = {i: self._compute_minhash_sketch(self.reads[i]) for i in read_indices}
        
        # Step 2: Identify contained reads in parallel
        self._find_containment_parallel()
        
        # Step 3: Compute overlaps using MinHash filtering in parallel
        print(f"Computing overlaps using MinHash filtering with {self.threads} threads...")
        
        # Generate valid read pairs for comparison
        valid_pairs = self._generate_valid_read_pairs()
        batches = self._split_work(valid_pairs, self.threads)
        
        # Compute overlaps in parallel
        if self.threads > 1:
            with multiprocessing.Pool(processes=self.threads) as pool:
                results = pool.map(self._compute_overlaps_batch, batches)
            
            # Flatten results and build the overlap graph
            all_overlaps = [overlap for batch_result in results for overlap in batch_result]
        else:
            # Single-threaded fallback
            all_overlaps = self._compute_overlaps_batch(valid_pairs)
        
        # Build the overlap graph from results
        for i, j, overlap_len, weight in all_overlaps:
            self.overlap_graph[i].append((j, overlap_len, weight))
        
        # Print some statistics
        edge_count = sum(len(edges) for edges in self.overlap_graph.values())
        print(f"Built overlap graph with {len(self.overlap_graph)} nodes and {edge_count} edges")
    
    def find_paths(self):
        """Identify paths in the overlap graph using greedy approach"""
        print("Finding paths in the overlap graph...")
        visited = set()
        paths = []
        
        # Find all start nodes (in-degree = 0)
        in_degree = defaultdict(int)
        for node in self.overlap_graph:
            for target, _, _ in self.overlap_graph[node]:
                in_degree[target] += 1
        
        start_nodes = [node for node in self.overlap_graph if in_degree[node] == 0]
        
        # If no start nodes found, use nodes with in-degree = 1 (may be circular)
        if not start_nodes:
            start_nodes = [node for node in self.overlap_graph if in_degree[node] <= 1]
        
        # For each start node, find the best path greedily
        for start in start_nodes:
            if start in visited:
                continue
            
            path = [start]
            visited.add(start)
            
            # Extend path greedily by choosing the best overlap at each step
            current = start
            while self.overlap_graph[current]:
                # Sort neighbors by overlap weight in descending order
                next_nodes = sorted(self.overlap_graph[current], key=lambda x: x[2], reverse=True)
                
                best_next = None
                for next_node, _, _ in next_nodes:
                    if next_node not in visited:
                        best_next = next_node
                        break
                
                if best_next is None:
                    break
                    
                path.append(best_next)
                visited.add(best_next)
                current = best_next
            
            if len(path) > 1:
                paths.append(path)
        
        print(f"Found {len(paths)} paths")
        return paths
    
    def generate_layout(self, paths):
        """Generate a layout of reads based on the paths"""
        print("Generating layout...")
        layouts = []
        
        for path in paths:
            layout = []
            
            # For each node in the path, get the read and overlap
            for i in range(len(path) - 1):
                current = path[i]
                next_node = path[i+1]
                
                # Find the overlap between current and next
                overlap_len = 0
                for target, o_len, _ in self.overlap_graph[current]:
                    if target == next_node:
                        overlap_len = o_len
                        break
                
                # Add to layout
                if i == 0:
                    layout.append((current, 0))  # First read starts at position 0
                
                # Calculate position of next read
                read_len = len(self.reads[current])
                next_pos = layout[-1][1] + read_len - overlap_len
                layout.append((next_node, next_pos))
            
            layouts.append(layout)
        
        return layouts
    
    def compute_consensus(self, layouts):
        """Compute a quality-aware consensus sequence for each layout"""
        print("Computing quality-aware consensus sequences...")
        
        for layout_idx, layout in enumerate(layouts):
            # Determine the total length of the contig
            max_pos = 0
            for node, pos in layout:
                max_pos = max(max_pos, pos + len(self.reads[node]))
            
            # Initialize arrays for consensus calculation
            base_votes = [{} for _ in range(max_pos)]
            base_weights = [{} for _ in range(max_pos)]
            
            # Add each read to the consensus, weighting by quality
            for node, pos in layout:
                read_seq = self.reads[node]
                read_qual = self.read_qualities[node]
                
                for i, (base, qual) in enumerate(zip(read_seq, read_qual)):
                    if pos + i < max_pos:
                        weight = 1 - qual  # Higher quality = lower error probability = higher weight
                        
                        if base not in base_votes[pos + i]:
                            base_votes[pos + i][base] = 0
                            base_weights[pos + i][base] = 0
                            
                        base_votes[pos + i][base] += 1
                        base_weights[pos + i][base] += weight
            
            # Determine consensus base at each position using weighted voting
            consensus = []
            for i in range(max_pos):
                if not base_weights[i]:
                    consensus.append('N')  # No coverage at this position
                else:
                    # Choose base with highest weighted vote
                    best_base = max(base_weights[i].items(), key=lambda x: x[1])[0]
                    consensus.append(best_base)
            
            # Create the consensus sequence
            consensus_seq = ''.join(consensus)
            self.contigs.append(consensus_seq)
            
            print(f"Contig {layout_idx+1}: {len(consensus_seq)} bp")
    
    def calculate_assembly_metrics(self):
        """Calculate assembly metrics"""
        print("Calculating assembly metrics...")
        
        # Total assembly length
        total_length = sum(len(contig) for contig in self.contigs)
        
        # Number of contigs
        num_contigs = len(self.contigs)
        
        # GC content
        gc_count = sum(contig.count('G') + contig.count('C') for contig in self.contigs)
        gc_content = (gc_count / total_length) * 100 if total_length > 0 else 0
        
        # Largest contig
        largest_contig = max(len(contig) for contig in self.contigs) if self.contigs else 0
        
        # N50, N90, L50 calculation
        sorted_contigs = sorted(self.contigs, key=len, reverse=True)
        sorted_lengths = [len(contig) for contig in sorted_contigs]
        
        cumulative_length = 0
        n50 = 0
        n90 = 0
        l50 = 0
        
        for i, length in enumerate(sorted_lengths):
            cumulative_length += length
            
            if not n50 and cumulative_length >= total_length * 0.5:
                n50 = length
                l50 = i + 1
                
            if not n90 and cumulative_length >= total_length * 0.9:
                n90 = length
                break
        
        self.assembly_metrics = {
            "Total assembly length": f"{total_length:,} bp",
            "Number of contigs": num_contigs,
            "GC content": f"{gc_content:.2f}%",
            "Largest contig": f"{largest_contig:,} bp",
            "N50": f"{n50:,} bp",
            "N90": f"{n90:,} bp",
            "L50": l50
        }
        
        return self.assembly_metrics
    
    def write_contigs(self, output_file):
        """Write contigs to a FASTA file"""
        records = []
        for i, contig in enumerate(self.contigs):
            record = SeqRecord(
                Seq(contig),
                id=f"contig_{i+1}",
                description=f"length={len(contig)}"
            )
            records.append(record)
        
        SeqIO.write(records, output_file, "fasta")
        print(f"Wrote {len(records)} contigs to {output_file}")
    
    def write_metrics(self, output_file):
        """Write assembly metrics to a text file"""
        with open(output_file, "w") as f:
            f.write("=" * 70 + "\n")
            f.write(f"{'ASSEMBLY METRICS':^70}\n")
            f.write("=" * 70 + "\n")
            
            for metric, value in self.assembly_metrics.items():
                f.write(f"{metric:<30}: {value}\n")
            
            f.write("\nNote: Metrics that require a reference genome are not included.\n")
        
        print(f"Wrote assembly metrics to {output_file}")
    
    def assemble(self, fastq_files, output_prefix):
        """Run the complete assembly pipeline"""
        output_fasta = f"{output_prefix}.fasta"
        output_metrics = f"{output_prefix}_metrics.txt"
        
        # Step 1: Load FASTQ files with quality scores
        self.load_fastq(fastq_files)
        
        # Step 2: Compute overlaps using MinHash and build graph
        self.compute_overlaps()
        
        # Step 3: Find paths
        paths = self.find_paths()
        
        # Step 4: Generate layout
        layouts = self.generate_layout(paths)
        
        # Step 5: Compute quality-aware consensus
        self.compute_consensus(layouts)
        
        # Step 6: Calculate assembly metrics
        self.calculate_assembly_metrics()
        
        # Step 7: Write outputs
        self.write_contigs(output_fasta)
        self.write_metrics(output_metrics)
        
        return self.assembly_metrics


def main():
    parser = argparse.ArgumentParser(description='OLC Assembly Algorithm')
    parser.add_argument('--fastq', nargs='+', required=True, help='Input FASTQ files')
    parser.add_argument('--min-overlap', type=int, default=20, help='Minimum overlap length')
    parser.add_argument('--min-identity', type=float, default=0.9, help='Minimum overlap identity')
    parser.add_argument('--kmer-size', type=int, default=15, help='K-mer size for MinHash')
    parser.add_argument('--sketch-size', type=int, default=100, help='Sketch size for MinHash')
    parser.add_argument('--threads', type=int, default=1, help='Number of CPU cores to use')
    parser.add_argument('--output', required=True, help='Output prefix for files')
    
    args = parser.parse_args()
    
    # Use all available cores if threads=0
    if args.threads == 0:
        args.threads = multiprocessing.cpu_count()
    
    assembler = OLCAssembler(
        min_overlap_length=args.min_overlap,
        min_overlap_identity=args.min_identity,
        k=args.kmer_size,
        sketch_size=args.sketch_size,
        threads=args.threads
    )
    
    metrics = assembler.assemble(args.fastq, args.output)
    
    # Print metrics to console
    print("\nASSEMBLY METRICS:")
    for metric, value in metrics.items():
        print(f"{metric}: {value}")


if __name__ == "__main__":
    main()