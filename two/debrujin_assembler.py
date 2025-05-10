#!/usr/bin/env python3

# This script implements a de Bruijn graph assembler for DNA sequences.
# It reads sequences from FASTQ files, constructs a de Bruijn graph,
# and finds Eulerian paths to generate contigs. The results can be output
# to FASTA format and the graph can be exported to GFA format for visualization
# with Bandage.
# The script uses the NetworkX library for graph operations and Biopython
# for sequence handling.

#usage: python debrujin_assembler.py -k 21 -i reads1.fastq reads2.fastq -o contigs.fasta -g graph.gfa

import os
import argparse
import networkx as nx
from collections import defaultdict, deque
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class DeBruijnGraphAssembler:
    def __init__(self, k):
        """Initialize the assembler with k-mer size."""
        self.k = k
        self.graph = nx.MultiDiGraph()
        self.edge_seqs = {}  # Maps edge IDs to sequences
        self.edge_counts = defaultdict(int)  # Tracks coverage
        self.contigs = []  # Store assembled contigs
        
    def add_read(self, read_seq):
        """Add a read to the de Bruijn graph."""
        # Skip reads shorter than k
        if len(read_seq) < self.k:
            return
            
        # Generate k-mers and build graph
        for i in range(len(read_seq) - self.k + 1):
            kmer = read_seq[i:i+self.k]
            # Create prefix and suffix (k-1 mers)
            prefix = kmer[:-1]
            suffix = kmer[1:]
            
            # Add nodes and edge to graph
            if not self.graph.has_edge(prefix, suffix):
                edge_id = len(self.edge_seqs)
                self.graph.add_edge(prefix, suffix, id=edge_id)
                self.edge_seqs[edge_id] = kmer[-1]  # Store the last base
            else:
                # Find the edge ID
                for edge_data in self.graph.get_edge_data(prefix, suffix).values():
                    edge_id = edge_data['id']
                    
            # Increment edge count (coverage)
            self.edge_counts[edge_id] += 1
    
    def load_fastq(self, fastq_file):
        """Load reads from FASTQ file."""
        with open(fastq_file, 'r') as f:
            for record in SeqIO.parse(f, 'fastq'):
                self.add_read(str(record.seq))
                
    def load_fastq_files(self, fastq_files):
        """Load reads from multiple FASTQ files."""
        for fastq_file in fastq_files:
            self.load_fastq(fastq_file)
        print(f"Loaded {len(self.edge_seqs)} unique k-mers from {len(fastq_files)} files")
        print(f"Graph has {self.graph.number_of_nodes()} nodes and {self.graph.number_of_edges()} edges")
    
    def hierholzer_algorithm(self, start_node=None):
        """Find Eulerian paths using Hierholzer's algorithm."""
        if not self.graph.nodes():
            return []
            
        # Make a copy of the graph to work with
        work_graph = self.graph.copy()
        
        paths = []
        
        # Process all connected components
        while work_graph.nodes():
            # If no start node provided or not in graph, choose one
            if start_node is None or start_node not in work_graph:
                # Find a node with unbalanced degree if possible
                for node in work_graph.nodes():
                    if work_graph.in_degree(node) != work_graph.out_degree(node):
                        if work_graph.out_degree(node) > work_graph.in_degree(node):
                            start_node = node
                            break
                
                # If all nodes are balanced, pick any node
                if start_node is None or start_node not in work_graph:
                    start_node = list(work_graph.nodes())[0]
            
            # Apply Hierholzer's algorithm
            path = []
            stack = [start_node]
            
            while stack:
                current = stack[-1]
                
                # If current node has no outgoing edges
                if work_graph.out_degree(current) == 0:
                    path.append(current)
                    stack.pop()
                else:
                    # Get the next edge and remove it from the working graph
                    for successor in list(work_graph.successors(current)):
                        edge_data = list(work_graph.get_edge_data(current, successor).values())[0]
                        edge_id = edge_data['id']
                        work_graph.remove_edge(current, successor, key=list(work_graph.get_edge_data(current, successor))[0])
                        stack.append(successor)
                        break
            
            # Reverse the path to get the correct order
            path.reverse()
            
            if len(path) > 1:  # Only keep non-trivial paths
                paths.append(path)
            
            # Remove isolated nodes
            work_graph.remove_nodes_from(list(nx.isolates(work_graph)))
            
            # Reset start_node for next component
            start_node = None
            
        return paths
    
    def paths_to_contigs(self, paths):
        """Convert paths to contigs."""
        contigs = []
        
        for i, path in enumerate(paths):
            if len(path) < 2:
                continue
                
            # Start with the first k-1 mer
            contig = path[0]
            
            # Add the last base of each subsequent k-mer
            for j in range(1, len(path)):
                # Find the edge connecting these nodes
                edge_datas = self.graph.get_edge_data(path[j-1], path[j])
                if edge_datas:  # There might be multiple edges
                    for key, edge_data in edge_datas.items():
                        edge_id = edge_data['id']
                        contig += self.edge_seqs[edge_id]
                        break
            
            # Create a SeqRecord for this contig
            record = SeqRecord(
                Seq(contig),
                id=f"contig_{i+1}",
                description=f"length={len(contig)}"
            )
            contigs.append(record)
            
        return contigs
    
    def calculate_assembly_metrics(self):
        """Calculate and return reference-free assembly metrics."""
        if not self.contigs:
            return None
            
        # Extract contig lengths and sort in descending order
        contig_lengths = sorted([len(contig.seq) for contig in self.contigs], reverse=True)
        
        # Total assembly length
        total_length = sum(contig_lengths)
        
        # Number of contigs
        num_contigs = len(contig_lengths)
        
        # Largest contig
        largest_contig = contig_lengths[0] if contig_lengths else 0
        
        # Calculate GC content
        gc_count = 0
        total_bases = 0
        for contig in self.contigs:
            sequence = str(contig.seq).upper()
            gc_count += sequence.count('G') + sequence.count('C')
            total_bases += len(sequence)
        gc_content = (gc_count / total_bases * 100) if total_bases > 0 else 0
        
        # Calculate N50 and N90
        cumulative_length = 0
        n50 = n90 = 0
        l50 = 0
        
        for i, length in enumerate(contig_lengths):
            cumulative_length += length
            
            # N50 calculation
            if not n50 and cumulative_length >= total_length * 0.5:
                n50 = length
                l50 = i + 1
                
            # N90 calculation
            if not n90 and cumulative_length >= total_length * 0.9:
                n90 = length
                
        metrics = {
            'total_length': total_length,
            'num_contigs': num_contigs,
            'gc_content': gc_content,
            'largest_contig': largest_contig,
            'n50': n50,
            'n90': n90,
            'l50': l50
        }
        
        return metrics
    
    def print_assembly_metrics(self):
        """Print assembly metrics to the terminal."""
        metrics = self.calculate_assembly_metrics()
        
        if not metrics:
            print("No contigs were assembled, cannot calculate metrics.")
            return
            
        print("\n" + "="*60)
        print(" "*20 + "ASSEMBLY METRICS")
        print("="*60)
        print(f"Total assembly length        : {metrics['total_length']:,} bp")
        print(f"Number of contigs            : {metrics['num_contigs']}")
        print(f"GC content                   : {metrics['gc_content']:.2f}%")
        print(f"Largest contig               : {metrics['largest_contig']:,} bp")
        print(f"N50                          : {metrics['n50']:,} bp")
        print(f"N90                          : {metrics['n90']:,} bp")
        print(f"L50                          : {metrics['l50']}")
        print("="*60)
        print("Note: Metrics that require a reference genome are not included.")
        print("="*60 + "\n")
        
    def assemble_and_output(self, output_fasta, output_gfa=None):
        """Run assembly and output contigs to FASTA and graph to GFA."""
        # Find Eulerian paths
        paths = self.hierholzer_algorithm()
        print(f"Found {len(paths)} potential contigs")
        
        # Convert paths to contigs
        self.contigs = self.paths_to_contigs(paths)
        
        # Calculate and display assembly metrics
        self.print_assembly_metrics()
        
        # Write contigs to FASTA
        with open(output_fasta, 'w') as f:
            SeqIO.write(self.contigs, f, 'fasta')
        print(f"Wrote {len(self.contigs)} contigs to {output_fasta}")
        
        # Export graph to GFA if requested
        if output_gfa:
            self.export_to_gfa(output_gfa, self.contigs)
            print(f"Wrote graph to {output_gfa}")
    
    def export_to_gfa(self, output_gfa, contigs):
        """Export the graph to GFA format for Bandage visualization."""
        with open(output_gfa, 'w') as f:
            # Write header
            f.write("H\tVN:Z:1.0\n")
            
            # Write segments (nodes)
            node_ids = {}
            for i, node in enumerate(self.graph.nodes()):
                node_ids[node] = i+1
                f.write(f"S\t{i+1}\t{node}\n")
            
            # Write links (edges)
            for u, v, data in self.graph.edges(data=True):
                edge_id = data['id']
                coverage = self.edge_counts[edge_id]
                f.write(f"L\t{node_ids[u]}\t+\t{node_ids[v]}\t+\t{self.k-1}M\tRC:i:{coverage}\n")
            
            # Write paths (contigs)
            for i, contig in enumerate(contigs):
                contig_id = contig.id
                contig_len = len(contig.seq)
                # Convert path to segments
                path_nodes = []
                for j in range(contig_len - self.k + 2):
                    if j < contig_len - self.k + 1:
                        kmer_prefix = str(contig.seq[j:j+self.k-1])
                        if kmer_prefix in node_ids:
                            path_nodes.append(f"{node_ids[kmer_prefix]}+")
                
                if path_nodes:
                    f.write(f"P\t{contig_id}\t{','.join(path_nodes)}\t*\n")


def main():
    parser = argparse.ArgumentParser(description='De Bruijn Graph Assembly')
    parser.add_argument('-k', type=int, required=True, help='k-mer size')
    parser.add_argument('-i', '--input', nargs='+', required=True, help='Input FASTQ file(s)')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    parser.add_argument('-g', '--gfa', help='Output GFA file for Bandage visualization')
    
    args = parser.parse_args()
    
    # Create assembler
    assembler = DeBruijnGraphAssembler(args.k)
    
    # Load reads
    assembler.load_fastq_files(args.input)
    
    # Assemble and output
    assembler.assemble_and_output(args.output, args.gfa)
    

if __name__ == "__main__":
    main()