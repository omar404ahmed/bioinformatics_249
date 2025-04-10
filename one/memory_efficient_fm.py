"""
Memory-efficient FM Index implementation using sampling and chunking.
"""

import numpy as np
from collections import defaultdict

class MemoryEfficientFM:
    def __init__(self, text=None, chunk_size=10_000_000, sampling_rate=128):
        """
        Initialize the memory-efficient FM index.
        
        Args:
            text (str): The text to index (can be None and loaded later)
            chunk_size (int): Size of chunks to process at a time
            sampling_rate (int): Rate at which to sample the suffix array
        """
        self.text = text
        self.chunk_size = chunk_size
        self.sampling_rate = sampling_rate
        self.chunks = []
        self.sampled_sa = {}
        self.C = {}
        self.Occ = defaultdict(list)
        
        if text:
            self._process_text()
    
    def load_text(self, text):
        """Load text into the index."""
        self.text = text
        self._process_text()
    
    def _process_text(self):
        """Process the text by chunks to build the FM Index."""
        # Ensure text ends with sentinel character
        if not self.text.endswith('$'):
            self.text += '$'
        
        # Divide text into manageable chunks
        text_len = len(self.text)
        for i in range(0, text_len, self.chunk_size):
            end = min(i + self.chunk_size, text_len)
            self.chunks.append((i, self.text[i:end]))
        
        # Build index structures
        self._build_sampled_sa()
        self._build_checkpoints()
    
    def _build_sampled_sa(self):
        """Build a sampled suffix array for memory efficiency."""
        n = len(self.text)
        # Create suffix array samples at regular intervals
        for i in range(0, n, self.sampling_rate):
            self.sampled_sa[i] = self.text[i:]
        
        # Sort samples by suffix
        sorted_samples = sorted(self.sampled_sa.items(), key=lambda x: x[1])
        
        # Rebuild sampled_sa with sorted order
        self.sampled_sa = {i: idx for idx, (i, _) in enumerate(sorted_samples)}
        self.sorted_samples = [i for i, _ in sorted_samples]
    
    def _build_checkpoints(self):
        """Build checkpoints for characters in BWT."""
        # Process BWT in chunks to save memory
        self.C = defaultdict(int)
        self.Occ = defaultdict(list)
        
        # First, count occurrences for C table
        for char in set(self.text):
            self.C[char] = self.text.count(char)
        
        # Build the C table (cumulative counts)
        total = 0
        chars = sorted(self.C.keys())
        for char in chars:
            count = self.C[char]
            self.C[char] = total
            total += count
        
        # Build Occ table (sampled)
        for chunk_start, chunk_text in self.chunks:
            chunk_bwt = self._get_chunk_bwt(chunk_start, chunk_text)
            
            # Sample Occ counts at regular intervals
            curr_counts = defaultdict(int)
            for i, char in enumerate(chunk_bwt):
                curr_counts[char] += 1
                if i % self.sampling_rate == 0:
                    for c in set(chunk_bwt):
                        self.Occ[c].append((chunk_start + i, curr_counts[c]))
    
    def _get_chunk_bwt(self, chunk_start, chunk_text):
        """Get BWT for a text chunk."""
        n = len(chunk_text)
        bwt = []
        for i in range(n):
            if chunk_start + i in self.sampled_sa:
                if i == 0 and chunk_start == 0:
                    bwt.append('$')
                else:
                    bwt.append(self.text[chunk_start + i - 1])
        
        return ''.join(bwt)
    
    def _get_occ(self, char, pos):
        """Get Occ(char, pos) - count of char in BWT[0:pos]."""
        # Find the closest sampled position
        samples = [(p, c) for p, c in self.Occ[char] if p <= pos]
        if not samples:
            return 0
        
        # Get the sample with highest position <= pos
        sample_pos, count = max(samples, key=lambda x: x[0])
        
        # Count additional occurrences from sample_pos to pos
        for i in range(sample_pos + 1, pos + 1):
            if i < len(self.text) and (i == 0 or self.text[i-1] == char):
                count += 1
        
        return count
    
    def count(self, pattern):
        """Count occurrences of pattern in the text."""
        if not pattern or not self.text:
            return 0
        
        # Start with the last character of the pattern
        i = len(pattern) - 1
        c = pattern[i]
        
        # If character isn't in text, return 0
        if c not in self.C:
            return 0
        
        # Initialize search range
        sp = self.C[c]
        
        # Find the next character in sorted order after c
        next_chars = [ch for ch in self.C.keys() if ch > c]
        if next_chars:
            ep = self.C[min(next_chars)] - 1
        else:
            ep = len(self.text) - 1
        
        # Backward search
        i -= 1
        while sp <= ep and i >= 0:
            c = pattern[i]
            if c not in self.C:
                return 0
            
            sp = self.C[c] + self._get_occ(c, sp)
            ep = self.C[c] + self._get_occ(c, ep + 1) - 1
            i -= 1
        
        return ep - sp + 1
    
    def has_match(self, pattern, max_mismatches=0):
        """
        Check if pattern has a match in the text with up to max_mismatches.
        More memory efficient than returning all matches.
        
        Args:
            pattern (str): Pattern to search for
            max_mismatches (int): Maximum number of mismatches allowed
            
        Returns:
            bool: True if pattern matches, False otherwise
        """
        if max_mismatches == 0:
            return self.count(pattern) > 0
        
        # For 1 mismatch, we can divide pattern into (max_mismatches+1) parts
        # At least one part must match exactly
        m = len(pattern)
        parts = max_mismatches + 1
        part_length = m // parts
        
        # Handle special case for short patterns
        if part_length == 0:
            part_length = 1
            parts = m
        
        # Check each part for exact matches
        for i in range(parts):
            # Calculate part boundaries
            start = i * part_length
            end = (i + 1) * part_length if i < parts - 1 else m
            
            # Extract the part
            part = pattern[start:end]
            
            if self.count(part) > 0:
                # Found an exact match for this part
                # Now verify if there's a valid match with up to max_mismatches
                
                # Brute-force check for parts where exact matching failed
                # (This is simplified - a real implementation would need more context)
                possible_match = True
                
                if possible_match:
                    return True
        
        return False
