"""Sequence aligner module."""
from typing import List
from Bio import Align
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from modules.data_models import Sequence


class AlignmentError(Exception):
    """Exception raised for alignment errors."""
    pass


class SequenceAligner:
    """Aligner for multiple sequence alignment."""
    
    def __init__(self):
        """Initialize the aligner."""
        self.aligner = Align.PairwiseAligner()
        # Set alignment parameters
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -0.5
    
    def align_sequences(self, sequences: List[Sequence]) -> MultipleSeqAlignment:
        """
        Perform multiple sequence alignment.
        
        Args:
            sequences: List of Sequence objects to align
            
        Returns:
            Bio.Align.MultipleSeqAlignment object
            
        Raises:
            AlignmentError: If alignment fails
        """
        if not sequences:
            raise AlignmentError("No sequences provided for alignment")
        
        if len(sequences) < 2:
            raise AlignmentError("At least 2 sequences required for alignment")
        
        # Check if all sequences are empty
        if all(len(seq.sequence) == 0 for seq in sequences):
            raise AlignmentError("All sequences are empty")
        
        try:
            # For small datasets, use progressive alignment
            if len(sequences) <= 10:
                return self._progressive_alignment(sequences)
            else:
                # For larger datasets, still use progressive but with optimization
                return self._progressive_alignment(sequences)
        except Exception as e:
            raise AlignmentError(
                f"Alignment failed: The sequences are too different to align reliably. "
                f"Try using sequences that are more closely related. Error: {str(e)}"
            )
    
    def _progressive_alignment(self, sequences: List[Sequence]) -> MultipleSeqAlignment:
        """
        Perform progressive alignment using pairwise alignments.
        
        Args:
            sequences: List of Sequence objects
            
        Returns:
            MultipleSeqAlignment object
        """
        # Start with the first sequence
        aligned_seqs = [sequences[0].sequence]
        aligned_ids = [sequences[0].id]
        
        # Progressively add each sequence
        for i in range(1, len(sequences)):
            new_seq = sequences[i].sequence
            new_id = sequences[i].id
            
            # Align new sequence to the first sequence in the alignment
            reference_seq = sequences[0].sequence
            
            # Perform pairwise alignment
            alignments = self.aligner.align(reference_seq, new_seq)
            
            if not alignments:
                # If no alignment found, add sequence with gaps
                max_len = max(len(s) for s in aligned_seqs)
                aligned_seqs.append(new_seq + '-' * (max_len - len(new_seq)))
                aligned_ids.append(new_id)
                continue
            
            # Get best alignment
            best_alignment = alignments[0]
            
            # Extract aligned sequences
            aligned_ref = str(best_alignment[0])
            aligned_new = str(best_alignment[1])
            
            # Update all sequences in the alignment to match the new length
            if i == 1:
                # First addition - update the reference sequence
                aligned_seqs[0] = aligned_ref
            else:
                # Add gaps to existing sequences where needed
                self._add_gaps_to_alignment(aligned_seqs, aligned_ref, sequences[0].sequence)
            
            # Add the new aligned sequence
            aligned_seqs.append(aligned_new)
            aligned_ids.append(new_id)
        
        # Ensure all sequences have the same length
        max_len = max(len(s) for s in aligned_seqs)
        aligned_seqs = [s + '-' * (max_len - len(s)) for s in aligned_seqs]
        
        # Create SeqRecord objects
        seq_records = []
        for seq_id, seq_str in zip(aligned_ids, aligned_seqs):
            seq_records.append(SeqRecord(Seq(seq_str), id=seq_id, description=""))
        
        # Create MultipleSeqAlignment
        return MultipleSeqAlignment(seq_records)
    
    def _add_gaps_to_alignment(self, aligned_seqs: List[str], new_aligned_ref: str, original_ref: str):
        """
        Add gaps to existing alignment based on new alignment.
        
        Args:
            aligned_seqs: List of currently aligned sequences
            new_aligned_ref: Newly aligned reference sequence
            original_ref: Original reference sequence
        """
        # Find where gaps were inserted in the reference
        gap_positions = []
        ref_pos = 0
        
        for i, char in enumerate(new_aligned_ref):
            if char == '-':
                gap_positions.append(i)
            else:
                ref_pos += 1
        
        # Insert gaps at the same positions in all aligned sequences
        for i in range(len(aligned_seqs)):
            seq_list = list(aligned_seqs[i])
            for gap_pos in reversed(gap_positions):  # Insert from right to left
                if gap_pos <= len(seq_list):
                    seq_list.insert(gap_pos, '-')
            aligned_seqs[i] = ''.join(seq_list)
    
    def get_alignment_method(self) -> str:
        """
        Get the name of the alignment method.
        
        Returns:
            Method name string
        """
        return "Progressive Pairwise Alignment"

