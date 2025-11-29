"""Sequence validator module."""
from typing import List
from modules.data_models import Sequence, ValidationResult


class SequenceValidator:
    """Validator for cleaning and validating biological sequences."""
    
    # Valid characters for different sequence types
    NUCLEOTIDE_CHARS = set('ACGTRYSWKMBDHVN')  # DNA with ambiguity codes
    PROTEIN_CHARS = set('ACDEFGHIKLMNPQRSTVWY')  # Standard amino acids
    
    def __init__(self, min_length: int = 10):
        """
        Initialize validator.
        
        Args:
            min_length: Minimum sequence length for valid sequences
        """
        self.min_length = min_length
    
    def validate_sequences(self, sequences: List[Sequence]) -> ValidationResult:
        """
        Validate and clean a list of sequences.
        
        Args:
            sequences: List of Sequence objects to validate
            
        Returns:
            ValidationResult with cleaned sequences and summary
        """
        if not sequences:
            return ValidationResult(
                valid_sequences=[],
                invalid_sequences=[],
                cleaning_summary={},
                sequence_type=""
            )
        
        # Detect sequence type
        seq_type = self.detect_sequence_type([seq.sequence for seq in sequences])
        
        valid_sequences = []
        invalid_sequences = []
        cleaning_summary = {
            'removed_spaces': 0,
            'removed_numbers': 0,
            'removed_special_chars': 0,
            'sequences_too_short': 0
        }
        
        for seq in sequences:
            # Clean the sequence
            cleaned_seq, clean_stats = self._clean_sequence_with_stats(seq.sequence, seq_type)
            
            # Update cleaning summary
            cleaning_summary['removed_spaces'] += clean_stats['spaces']
            cleaning_summary['removed_numbers'] += clean_stats['numbers']
            cleaning_summary['removed_special_chars'] += clean_stats['special_chars']
            
            # Check minimum length
            if self.check_minimum_length(cleaned_seq, self.min_length):
                valid_sequences.append(Sequence(
                    id=seq.id,
                    sequence=cleaned_seq,
                    description=seq.description
                ))
            else:
                cleaning_summary['sequences_too_short'] += 1
                invalid_sequences.append(seq)
        
        return ValidationResult(
            valid_sequences=valid_sequences,
            invalid_sequences=invalid_sequences,
            cleaning_summary=cleaning_summary,
            sequence_type=seq_type
        )
    
    def clean_sequence(self, seq: str, seq_type: str = None) -> str:
        """
        Clean a sequence by removing invalid characters.
        
        Args:
            seq: Sequence string to clean
            seq_type: Type of sequence ('nucleotide' or 'protein'), auto-detect if None
            
        Returns:
            Cleaned sequence string
        """
        if seq_type is None:
            seq_type = self.detect_sequence_type([seq])
        
        return self.remove_invalid_characters(seq, seq_type)
    
    def _clean_sequence_with_stats(self, seq: str, seq_type: str):
        """
        Clean sequence and return statistics.
        
        Args:
            seq: Sequence string to clean
            seq_type: Type of sequence
            
        Returns:
            Tuple of (cleaned_sequence, stats_dict)
        """
        stats = {'spaces': 0, 'numbers': 0, 'special_chars': 0}
        
        # Count removals
        stats['spaces'] = seq.count(' ') + seq.count('\t') + seq.count('\n') + seq.count('\r')
        stats['numbers'] = sum(1 for c in seq if c.isdigit())
        
        # Clean the sequence
        cleaned = self.remove_invalid_characters(seq, seq_type)
        
        # Calculate special chars removed (total removed - spaces - numbers)
        total_removed = len(seq) - len(cleaned)
        stats['special_chars'] = total_removed - stats['spaces'] - stats['numbers']
        
        return cleaned, stats
    
    def detect_sequence_type(self, sequences: List[str]) -> str:
        """
        Detect whether sequences are nucleotide or protein.
        
        Args:
            sequences: List of sequence strings
            
        Returns:
            'nucleotide' or 'protein'
        """
        if not sequences:
            return 'nucleotide'
        
        # Sample sequences to determine type
        sample_size = min(len(sequences), 5)
        sample_seqs = sequences[:sample_size]
        
        nucleotide_score = 0
        protein_score = 0
        
        for seq in sample_seqs:
            seq_upper = seq.upper()
            # Remove whitespace and numbers for analysis
            seq_clean = ''.join(c for c in seq_upper if c.isalpha())
            
            if not seq_clean:
                continue
            
            # Count nucleotide-specific characters
            nucleotide_count = sum(1 for c in seq_clean if c in self.NUCLEOTIDE_CHARS)
            # Count protein-specific characters (not in nucleotide set)
            protein_specific = sum(1 for c in seq_clean if c in self.PROTEIN_CHARS and c not in self.NUCLEOTIDE_CHARS)
            
            # Calculate percentages
            if len(seq_clean) > 0:
                nuc_percent = nucleotide_count / len(seq_clean)
                
                # If mostly nucleotide characters, likely DNA
                if nuc_percent > 0.9:
                    nucleotide_score += 1
                # If has protein-specific characters, likely protein
                elif protein_specific > 0:
                    protein_score += 1
                # If moderate nucleotide content, could be either
                elif nuc_percent > 0.5:
                    nucleotide_score += 0.5
        
        # Determine type based on scores
        if protein_score > nucleotide_score:
            return 'protein'
        else:
            return 'nucleotide'
    
    def remove_invalid_characters(self, seq: str, seq_type: str) -> str:
        """
        Remove invalid characters from a sequence.
        
        Args:
            seq: Sequence string
            seq_type: 'nucleotide' or 'protein'
            
        Returns:
            Cleaned sequence with only valid characters
        """
        seq_upper = seq.upper()
        
        if seq_type == 'nucleotide':
            valid_chars = self.NUCLEOTIDE_CHARS
        else:
            valid_chars = self.PROTEIN_CHARS
        
        # Keep only valid characters, preserving order
        cleaned = ''.join(c for c in seq_upper if c in valid_chars)
        
        return cleaned
    
    def check_minimum_length(self, seq: str, min_length: int = 10) -> bool:
        """
        Check if sequence meets minimum length requirement.
        
        Args:
            seq: Sequence string
            min_length: Minimum required length
            
        Returns:
            True if sequence is long enough, False otherwise
        """
        return len(seq) >= min_length

