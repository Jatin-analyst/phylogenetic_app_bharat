"""Property-based tests for aligner module."""
from hypothesis import given, strategies as st, assume, settings
from modules.aligner import SequenceAligner, AlignmentError
from modules.data_models import Sequence
from Bio.Align import MultipleSeqAlignment


# Custom strategies
@st.composite
def sequence_list(draw, min_seqs=2, max_seqs=5):
    """Generate a list of related sequences."""
    num_seqs = draw(st.integers(min_value=min_seqs, max_value=max_seqs))
    base_length = draw(st.integers(min_value=20, max_value=50))
    
    # Generate a base sequence
    base_seq = ''.join(draw(st.lists(
        st.sampled_from('ACGT'),
        min_size=base_length,
        max_size=base_length
    )))
    
    sequences = [Sequence(id="Seq_0", sequence=base_seq)]
    
    # Generate related sequences with mutations
    for i in range(1, num_seqs):
        # Create a variant of the base sequence
        variant = list(base_seq)
        
        # Add some mutations (substitutions)
        num_mutations = draw(st.integers(min_value=0, max_value=min(5, len(variant) // 4)))
        for _ in range(num_mutations):
            pos = draw(st.integers(min_value=0, max_value=len(variant) - 1))
            variant[pos] = draw(st.sampled_from('ACGT'))
        
        # Optionally add/remove a few characters
        if draw(st.booleans()) and len(variant) > 10:
            # Remove a character
            pos = draw(st.integers(min_value=0, max_value=len(variant) - 1))
            variant.pop(pos)
        
        if draw(st.booleans()):
            # Add a character
            pos = draw(st.integers(min_value=0, max_value=len(variant)))
            variant.insert(pos, draw(st.sampled_from('ACGT')))
        
        sequences.append(Sequence(id=f"Seq_{i}", sequence=''.join(variant)))
    
    return sequences


# Feature: phylogenetics-app, Property 6: Alignment preserves sequences and identifiers
# Validates: Requirements 3.1, 3.4
@given(sequence_list(min_seqs=2, max_seqs=5))
@settings(deadline=None)  # Disable deadline for alignment operations
def test_alignment_preserves_count_and_ids(sequences):
    """
    Property: For any set of valid sequences, alignment should produce
    an output with the same number of sequences and identifiers.
    """
    aligner = SequenceAligner()
    
    alignment = aligner.align_sequences(sequences)
    
    # Check sequence count is preserved
    assert len(alignment) == len(sequences), \
        f"Expected {len(sequences)} sequences, got {len(alignment)}"
    
    # Check IDs are preserved
    input_ids = [seq.id for seq in sequences]
    output_ids = [record.id for record in alignment]
    
    assert output_ids == input_ids, \
        f"IDs not preserved. Expected {input_ids}, got {output_ids}"


@given(sequence_list(min_seqs=2, max_seqs=5))
@settings(deadline=None)
def test_alignment_produces_equal_length_sequences(sequences):
    """
    Property: For any alignment, all sequences should have the same length
    (due to gap insertion).
    """
    aligner = SequenceAligner()
    
    alignment = aligner.align_sequences(sequences)
    
    # Get all sequence lengths
    lengths = [len(record.seq) for record in alignment]
    
    # All lengths should be equal
    assert len(set(lengths)) == 1, \
        f"Not all sequences have same length: {lengths}"
    
    # Length should be at least as long as the longest input sequence
    max_input_len = max(len(seq.sequence) for seq in sequences)
    assert lengths[0] >= max_input_len, \
        f"Aligned length {lengths[0]} is shorter than longest input {max_input_len}"


@given(st.lists(st.text(alphabet='ACGT', min_size=15, max_size=30), min_size=2, max_size=5))
@settings(deadline=None)
def test_alignment_preserves_sequence_content(sequences):
    """
    Property: Alignment should preserve the original sequence characters
    (gaps may be added but no characters should be changed).
    """
    # Create Sequence objects
    seq_objects = [Sequence(id=f"Seq_{i}", sequence=seq) for i, seq in enumerate(sequences)]
    
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(seq_objects)
    
    # Check each aligned sequence
    for i, record in enumerate(alignment):
        aligned_seq = str(record.seq)
        original_seq = sequences[i]
        
        # Remove gaps from aligned sequence
        ungapped = aligned_seq.replace('-', '')
        
        # Should match original sequence
        assert ungapped == original_seq, \
            f"Sequence {i}: aligned (ungapped) '{ungapped}' != original '{original_seq}'"


@given(st.just([
    Sequence(id="Seq_1", sequence="ACGT"),
    Sequence(id="Seq_2", sequence="ACGT")
]))
@settings(deadline=None)
def test_identical_sequences_align_without_gaps(identical_seqs):
    """
    Property: Identical sequences should align without gaps.
    """
    aligner = SequenceAligner()
    
    alignment = aligner.align_sequences(identical_seqs)
    
    # Check no gaps were added
    for record in alignment:
        assert '-' not in str(record.seq), \
            "Identical sequences should not have gaps"



# Feature: phylogenetics-app, Property 7: Aligned output is valid FASTA
# Validates: Requirements 3.2
@given(sequence_list(min_seqs=2, max_seqs=5))
@settings(deadline=None)
def test_aligned_fasta_export_is_valid(sequences):
    """
    Property: For any alignment result, exporting to FASTA format should
    produce a valid FASTA file that can be parsed back.
    """
    from io import StringIO
    from Bio import AlignIO, SeqIO
    
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    # Export to FASTA format
    output = StringIO()
    AlignIO.write(alignment, output, "fasta")
    fasta_content = output.getvalue()
    
    # Parse it back
    input_stream = StringIO(fasta_content)
    parsed_records = list(SeqIO.parse(input_stream, "fasta"))
    
    # Check same number of sequences
    assert len(parsed_records) == len(sequences), \
        f"Expected {len(sequences)} sequences after round-trip, got {len(parsed_records)}"
    
    # Check IDs are preserved
    original_ids = [seq.id for seq in sequences]
    parsed_ids = [record.id for record in parsed_records]
    assert parsed_ids == original_ids, \
        "IDs not preserved in FASTA round-trip"


@given(sequence_list(min_seqs=2, max_seqs=3))
@settings(deadline=None)
def test_fasta_export_preserves_alignment_length(sequences):
    """
    Property: FASTA export should preserve the aligned sequence lengths.
    """
    from io import StringIO
    from Bio import AlignIO, SeqIO
    
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    # Get original aligned length
    original_length = len(alignment[0].seq)
    
    # Export and re-import
    output = StringIO()
    AlignIO.write(alignment, output, "fasta")
    fasta_content = output.getvalue()
    
    input_stream = StringIO(fasta_content)
    parsed_records = list(SeqIO.parse(input_stream, "fasta"))
    
    # Check all sequences have the same length as original
    for record in parsed_records:
        assert len(record.seq) == original_length, \
            f"Length mismatch: expected {original_length}, got {len(record.seq)}"
