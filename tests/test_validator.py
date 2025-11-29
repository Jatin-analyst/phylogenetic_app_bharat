"""Property-based tests for validator module."""
from hypothesis import given, strategies as st, assume
from modules.validator import SequenceValidator
from modules.data_models import Sequence


# Custom strategies
@st.composite
def sequence_with_invalid_chars(draw, seq_type='nucleotide'):
    """Generate a sequence with mix of valid and invalid characters."""
    if seq_type == 'nucleotide':
        valid_chars = 'ACGTRYSWKMBDHVN'
    else:
        valid_chars = 'ACDEFGHIKLMNPQRSTVWY'
    
    # Generate valid sequence
    valid_length = draw(st.integers(min_value=5, max_value=50))
    valid_seq = ''.join(draw(st.lists(
        st.sampled_from(valid_chars),
        min_size=valid_length,
        max_size=valid_length
    )))
    
    # Add invalid characters (spaces, numbers, special chars)
    invalid_chars = [' ', '\t', '\n', '0', '1', '2', '-', '*', '.']
    num_invalid = draw(st.integers(min_value=0, max_value=10))
    
    # Insert invalid characters at random positions
    seq_list = list(valid_seq)
    for _ in range(num_invalid):
        pos = draw(st.integers(min_value=0, max_value=len(seq_list)))
        char = draw(st.sampled_from(invalid_chars))
        seq_list.insert(pos, char)
    
    return ''.join(seq_list), valid_seq


# Feature: phylogenetics-app, Property 3: Sequence cleaning preserves valid characters
# Validates: Requirements 2.1
@given(sequence_with_invalid_chars('nucleotide'))
def test_nucleotide_cleaning_preserves_valid_chars(seq_data):
    """
    Property: For any nucleotide sequence with invalid characters,
    cleaning should preserve all valid characters in their original order.
    """
    dirty_seq, expected_clean = seq_data
    validator = SequenceValidator()
    
    cleaned = validator.clean_sequence(dirty_seq, 'nucleotide')
    
    # Check that cleaned sequence matches expected (all valid chars preserved)
    assert cleaned == expected_clean.upper(), \
        f"Expected '{expected_clean.upper()}', got '{cleaned}'"
    
    # Verify no invalid characters remain
    valid_chars = set('ACGTRYSWKMBDHVN')
    assert all(c in valid_chars for c in cleaned), \
        "Cleaned sequence contains invalid characters"


@given(sequence_with_invalid_chars('protein'))
def test_protein_cleaning_preserves_valid_chars(seq_data):
    """
    Property: For any protein sequence with invalid characters,
    cleaning should preserve all valid characters in their original order.
    """
    dirty_seq, expected_clean = seq_data
    validator = SequenceValidator()
    
    cleaned = validator.clean_sequence(dirty_seq, 'protein')
    
    # Check that cleaned sequence matches expected
    assert cleaned == expected_clean.upper(), \
        f"Expected '{expected_clean.upper()}', got '{cleaned}'"
    
    # Verify no invalid characters remain
    valid_chars = set('ACDEFGHIKLMNPQRSTVWY')
    assert all(c in valid_chars for c in cleaned), \
        "Cleaned sequence contains invalid characters"


@given(st.text(alphabet='ACGT', min_size=10, max_size=100))
def test_cleaning_pure_sequence_unchanged(pure_seq):
    """
    Property: Cleaning a sequence with only valid characters should
    return the same sequence (uppercased).
    """
    validator = SequenceValidator()
    
    cleaned = validator.clean_sequence(pure_seq, 'nucleotide')
    
    assert cleaned == pure_seq.upper(), \
        "Pure sequence should remain unchanged (except case)"



# Feature: phylogenetics-app, Property 4: Sequence type detection consistency
# Validates: Requirements 2.2
@given(st.lists(st.text(alphabet='ACGT', min_size=20, max_size=100), min_size=1, max_size=10))
def test_nucleotide_detection_consistency(sequences):
    """
    Property: For any set of nucleotide sequences, the validator should
    consistently identify them as nucleotide type.
    """
    validator = SequenceValidator()
    
    seq_type = validator.detect_sequence_type(sequences)
    
    assert seq_type == 'nucleotide', \
        f"Expected 'nucleotide', got '{seq_type}'"


@given(st.lists(st.text(alphabet='ACDEFGHIKLMNPQRSTVWY', min_size=20, max_size=100), min_size=1, max_size=10))
def test_protein_detection_consistency(sequences):
    """
    Property: For any set of protein sequences (with protein-specific amino acids),
    the validator should consistently identify them as protein type.
    """
    # Ensure at least some sequences have protein-specific characters
    assume(any(any(c in 'EFIKLMPQVWY' for c in seq) for seq in sequences))
    
    validator = SequenceValidator()
    
    seq_type = validator.detect_sequence_type(sequences)
    
    assert seq_type == 'protein', \
        f"Expected 'protein', got '{seq_type}'"


@given(st.lists(st.text(alphabet='ACGT', min_size=20, max_size=100), min_size=1, max_size=5))
def test_same_type_sequences_detected_consistently(sequences):
    """
    Property: For any set of sequences of the same type, detection should
    return consistent results.
    """
    validator = SequenceValidator()
    
    # Detect type multiple times
    type1 = validator.detect_sequence_type(sequences)
    type2 = validator.detect_sequence_type(sequences)
    
    assert type1 == type2, \
        "Type detection should be consistent for same input"



# Feature: phylogenetics-app, Property 5: Short sequence filtering
# Validates: Requirements 2.3
@given(
    st.lists(st.text(alphabet='ACGT', min_size=1, max_size=50), min_size=1, max_size=20),
    st.integers(min_value=5, max_value=30)
)
def test_short_sequence_filtering(sequences, min_length):
    """
    Property: For any set of sequences and minimum length threshold,
    all sequences shorter than the threshold should be excluded,
    and all sequences meeting the threshold should be included.
    """
    validator = SequenceValidator(min_length=min_length)
    
    # Create Sequence objects
    seq_objects = [Sequence(id=f"Seq_{i}", sequence=seq) for i, seq in enumerate(sequences)]
    
    # Validate
    result = validator.validate_sequences(seq_objects)
    
    # Check that all valid sequences meet minimum length
    for seq in result.valid_sequences:
        assert len(seq.sequence) >= min_length, \
            f"Valid sequence '{seq.id}' has length {len(seq.sequence)}, expected >= {min_length}"
    
    # Check that all invalid sequences are too short
    for seq in result.invalid_sequences:
        assert len(seq.sequence) < min_length, \
            f"Invalid sequence '{seq.id}' has length {len(seq.sequence)}, expected < {min_length}"
    
    # Check that total count is preserved
    assert len(result.valid_sequences) + len(result.invalid_sequences) == len(sequences), \
        "Total sequence count should be preserved"


@given(st.lists(st.text(alphabet='ACGT', min_size=15, max_size=50), min_size=1, max_size=10))
def test_all_long_sequences_pass_validation(sequences):
    """
    Property: For sequences all longer than minimum length,
    all should pass validation.
    """
    min_length = 10
    validator = SequenceValidator(min_length=min_length)
    
    # Create Sequence objects
    seq_objects = [Sequence(id=f"Seq_{i}", sequence=seq) for i, seq in enumerate(sequences)]
    
    # Validate
    result = validator.validate_sequences(seq_objects)
    
    # All should be valid
    assert len(result.valid_sequences) == len(sequences), \
        "All sequences should pass validation"
    assert len(result.invalid_sequences) == 0, \
        "No sequences should be invalid"


@given(st.lists(st.text(alphabet='ACGT', min_size=1, max_size=5), min_size=1, max_size=10))
def test_all_short_sequences_fail_validation(sequences):
    """
    Property: For sequences all shorter than minimum length,
    all should fail validation.
    """
    min_length = 10
    validator = SequenceValidator(min_length=min_length)
    
    # Create Sequence objects
    seq_objects = [Sequence(id=f"Seq_{i}", sequence=seq) for i, seq in enumerate(sequences)]
    
    # Validate
    result = validator.validate_sequences(seq_objects)
    
    # All should be invalid
    assert len(result.valid_sequences) == 0, \
        "No sequences should pass validation"
    assert len(result.invalid_sequences) == len(sequences), \
        "All sequences should be invalid"
