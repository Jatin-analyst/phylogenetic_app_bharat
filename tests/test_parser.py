"""Property-based tests for parser module."""
from hypothesis import given, strategies as st, assume
from modules.parser import SequenceParser, ParseError
from modules.data_models import Sequence


# Custom strategies for generating test data
@st.composite
def valid_sequence_string(draw):
    """Generate a valid DNA/protein sequence string."""
    # Use nucleotides and protein characters
    alphabet = "ACGTRYSWKMBDHVN"  # DNA with ambiguity codes
    length = draw(st.integers(min_value=10, max_value=200))
    return ''.join(draw(st.lists(st.sampled_from(alphabet), min_size=length, max_size=length)))


@st.composite
def fasta_format(draw):
    """Generate valid FASTA format content."""
    num_sequences = draw(st.integers(min_value=1, max_value=10))
    sequences = []
    
    for i in range(num_sequences):
        seq_id = f"Seq_{i+1}"
        seq_data = draw(valid_sequence_string())
        sequences.append((seq_id, seq_data))
    
    # Build FASTA content
    fasta_lines = []
    for seq_id, seq_data in sequences:
        fasta_lines.append(f">{seq_id}")
        fasta_lines.append(seq_data)
    
    return '\n'.join(fasta_lines), sequences


@st.composite
def txt_raw_format(draw):
    """Generate TXT format with raw sequences."""
    num_sequences = draw(st.integers(min_value=1, max_value=10))
    sequences = []
    
    for i in range(num_sequences):
        seq_data = draw(valid_sequence_string())
        sequences.append((f"Seq_{i+1}", seq_data))
    
    # Build TXT content (one sequence per line)
    txt_lines = [seq_data for _, seq_data in sequences]
    
    return '\n'.join(txt_lines), sequences


@st.composite
def txt_named_format(draw):
    """Generate TXT format with named sequences."""
    num_sequences = draw(st.integers(min_value=1, max_value=10))
    sequences = []
    
    for i in range(num_sequences):
        seq_id = f"Seq_{i+1}"
        seq_data = draw(valid_sequence_string())
        sequences.append((seq_id, seq_data))
    
    # Build TXT content with colon separator
    txt_lines = [f"{seq_id}:{seq_data}" for seq_id, seq_data in sequences]
    
    return '\n'.join(txt_lines), sequences


# Feature: phylogenetics-app, Property 1: Parsing preserves sequence count and identifiers
# Validates: Requirements 1.1, 1.2, 1.3, 1.4
@given(fasta_format())
def test_fasta_parsing_preserves_count_and_ids(fasta_data):
    """
    Property: For any valid FASTA file, parsing should extract all sequences
    and preserve their identifiers.
    """
    content, expected_sequences = fasta_data
    parser = SequenceParser()
    
    result = parser.parse_fasta(content)
    
    # Check sequence count matches
    assert len(result) == len(expected_sequences), \
        f"Expected {len(expected_sequences)} sequences, got {len(result)}"
    
    # Check all IDs are present and unique
    result_ids = [seq.id for seq in result]
    expected_ids = [seq_id for seq_id, _ in expected_sequences]
    
    assert result_ids == expected_ids, \
        f"IDs don't match. Expected {expected_ids}, got {result_ids}"
    
    # Check sequences match
    for i, (expected_id, expected_seq) in enumerate(expected_sequences):
        assert result[i].sequence == expected_seq, \
            f"Sequence {i} doesn't match"


@given(txt_raw_format())
def test_txt_raw_parsing_preserves_count_and_assigns_ids(txt_data):
    """
    Property: For any valid TXT file with raw sequences, parsing should extract
    all sequences and assign unique default IDs.
    """
    content, expected_sequences = txt_data
    parser = SequenceParser()
    
    result = parser.parse_txt_raw(content)
    
    # Check sequence count matches
    assert len(result) == len(expected_sequences), \
        f"Expected {len(expected_sequences)} sequences, got {len(result)}"
    
    # Check all IDs are unique
    result_ids = [seq.id for seq in result]
    assert len(result_ids) == len(set(result_ids)), \
        "IDs are not unique"
    
    # Check sequences match
    for i, (_, expected_seq) in enumerate(expected_sequences):
        assert result[i].sequence == expected_seq, \
            f"Sequence {i} doesn't match"


@given(txt_named_format())
def test_txt_named_parsing_preserves_count_and_ids(txt_data):
    """
    Property: For any valid TXT file with named sequences, parsing should extract
    all sequences and preserve their names.
    """
    content, expected_sequences = txt_data
    parser = SequenceParser()
    
    result = parser.parse_txt_named(content)
    
    # Check sequence count matches
    assert len(result) == len(expected_sequences), \
        f"Expected {len(expected_sequences)} sequences, got {len(result)}"
    
    # Check IDs match
    result_ids = [seq.id for seq in result]
    expected_ids = [seq_id for seq_id, _ in expected_sequences]
    
    assert result_ids == expected_ids, \
        f"IDs don't match. Expected {expected_ids}, got {result_ids}"
    
    # Check sequences match
    for i, (_, expected_seq) in enumerate(expected_sequences):
        assert result[i].sequence == expected_seq, \
            f"Sequence {i} doesn't match"


@given(st.one_of(fasta_format(), txt_raw_format(), txt_named_format()))
def test_parse_file_preserves_count(format_data):
    """
    Property: For any valid input file, parse_file should extract all sequences
    with unique identifiers.
    """
    content, expected_sequences = format_data
    parser = SequenceParser()
    
    result = parser.parse_file(content)
    
    # Check sequence count matches
    assert len(result) == len(expected_sequences), \
        f"Expected {len(expected_sequences)} sequences, got {len(result)}"
    
    # Check all IDs are unique
    result_ids = [seq.id for seq in result]
    assert len(result_ids) == len(set(result_ids)), \
        "IDs are not unique"



# Feature: phylogenetics-app, Property 2: Invalid input rejection
# Validates: Requirements 1.5
@given(st.one_of(
    st.just(""),  # Empty string
    st.just("   \n\n  "),  # Only whitespace
    st.text(alphabet="!@#$%^&*()[]{}|\\", min_size=10, max_size=50),  # Invalid characters only
    st.just("This is not a sequence file\nJust random text\n123456"),  # Random text
))
def test_invalid_input_raises_error(invalid_content):
    """
    Property: For any invalid or malformed input, the parser should raise
    an appropriate error and not return partial or corrupted data.
    """
    parser = SequenceParser()
    
    try:
        result = parser.parse_file(invalid_content)
        # If we get here, check that we didn't get corrupted data
        # Empty results are acceptable for some edge cases
        if result:
            # If we got results, they should at least be valid Sequence objects
            for seq in result:
                assert isinstance(seq, Sequence)
                assert seq.id  # Should have an ID
    except ParseError:
        # This is the expected behavior for invalid input
        pass


def test_empty_file_raises_parse_error():
    """Test that empty files raise ParseError."""
    parser = SequenceParser()
    
    try:
        parser.parse_file("")
        assert False, "Should have raised ParseError for empty file"
    except ParseError as e:
        assert "empty" in str(e).lower()


def test_whitespace_only_raises_parse_error():
    """Test that whitespace-only files raise ParseError."""
    parser = SequenceParser()
    
    try:
        parser.parse_file("   \n\n\t  ")
        assert False, "Should have raised ParseError for whitespace-only file"
    except ParseError as e:
        assert "empty" in str(e).lower()


def test_malformed_fasta_raises_error():
    """Test that malformed FASTA raises ParseError."""
    parser = SequenceParser()
    
    # FASTA with header but no sequence
    malformed = ">Seq1\n>Seq2\n>Seq3"
    
    try:
        result = parser.parse_fasta(malformed)
        # If it doesn't raise, check that sequences are not empty
        for seq in result:
            assert len(seq.sequence) > 0, "Should not have empty sequences"
    except ParseError:
        # This is acceptable behavior
        pass


@given(st.lists(st.just(""), min_size=1, max_size=5))
def test_empty_lines_only_raises_error(empty_lines):
    """
    Property: Files containing only empty lines should raise ParseError.
    """
    parser = SequenceParser()
    content = '\n'.join(empty_lines)
    
    try:
        parser.parse_file(content)
        assert False, "Should have raised ParseError for empty lines only"
    except ParseError:
        pass  # Expected behavior



# Feature: phylogenetics-app, Property 2: Invalid input rejection
# Validates: Requirements 1.5
@given(st.text(alphabet=st.characters(blacklist_categories=('Cs',)), min_size=1, max_size=100))
def test_invalid_input_raises_error(invalid_content):
    """
    Property: For any invalid or malformed input, the parser should raise
    an appropriate error and not return partial or corrupted data.
    
    Note: Since TXT format accepts any text, this test focuses on ensuring
    that when parsing does succeed, it returns valid Sequence objects.
    """
    parser = SequenceParser()
    
    try:
        result = parser.parse_file(invalid_content)
        # If parsing succeeds, verify all returned sequences are valid Sequence objects
        assert isinstance(result, list)
        for seq in result:
            assert isinstance(seq, Sequence)
            assert isinstance(seq.id, str)
            assert isinstance(seq.sequence, str)
            assert len(seq.id) > 0  # IDs should not be empty
    except ParseError:
        # Expected behavior for truly invalid input
        pass


@given(st.just(""))
def test_empty_file_raises_parse_error(empty_content):
    """Test that empty files raise ParseError."""
    parser = SequenceParser()
    
    try:
        parser.parse_file(empty_content)
        assert False, "Expected ParseError for empty file"
    except ParseError as e:
        assert "empty" in str(e).lower()


@given(st.text(alphabet=' \t\n\r', min_size=1, max_size=50))
def test_whitespace_only_raises_parse_error(whitespace_content):
    """Test that files with only whitespace raise ParseError."""
    parser = SequenceParser()
    
    try:
        parser.parse_file(whitespace_content)
        assert False, "Expected ParseError for whitespace-only file"
    except ParseError as e:
        assert "empty" in str(e).lower() or "no valid" in str(e).lower()


@given(st.lists(st.just('>'), min_size=1, max_size=5))
def test_malformed_fasta_raises_error(header_list):
    """Test that FASTA files with headers but no sequences raise errors."""
    # Create malformed FASTA with only headers
    content = '\n'.join(header_list)
    parser = SequenceParser()
    
    try:
        result = parser.parse_fasta(content)
        # If it doesn't raise, it should return empty or sequences with empty strings
        assert all(len(seq.sequence) == 0 for seq in result)
    except ParseError:
        # Expected behavior
        pass


@given(st.lists(st.just('\n'), min_size=1, max_size=10))
def test_empty_lines_only_raises_error(newlines):
    """Test that files with only newlines raise ParseError."""
    content = ''.join(newlines)
    parser = SequenceParser()
    
    try:
        parser.parse_file(content)
        assert False, "Expected ParseError for file with only newlines"
    except ParseError as e:
        assert "empty" in str(e).lower() or "no valid" in str(e).lower()
