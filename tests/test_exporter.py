"""Property-based tests for exporter module."""
import json
from hypothesis import given, strategies as st, settings
from io import StringIO
from Bio import Phylo
from modules.exporter import OutputExporter
from modules.tree_builder import TreeBuilder
from modules.aligner import SequenceAligner
from modules.data_models import Sequence


# Custom strategies
@st.composite
def tree_from_sequences(draw, min_seqs=3, max_seqs=5):
    """Generate a tree from aligned sequences."""
    num_seqs = draw(st.integers(min_value=min_seqs, max_value=max_seqs))
    base_length = draw(st.integers(min_value=20, max_value=40))
    
    # Generate a base sequence
    base_seq = ''.join(draw(st.lists(
        st.sampled_from('ACGT'),
        min_size=base_length,
        max_size=base_length
    )))
    
    sequences = [Sequence(id="Seq_0", sequence=base_seq)]
    
    # Generate related sequences
    for i in range(1, num_seqs):
        variant = list(base_seq)
        num_mutations = draw(st.integers(min_value=1, max_value=min(5, len(variant) // 4)))
        for _ in range(num_mutations):
            pos = draw(st.integers(min_value=0, max_value=len(variant) - 1))
            variant[pos] = draw(st.sampled_from('ACGT'))
        sequences.append(Sequence(id=f"Seq_{i}", sequence=''.join(variant)))
    
    # Align and build tree
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="nj")
    
    return tree, sequences


# Feature: phylogenetics-app, Property 9: Newick round-trip consistency
# Validates: Requirements 4.2
@given(tree_from_sequences(min_seqs=3, max_seqs=5))
@settings(deadline=None)
def test_newick_round_trip_preserves_topology(tree_data):
    """
    Property: For any tree, exporting to Newick and parsing back should
    preserve the tree topology and terminal names.
    """
    tree, sequences = tree_data
    exporter = OutputExporter()
    
    # Export to Newick
    newick_str = exporter.export_newick(tree)
    
    # Parse back
    input_stream = StringIO(newick_str)
    parsed_tree = Phylo.read(input_stream, 'newick')
    
    # Check same number of terminals
    original_terminals = tree.get_terminals()
    parsed_terminals = parsed_tree.get_terminals()
    
    assert len(parsed_terminals) == len(original_terminals), \
        f"Terminal count mismatch: {len(parsed_terminals)} vs {len(original_terminals)}"
    
    # Check terminal names are preserved
    original_names = sorted([t.name for t in original_terminals])
    parsed_names = sorted([t.name for t in parsed_terminals])
    
    assert parsed_names == original_names, \
        f"Terminal names not preserved: {parsed_names} vs {original_names}"


@given(tree_from_sequences(min_seqs=3, max_seqs=5))
@settings(deadline=None)
def test_newick_round_trip_preserves_branch_lengths(tree_data):
    """
    Property: Newick round-trip should preserve branch lengths within
    numerical precision.
    """
    tree, sequences = tree_data
    exporter = OutputExporter()
    
    # Export to Newick
    newick_str = exporter.export_newick(tree)
    
    # Parse back
    input_stream = StringIO(newick_str)
    parsed_tree = Phylo.read(input_stream, 'newick')
    
    # Get terminal names and their branch lengths
    original_terminals = {t.name: t.branch_length for t in tree.get_terminals()}
    parsed_terminals = {t.name: t.branch_length for t in parsed_tree.get_terminals()}
    
    # Check branch lengths are approximately equal
    # Newick format has limited precision (typically 5-6 decimal places)
    epsilon = 1e-4
    for name in original_terminals:
        if name in parsed_terminals:
            orig_len = original_terminals[name]
            parsed_len = parsed_terminals[name]
            
            if orig_len is not None and parsed_len is not None:
                assert abs(orig_len - parsed_len) < epsilon, \
                    f"Branch length mismatch for {name}: {orig_len} vs {parsed_len}"


@given(tree_from_sequences(min_seqs=3, max_seqs=4))
@settings(deadline=None)
def test_newick_export_is_valid_format(tree_data):
    """
    Property: Exported Newick string should be parseable by Biopython.
    """
    tree, sequences = tree_data
    exporter = OutputExporter()
    
    # Export to Newick
    newick_str = exporter.export_newick(tree)
    
    # Should be non-empty
    assert len(newick_str) > 0, "Newick export is empty"
    
    # Should be parseable
    try:
        input_stream = StringIO(newick_str)
        parsed_tree = Phylo.read(input_stream, 'newick')
        assert parsed_tree is not None
    except Exception as e:
        assert False, f"Failed to parse Newick string: {e}"



# Feature: phylogenetics-app, Property 13: Closest relatives completeness
# Validates: Requirements 7.1, 7.2, 7.5
@given(tree_from_sequences(min_seqs=3, max_seqs=6))
@settings(deadline=None)
def test_closest_relatives_has_all_sequences(tree_data):
    """
    Property: For any tree with N sequences, the closest relatives JSON
    should contain exactly N entries.
    """
    tree, sequences = tree_data
    exporter = OutputExporter()
    
    # Calculate closest relatives
    relatives_dict = exporter.calculate_closest_relatives(tree)
    
    # Check count matches
    assert len(relatives_dict) == len(sequences), \
        f"Expected {len(sequences)} entries, got {len(relatives_dict)}"
    
    # Check all sequence IDs are present
    seq_ids = [seq.id for seq in sequences]
    for seq_id in seq_ids:
        assert seq_id in relatives_dict, \
            f"Sequence ID '{seq_id}' not found in relatives dict"


@given(tree_from_sequences(min_seqs=3, max_seqs=6))
@settings(deadline=None)
def test_closest_relatives_have_non_negative_distances(tree_data):
    """
    Property: All distance values in closest relatives should be non-negative.
    """
    tree, sequences = tree_data
    exporter = OutputExporter()
    
    # Calculate closest relatives
    relatives_dict = exporter.calculate_closest_relatives(tree)
    
    # Check all distances are non-negative (within numerical precision)
    epsilon = 1e-10
    for seq_id, rel_info in relatives_dict.items():
        assert rel_info.distance >= -epsilon, \
            f"Significantly negative distance for {seq_id}: {rel_info.distance}"
        assert isinstance(rel_info.distance, (int, float)), \
            f"Distance is not a number for {seq_id}"


@given(tree_from_sequences(min_seqs=3, max_seqs=5))
@settings(deadline=None)
def test_closest_relatives_json_is_valid(tree_data):
    """
    Property: Exported JSON should be valid and parseable.
    """
    tree, sequences = tree_data
    exporter = OutputExporter()
    
    # Export to JSON
    json_str = exporter.export_closest_relatives(tree)
    
    # Should be parseable
    try:
        parsed = json.loads(json_str)
        assert isinstance(parsed, dict)
    except json.JSONDecodeError as e:
        assert False, f"Invalid JSON: {e}"
    
    # Check structure
    for seq_id, entry in parsed.items():
        assert 'closest_relative' in entry, \
            f"Missing 'closest_relative' for {seq_id}"
        assert 'distance' in entry, \
            f"Missing 'distance' for {seq_id}"
        assert 'explanation' in entry, \
            f"Missing 'explanation' for {seq_id}"


@given(tree_from_sequences(min_seqs=3, max_seqs=5))
@settings(deadline=None)
def test_closest_relatives_have_explanations(tree_data):
    """
    Property: Each entry should have a non-empty explanation string.
    """
    tree, sequences = tree_data
    exporter = OutputExporter()
    
    # Calculate closest relatives
    relatives_dict = exporter.calculate_closest_relatives(tree)
    
    # Check all have explanations
    for seq_id, rel_info in relatives_dict.items():
        assert isinstance(rel_info.explanation, str), \
            f"Explanation is not a string for {seq_id}"
        assert len(rel_info.explanation) > 0, \
            f"Empty explanation for {seq_id}"
        # Explanation should mention the sequence ID
        assert seq_id in rel_info.explanation, \
            f"Explanation doesn't mention {seq_id}"


@given(tree_from_sequences(min_seqs=3, max_seqs=5))
@settings(deadline=None)
def test_closest_relative_is_different_sequence(tree_data):
    """
    Property: A sequence's closest relative should be a different sequence.
    """
    tree, sequences = tree_data
    exporter = OutputExporter()
    
    # Calculate closest relatives
    relatives_dict = exporter.calculate_closest_relatives(tree)
    
    # Check closest relative is not self
    for seq_id, rel_info in relatives_dict.items():
        if rel_info.closest_relative != "None":
            assert rel_info.closest_relative != seq_id, \
                f"Sequence {seq_id} is its own closest relative"
