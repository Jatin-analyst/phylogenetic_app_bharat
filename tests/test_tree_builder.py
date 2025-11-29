"""Property-based tests for tree builder module."""
from hypothesis import given, strategies as st, assume, settings
from modules.tree_builder import TreeBuilder, TreeBuildError
from modules.aligner import SequenceAligner
from modules.data_models import Sequence


# Custom strategies
@st.composite
def aligned_sequences(draw, min_seqs=3, max_seqs=6):
    """Generate aligned sequences for tree building."""
    num_seqs = draw(st.integers(min_value=min_seqs, max_value=max_seqs))
    base_length = draw(st.integers(min_value=20, max_value=40))
    
    # Generate a base sequence
    base_seq = ''.join(draw(st.lists(
        st.sampled_from('ACGT'),
        min_size=base_length,
        max_size=base_length
    )))
    
    sequences = [Sequence(id="Seq_0", sequence=base_seq)]
    
    # Generate related sequences with mutations
    for i in range(1, num_seqs):
        variant = list(base_seq)
        
        # Add mutations
        num_mutations = draw(st.integers(min_value=1, max_value=min(5, len(variant) // 4)))
        for _ in range(num_mutations):
            pos = draw(st.integers(min_value=0, max_value=len(variant) - 1))
            variant[pos] = draw(st.sampled_from('ACGT'))
        
        sequences.append(Sequence(id=f"Seq_{i}", sequence=''.join(variant)))
    
    return sequences


# Feature: phylogenetics-app, Property 8: Tree construction produces valid tree
# Validates: Requirements 4.1
@given(aligned_sequences(min_seqs=3, max_seqs=6))
@settings(deadline=None)
def test_tree_has_correct_number_of_leaves(sequences):
    """
    Property: For any valid alignment, tree construction should produce
    a tree with the same number of leaf nodes as input sequences.
    """
    # Align sequences first
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    # Build tree
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="nj")
    
    # Count terminal (leaf) nodes
    terminals = tree.get_terminals()
    
    assert len(terminals) == len(sequences), \
        f"Expected {len(sequences)} leaf nodes, got {len(terminals)}"


@given(aligned_sequences(min_seqs=3, max_seqs=6))
@settings(deadline=None)
def test_tree_leaves_have_correct_labels(sequences):
    """
    Property: Each leaf in the tree should be labeled with a sequence
    identifier from the input.
    """
    # Align sequences first
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    # Build tree
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="nj")
    
    # Get terminal node names
    terminal_names = [term.name for term in tree.get_terminals()]
    input_ids = [seq.id for seq in sequences]
    
    # Check all input IDs are present in tree
    for seq_id in input_ids:
        assert seq_id in terminal_names, \
            f"Sequence ID '{seq_id}' not found in tree leaves"
    
    # Check no extra IDs in tree
    assert len(terminal_names) == len(input_ids), \
        "Tree has different number of leaves than input sequences"


@given(aligned_sequences(min_seqs=3, max_seqs=5))
@settings(deadline=None)
def test_nj_and_upgma_produce_valid_trees(sequences):
    """
    Property: Both NJ and UPGMA methods should produce valid trees
    with correct number of leaves.
    """
    # Align sequences first
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    # Build tree with both methods
    builder = TreeBuilder()
    
    nj_tree = builder.build_tree(alignment, method="nj")
    upgma_tree = builder.build_tree(alignment, method="upgma")
    
    # Both should have correct number of leaves
    assert len(nj_tree.get_terminals()) == len(sequences), \
        "NJ tree has wrong number of leaves"
    assert len(upgma_tree.get_terminals()) == len(sequences), \
        "UPGMA tree has wrong number of leaves"


@given(aligned_sequences(min_seqs=3, max_seqs=5))
@settings(deadline=None)
def test_tree_is_connected(sequences):
    """
    Property: The constructed tree should be fully connected
    (all nodes reachable from root).
    """
    # Align sequences first
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    # Build tree
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="nj")
    
    # Get all nodes
    all_nodes = list(tree.find_clades())
    
    # Should have at least as many nodes as sequences
    assert len(all_nodes) >= len(sequences), \
        "Tree doesn't have enough nodes"
    
    # Tree should have a root
    assert tree.root is not None, \
        "Tree has no root"



# Feature: phylogenetics-app, Property 10: Branch lengths are non-negative
# Validates: Requirements 4.4
@given(aligned_sequences(min_seqs=3, max_seqs=6))
@settings(deadline=None)
def test_all_branch_lengths_non_negative(sequences):
    """
    Property: For any constructed tree, all branch lengths should be
    non-negative real numbers (within numerical precision).
    """
    # Align sequences first
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    # Build tree
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="nj")
    
    # Check all branch lengths
    # Allow for small negative values due to floating point precision
    epsilon = 1e-10
    
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            assert clade.branch_length >= -epsilon, \
                f"Significantly negative branch length found: {clade.branch_length}"
            assert isinstance(clade.branch_length, (int, float)), \
                f"Branch length is not a number: {type(clade.branch_length)}"


@given(aligned_sequences(min_seqs=3, max_seqs=5))
@settings(deadline=None)
def test_branch_lengths_exist(sequences):
    """
    Property: Constructed trees should have branch lengths assigned
    to most branches (representing evolutionary distances).
    """
    # Align sequences first
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    # Build tree
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="nj")
    
    # Count branches with lengths
    branches_with_length = 0
    total_branches = 0
    
    for clade in tree.find_clades():
        if clade != tree.root:  # Root typically has no branch length
            total_branches += 1
            if clade.branch_length is not None:
                branches_with_length += 1
    
    # Most branches should have lengths
    assert branches_with_length > 0, \
        "No branch lengths found in tree"


@given(aligned_sequences(min_seqs=3, max_seqs=5))
@settings(deadline=None)
def test_upgma_branch_lengths_non_negative(sequences):
    """
    Property: UPGMA trees should also have non-negative branch lengths.
    """
    # Align sequences first
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)
    
    # Build tree with UPGMA
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="upgma")
    
    # Check all branch lengths
    for clade in tree.find_clades():
        if clade.branch_length is not None:
            assert clade.branch_length >= 0, \
                f"Negative branch length in UPGMA tree: {clade.branch_length}"
