"""Integration tests for the complete phylogenetics workflow."""
from modules.parser import SequenceParser
from modules.validator import SequenceValidator
from modules.aligner import SequenceAligner
from modules.tree_builder import TreeBuilder
from modules.visualizer import TreeVisualizer
from modules.exporter import OutputExporter


def test_complete_workflow_with_fasta():
    """Test the complete workflow from FASTA input to all outputs."""
    # Sample FASTA content
    fasta_content = """>Seq_1
ACGTACGTACGTACGT
>Seq_2
ACGTACGTACGTACGG
>Seq_3
ACGTACGTACGTACGA
"""
    
    # Step 1: Parse
    parser = SequenceParser()
    sequences = parser.parse_file(fasta_content)
    assert len(sequences) == 3
    
    # Step 2: Validate
    validator = SequenceValidator(min_length=10)
    validation_result = validator.validate_sequences(sequences)
    assert len(validation_result.valid_sequences) == 3
    
    # Step 3: Align
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(validation_result.valid_sequences)
    assert len(alignment) == 3
    
    # Step 4: Build tree
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="nj")
    assert len(tree.get_terminals()) == 3
    
    # Step 5: Visualize
    visualizer = TreeVisualizer()
    fig = visualizer.draw_tree(tree, layout='rectangular', theme='default')
    assert fig is not None
    
    # Step 6: Export
    exporter = OutputExporter()
    
    # Export Newick
    newick_str = exporter.export_newick(tree)
    assert len(newick_str) > 0
    assert 'Seq_1' in newick_str
    
    # Export aligned FASTA
    aligned_fasta = exporter.export_aligned_fasta(alignment)
    assert len(aligned_fasta) > 0
    assert '>Seq_1' in aligned_fasta
    
    # Export closest relatives
    relatives_json = exporter.export_closest_relatives(tree)
    assert len(relatives_json) > 0
    assert 'Seq_1' in relatives_json
    
    # Clean up
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_complete_workflow_with_txt():
    """Test the complete workflow with TXT input."""
    # Sample TXT content with named sequences
    txt_content = """Seq_1:ACGTACGTACGTACGT
Seq_2:ACGTACGTACGTACGG
Seq_3:ACGTACGTACGTACGA
"""
    
    # Parse
    parser = SequenceParser()
    sequences = parser.parse_file(txt_content)
    assert len(sequences) == 3
    
    # Validate
    validator = SequenceValidator(min_length=10)
    validation_result = validator.validate_sequences(sequences)
    assert len(validation_result.valid_sequences) == 3
    
    # Align
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(validation_result.valid_sequences)
    assert len(alignment) == 3
    
    # Build tree
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="upgma")
    assert len(tree.get_terminals()) == 3
    
    # Export
    exporter = OutputExporter()
    newick_str = exporter.export_newick(tree)
    assert len(newick_str) > 0


def test_workflow_with_different_settings():
    """Test workflow with different visualization settings."""
    fasta_content = """>A
AAAAAAAAAAAAAAAA
>B
AAAAAAAAAAAAAAAC
>C
AAAAAAAAAAAAAAAG
>D
AAAAAAAAAAAAAAAA
"""
    
    # Parse and process
    parser = SequenceParser()
    sequences = parser.parse_file(fasta_content)
    
    validator = SequenceValidator()
    validation_result = validator.validate_sequences(sequences)
    
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(validation_result.valid_sequences)
    
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="nj")
    
    # Test different layouts and themes
    visualizer = TreeVisualizer()
    
    layouts = ['rectangular', 'circular']
    themes = ['default', 'dark', 'colorful', 'minimal']
    
    for layout in layouts:
        for theme in themes:
            fig = visualizer.draw_tree(tree, layout=layout, theme=theme)
            assert fig is not None
            
            # Clean up
            import matplotlib.pyplot as plt
            plt.close(fig)


def test_workflow_with_sequence_cleaning():
    """Test workflow with sequences that need cleaning."""
    # Sequences with invalid characters
    fasta_content = """>Seq1
ACGT 123 ACGT
>Seq2
ACGT-ACGT
>Seq3
ACGT*ACGT
"""
    
    parser = SequenceParser()
    sequences = parser.parse_file(fasta_content)
    
    validator = SequenceValidator(min_length=5)
    validation_result = validator.validate_sequences(sequences)
    
    # Should have cleaned the sequences
    assert len(validation_result.valid_sequences) == 3
    assert validation_result.cleaning_summary['removed_spaces'] > 0 or \
           validation_result.cleaning_summary['removed_numbers'] > 0 or \
           validation_result.cleaning_summary['removed_special_chars'] > 0
    
    # Continue workflow
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(validation_result.valid_sequences)
    
    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method="nj")
    
    assert len(tree.get_terminals()) == 3


def test_error_handling_insufficient_sequences():
    """Test that workflow properly handles insufficient sequences."""
    fasta_content = """>Seq1
ACGTACGTACGT
>Seq2
ACGTACGGACGT
"""
    
    parser = SequenceParser()
    sequences = parser.parse_file(fasta_content)
    
    validator = SequenceValidator()
    validation_result = validator.validate_sequences(sequences)
    
    # Should have 2 valid sequences
    assert len(validation_result.valid_sequences) == 2
    
    aligner = SequenceAligner()
    alignment = aligner.align_sequences(validation_result.valid_sequences)
    
    # Should raise error with < 3 sequences
    builder = TreeBuilder()
    try:
        tree = builder.build_tree(alignment, method="nj")
        assert False, "Should have raised TreeBuildError"
    except Exception as e:
        assert "at least 3 sequences" in str(e).lower()
