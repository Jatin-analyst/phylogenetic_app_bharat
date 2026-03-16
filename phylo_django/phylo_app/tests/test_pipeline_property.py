"""Pipeline property tests using Hypothesis.
# Feature: django-deployment
"""
import pytest
from io import StringIO
from hypothesis import given, settings as h_settings, HealthCheck
from hypothesis import strategies as st
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from modules.parser import SequenceParser
from modules.aligner import SequenceAligner
from modules.tree_builder import TreeBuilder
from modules.exporter import OutputExporter


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

DNA_CHARS = st.sampled_from('ACGT')

@st.composite
def fasta_sequence_strategy(draw):
    """Generate a single valid FASTA sequence (id + nucleotide string ≥ 10 bp)."""
    seq_id = draw(st.text(alphabet='abcdefghijklmnopqrstuvwxyz0123456789', min_size=1, max_size=10))
    length = draw(st.integers(min_value=30, max_value=80))
    seq = ''.join(draw(st.lists(DNA_CHARS, min_size=length, max_size=length)))
    return f'>{seq_id}\n{seq}\n'


@st.composite
def valid_fasta_content_strategy(draw):
    """Generate a FASTA string with 3–6 sequences."""
    n = draw(st.integers(min_value=3, max_value=6))
    blocks = [draw(fasta_sequence_strategy()) for _ in range(n)]
    # Ensure unique IDs by prefixing with index
    lines = []
    for i, block in enumerate(blocks):
        header, seq = block.split('\n', 1)
        lines.append(f'{header}_{i}\n{seq}')
    return ''.join(lines)


# ---------------------------------------------------------------------------
# Property 1 — FASTA round-trip preservation
# Feature: django-deployment, Property 1: FASTA round-trip preservation
# ---------------------------------------------------------------------------

@given(content=valid_fasta_content_strategy())
@h_settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
def test_fasta_roundtrip_preservation(content):
    """Parse → align → export_aligned_fasta → re-parse; IDs must be set-equal."""
    parser = SequenceParser()
    sequences = parser.parse_file(content)
    if len(sequences) < 3:
        return  # skip degenerate cases

    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)

    exporter = OutputExporter()
    fasta_out = exporter.export_aligned_fasta(alignment)

    re_parsed = parser.parse_file(fasta_out)
    original_ids = {s.id for s in sequences}
    roundtrip_ids = {s.id for s in re_parsed}
    assert original_ids == roundtrip_ids, (
        f"ID mismatch after round-trip: original={original_ids}, got={roundtrip_ids}"
    )


# ---------------------------------------------------------------------------
# Property 2 — Newick export is re-parseable
# Feature: django-deployment, Property 2: Newick export is re-parseable
# ---------------------------------------------------------------------------

@given(content=valid_fasta_content_strategy())
@h_settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
def test_newick_export_reparseable(content):
    """Build tree → export_newick → re-parse; terminal count must match."""
    parser = SequenceParser()
    sequences = parser.parse_file(content)
    if len(sequences) < 3:
        return

    aligner = SequenceAligner()
    alignment = aligner.align_sequences(sequences)

    builder = TreeBuilder()
    tree = builder.build_tree(alignment, method='nj')

    exporter = OutputExporter()
    newick = exporter.export_newick(tree)

    reparsed = Phylo.read(StringIO(newick), 'newick')
    original_count = len(tree.get_terminals())
    reparsed_count = len(reparsed.get_terminals())
    assert original_count == reparsed_count, (
        f"Terminal count mismatch: original={original_count}, reparsed={reparsed_count}"
    )
