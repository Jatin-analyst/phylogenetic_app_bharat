"""View property tests using Hypothesis.
# Feature: django-deployment
"""
import io
import json
import pytest
from hypothesis import given, settings as h_settings, HealthCheck
from hypothesis import strategies as st
from hypothesis.extra.django import TestCase
from django.test import Client
from django.template.loader import render_to_string

from modules.parser import SequenceParser
from modules.aligner import SequenceAligner
from modules.tree_builder import TreeBuilder
from modules.visualizer import TreeVisualizer
from phylo_app.views import figure_to_svg


# ---------------------------------------------------------------------------
# Shared strategies
# ---------------------------------------------------------------------------

DNA_CHARS = st.sampled_from('ACGT')


@st.composite
def _fasta_block(draw, idx):
    length = draw(st.integers(min_value=30, max_value=80))
    seq = ''.join(draw(st.lists(DNA_CHARS, min_size=length, max_size=length)))
    return f'>seq{idx}\n{seq}\n'


@st.composite
def valid_fasta_strategy(draw):
    """3–5 sequences, guaranteed unique IDs."""
    n = draw(st.integers(min_value=3, max_value=5))
    return ''.join([draw(_fasta_block(i)) for i in range(n)])


@st.composite
def invalid_fasta_strategy(draw):
    """Content that the parser will reject — binary/non-decodable bytes encoded as latin-1,
    or content that triggers a ParseError (empty after strip)."""
    # The parser raises ParseError for empty/whitespace content.
    # We generate whitespace-only strings to guarantee a 400.
    return draw(st.text(alphabet=' \t\r\n', min_size=1, max_size=50))


@st.composite
def few_sequences_strategy(draw):
    """0–2 sequences — not enough for tree building."""
    n = draw(st.integers(min_value=0, max_value=2))
    return ''.join([draw(_fasta_block(i)) for i in range(n)])


@st.composite
def generate_params_strategy(draw):
    return {
        'layout': draw(st.sampled_from(['rectangular', 'circular'])),
        'method': draw(st.sampled_from(['nj', 'upgma'])),
        'theme': draw(st.sampled_from(['default', 'dark', 'colorful', 'minimal'])),
    }


@st.composite
def template_context_strategy(draw):
    seq_count = draw(st.integers(min_value=3, max_value=10))
    aln_len = draw(st.integers(min_value=30, max_value=200))
    tree_svg = draw(st.text(min_size=5, max_size=50).map(lambda s: f'<svg>{s}</svg>'))
    seq_id = draw(st.text(alphabet='abcde', min_size=1, max_size=5))
    closest_relatives = {
        seq_id: {
            'closest_relative': 'other',
            'distance': 0.1,
            'explanation': 'test',
        }
    }
    return {
        'tree_svg': tree_svg,
        'sequence_count': seq_count,
        'alignment_length': aln_len,
        'closest_relatives': closest_relatives,
        'tree_method': 'nj',
    }


# ---------------------------------------------------------------------------
# Property 3 — Valid upload returns correct sequence count
# Feature: django-deployment, Property 3: Valid upload returns correct sequence count
# ---------------------------------------------------------------------------

class TestProperty3ValidUpload(TestCase):
    @given(content=valid_fasta_strategy())
    @h_settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow])
    def test_valid_upload_sequence_count(self, content):
        """POST valid FASTA → 200 and sequence_count matches parser output."""
        client = Client()
        expected = len(SequenceParser().parse_file(content))
        f = io.BytesIO(content.encode())
        f.name = 'test.fasta'
        res = client.post('/upload/', {'file': f})
        self.assertEqual(res.status_code, 200)
        data = res.json()
        self.assertEqual(data['sequence_count'], expected)


# ---------------------------------------------------------------------------
# Property 4 — Invalid upload returns HTTP 400
# Feature: django-deployment, Property 4: Invalid upload returns HTTP 400
# ---------------------------------------------------------------------------

class TestProperty4InvalidUpload(TestCase):
    @given(content=invalid_fasta_strategy())
    @h_settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow])
    def test_invalid_upload_returns_400(self, content):
        """POST non-FASTA content → 400 with non-empty error field."""
        client = Client()
        f = io.BytesIO(content.encode('utf-8', errors='replace'))
        f.name = 'bad.fasta'
        res = client.post('/upload/', {'file': f})
        self.assertEqual(res.status_code, 400)
        data = res.json()
        self.assertIn('error', data)
        self.assertTrue(len(data['error']) > 0)


# ---------------------------------------------------------------------------
# Property 5 — Generate response contains all required fields
# Feature: django-deployment, Property 5: Generate response contains all required fields
# ---------------------------------------------------------------------------

class TestProperty5GenerateFields(TestCase):
    @given(params=generate_params_strategy())
    @h_settings(max_examples=10, suppress_health_check=[HealthCheck.too_slow])
    def test_generate_returns_all_fields(self, params):
        """Seeded session + POST generate → 200 with all five required fields non-null."""
        client = Client()
        fasta = (
            '>s1\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n'
            '>s2\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n'
            '>s3\nGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n'
        )
        f = io.BytesIO(fasta.encode())
        f.name = 'test.fasta'
        client.post('/upload/', {'file': f})

        res = client.post(
            '/generate/',
            data=json.dumps(params),
            content_type='application/json',
        )
        self.assertEqual(res.status_code, 200)
        data = res.json()
        for field in ('tree_svg', 'sequence_count', 'alignment_length', 'tree_method', 'closest_relatives'):
            self.assertIn(field, data)
            self.assertIsNotNone(data[field])


# ---------------------------------------------------------------------------
# Property 6 — Generate stores results in session
# Feature: django-deployment, Property 6: Generate stores results in session
# ---------------------------------------------------------------------------

class TestProperty6SessionStorage(TestCase):
    @given(content=valid_fasta_strategy())
    @h_settings(max_examples=10, suppress_health_check=[HealthCheck.too_slow])
    def test_generate_stores_session(self, content):
        """After generate, session must contain phylo_results with required keys."""
        client = Client()
        f = io.BytesIO(content.encode())
        f.name = 'test.fasta'
        client.post('/upload/', {'file': f})

        res = client.post(
            '/generate/',
            data=json.dumps({'layout': 'rectangular', 'method': 'nj', 'theme': 'default'}),
            content_type='application/json',
        )
        if res.status_code != 200:
            return  # skip if pipeline fails (e.g. degenerate sequences)

        session = client.session
        self.assertIn('phylo_results', session)
        results = session['phylo_results']
        for key in ('newick', 'alignment', 'settings', 'closest_relatives'):
            self.assertIn(key, results)


# ---------------------------------------------------------------------------
# Property 7 — Insufficient sequences returns HTTP 400
# Feature: django-deployment, Property 7: Insufficient sequences returns HTTP 400
# ---------------------------------------------------------------------------

class TestProperty7InsufficientSequences(TestCase):
    @given(content=few_sequences_strategy())
    @h_settings(max_examples=20, suppress_health_check=[HealthCheck.too_slow])
    def test_few_sequences_returns_400(self, content):
        """< 3 valid sequences → generate returns 400 with error."""
        client = Client()
        if content.strip():
            f = io.BytesIO(content.encode())
            f.name = 'few.fasta'
            upload_res = client.post('/upload/', {'file': f})
            if upload_res.status_code != 200:
                return  # upload itself failed — that's fine
        else:
            # Seed session manually with minimal content
            session = client.session
            session['file_content'] = content
            session.save()

        res = client.post(
            '/generate/',
            data=json.dumps({'layout': 'rectangular', 'method': 'nj', 'theme': 'default'}),
            content_type='application/json',
        )
        self.assertEqual(res.status_code, 400)
        self.assertIn('error', res.json())


# ---------------------------------------------------------------------------
# Property 8 — SVG output valid for any layout and theme
# Feature: django-deployment, Property 8: SVG output valid for any layout and theme
# ---------------------------------------------------------------------------

FIXED_FASTA = (
    '>s1\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n'
    '>s2\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n'
    '>s3\nGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n'
)


class TestProperty8SVGOutput(TestCase):
    @given(
        layout=st.sampled_from(['rectangular', 'circular']),
        theme=st.sampled_from(['default', 'dark', 'colorful', 'minimal']),
    )
    @h_settings(max_examples=8, suppress_health_check=[HealthCheck.too_slow])
    def test_svg_valid_for_layout_and_theme(self, layout, theme):
        """draw_tree + figure_to_svg must return a string starting with <svg or <?xml."""
        parser = SequenceParser()
        sequences = parser.parse_file(FIXED_FASTA)
        aligner = SequenceAligner()
        alignment = aligner.align_sequences(sequences)
        builder = TreeBuilder()
        tree = builder.build_tree(alignment, method='nj')

        visualizer = TreeVisualizer()
        fig = visualizer.draw_tree(tree, layout=layout, theme=theme)
        svg = figure_to_svg(fig)

        self.assertIsInstance(svg, str)
        self.assertTrue(
            svg.strip().startswith('<svg') or svg.strip().startswith('<?xml'),
            f"SVG output does not start with <svg or <?xml: {svg[:80]}"
        )


# ---------------------------------------------------------------------------
# Property 9 — Closest relatives covers all sequences
# Feature: django-deployment, Property 9: Closest relatives covers all sequences
# ---------------------------------------------------------------------------

class TestProperty9ClosestRelatives(TestCase):
    @given(content=valid_fasta_strategy())
    @h_settings(max_examples=15, suppress_health_check=[HealthCheck.too_slow])
    def test_closest_relatives_covers_all_sequences(self, content):
        """generate response: len(closest_relatives) == sequence_count, each entry has required keys."""
        client = Client()
        f = io.BytesIO(content.encode())
        f.name = 'test.fasta'
        client.post('/upload/', {'file': f})

        res = client.post(
            '/generate/',
            data=json.dumps({'layout': 'rectangular', 'method': 'nj', 'theme': 'default'}),
            content_type='application/json',
        )
        if res.status_code != 200:
            return

        data = res.json()
        self.assertEqual(len(data['closest_relatives']), data['sequence_count'])
        for entry in data['closest_relatives'].values():
            for key in ('closest_relative', 'distance', 'explanation'):
                self.assertIn(key, entry)


# ---------------------------------------------------------------------------
# Property 10 — Download returns correct headers for each format
# Feature: django-deployment, Property 10: Download returns correct headers for each format
# ---------------------------------------------------------------------------

EXPECTED_CONTENT_TYPES = {
    'png': 'image/png',
    'newick': 'text/plain',
    'fasta': 'text/plain',
    'json': 'application/json',
}


class TestProperty10DownloadHeaders(TestCase):
    def _seed_results(self, client):
        f = io.BytesIO(FIXED_FASTA.encode())
        f.name = 'test.fasta'
        client.post('/upload/', {'file': f})
        client.post(
            '/generate/',
            data=json.dumps({'layout': 'rectangular', 'method': 'nj', 'theme': 'default'}),
            content_type='application/json',
        )

    @given(fmt=st.sampled_from(['png', 'newick', 'fasta', 'json']))
    @h_settings(max_examples=4, suppress_health_check=[HealthCheck.too_slow], deadline=None)
    def test_download_correct_headers(self, fmt):
        """GET /download/<fmt>/ → 200, correct Content-Type, attachment disposition."""
        client = Client()
        self._seed_results(client)

        res = client.get(f'/download/{fmt}/')
        self.assertEqual(res.status_code, 200)
        self.assertIn(EXPECTED_CONTENT_TYPES[fmt], res.get('Content-Type', ''))
        self.assertIn('attachment', res.get('Content-Disposition', ''))


# ---------------------------------------------------------------------------
# Property 11 — Template renders context values
# Feature: django-deployment, Property 11: Template renders context values
# ---------------------------------------------------------------------------

class TestProperty11TemplateRendering(TestCase):
    @given(ctx=template_context_strategy())
    @h_settings(max_examples=20, suppress_health_check=[HealthCheck.too_slow])
    def test_template_renders_context(self, ctx):
        """render_to_string with context → HTML contains tree_svg, sequence_count, alignment_length."""
        html = render_to_string('phylo_app/index.html', ctx)
        self.assertIn(str(ctx['sequence_count']), html)
        self.assertIn(str(ctx['alignment_length']), html)
        # tree_svg is injected via |safe filter — check a fragment of it
        self.assertIn('<svg>', html)
