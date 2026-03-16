"""Unit tests for phylo_app views, URL routing, and configuration."""
import io
import json
import yaml
import pytest
from django.test import TestCase, Client
from django.urls import resolve, reverse

from phylo_app import views


# ---------------------------------------------------------------------------
# 15.1 URL routing
# ---------------------------------------------------------------------------

class TestURLRouting(TestCase):
    def test_index_resolves(self):
        match = resolve('/')
        self.assertEqual(match.func, views.index_view)

    def test_upload_resolves(self):
        match = resolve('/upload/')
        self.assertEqual(match.func, views.upload_view)

    def test_generate_resolves(self):
        match = resolve('/generate/')
        self.assertEqual(match.func, views.generate_view)

    def test_download_resolves(self):
        match = resolve('/download/png/')
        self.assertEqual(match.func, views.download_view)


# ---------------------------------------------------------------------------
# 15.2 Upload edge cases
# ---------------------------------------------------------------------------

class TestUploadEdgeCases(TestCase):
    def setUp(self):
        self.client = Client()

    def test_empty_file_returns_400(self):
        f = io.BytesIO(b'')
        f.name = 'empty.fasta'
        res = self.client.post('/upload/', {'file': f})
        self.assertEqual(res.status_code, 400)
        self.assertIn('error', res.json())

    def test_non_utf8_returns_400(self):
        f = io.BytesIO(b'\xff\xfe invalid bytes')
        f.name = 'bad.fasta'
        res = self.client.post('/upload/', {'file': f})
        self.assertEqual(res.status_code, 400)
        self.assertIn('error', res.json())

    def test_unsupported_extension_returns_400(self):
        f = io.BytesIO(b'>seq1\nACGT\n')
        f.name = 'sequences.csv'
        res = self.client.post('/upload/', {'file': f})
        self.assertEqual(res.status_code, 400)
        self.assertIn('error', res.json())

    def test_large_file_returns_400(self):
        big = io.BytesIO(b'A' * (11 * 1024 * 1024))
        big.name = 'big.fasta'
        res = self.client.post('/upload/', {'file': big})
        self.assertEqual(res.status_code, 400)
        self.assertIn('error', res.json())


# ---------------------------------------------------------------------------
# 15.3 Download before generate
# ---------------------------------------------------------------------------

class TestDownloadBeforeGenerate(TestCase):
    def test_download_png_without_results_returns_400(self):
        client = Client()
        res = client.get('/download/png/')
        self.assertEqual(res.status_code, 400)
        self.assertIn('error', res.json())


# ---------------------------------------------------------------------------
# 15.4 Unknown download format
# ---------------------------------------------------------------------------

class TestUnknownDownloadFormat(TestCase):
    def setUp(self):
        self.client = Client()
        # Seed session with fake results so we reach the format check
        session = self.client.session
        session['phylo_results'] = {
            'newick': '(A:0.1,B:0.2);',
            'alignment': [
                {'id': 'A', 'sequence': 'ACGT', 'description': ''},
                {'id': 'B', 'sequence': 'ACGT', 'description': ''},
            ],
            'settings': {'layout': 'rectangular', 'theme': 'default', 'tree_method': 'nj', 'dpi': 300},
            'closest_relatives': {},
            'tree_svg': '<svg></svg>',
            'sequence_count': 2,
            'alignment_length': 4,
            'tree_method': 'nj',
        }
        session.save()

    def test_unknown_format_returns_404(self):
        res = self.client.get('/download/xyz/')
        self.assertEqual(res.status_code, 404)
        self.assertIn('error', res.json())


# ---------------------------------------------------------------------------
# 15.5 Settings configuration
# ---------------------------------------------------------------------------

class TestSettingsConfig(TestCase):
    def test_secret_key_present(self):
        from django.conf import settings
        self.assertTrue(hasattr(settings, 'SECRET_KEY'))
        self.assertIsInstance(settings.SECRET_KEY, str)
        self.assertTrue(len(settings.SECRET_KEY) > 0)

    def test_debug_is_bool(self):
        from django.conf import settings
        self.assertIsInstance(settings.DEBUG, bool)

    def test_allowed_hosts_is_list(self):
        from django.conf import settings
        self.assertIsInstance(settings.ALLOWED_HOSTS, list)
        self.assertIn('localhost', settings.ALLOWED_HOSTS)

    def test_session_engine_is_file(self):
        from django.conf import settings
        self.assertEqual(settings.SESSION_ENGINE, 'django.contrib.sessions.backends.file')

    def test_static_root_configured(self):
        from django.conf import settings
        self.assertTrue(hasattr(settings, 'STATIC_ROOT'))
        self.assertIsNotNone(settings.STATIC_ROOT)

    def test_whitenoise_in_middleware(self):
        from django.conf import settings
        middleware_str = ' '.join(settings.MIDDLEWARE)
        self.assertIn('whitenoise', middleware_str.lower())


# ---------------------------------------------------------------------------
# 15.6 requirements.txt contents
# ---------------------------------------------------------------------------

class TestRequirementsTxt(TestCase):
    def _read_requirements(self):
        import pathlib
        req_path = pathlib.Path(__file__).resolve().parents[2] / 'requirements.txt'
        return req_path.read_text().lower()

    def test_django_present(self):
        self.assertIn('django', self._read_requirements())

    def test_gunicorn_present(self):
        self.assertIn('gunicorn', self._read_requirements())

    def test_whitenoise_present(self):
        self.assertIn('whitenoise', self._read_requirements())

    def test_biopython_present(self):
        self.assertIn('biopython', self._read_requirements())

    def test_hypothesis_present(self):
        self.assertIn('hypothesis', self._read_requirements())

    def test_pytest_django_present(self):
        self.assertIn('pytest-django', self._read_requirements())

    def test_streamlit_absent(self):
        self.assertNotIn('streamlit', self._read_requirements())


# ---------------------------------------------------------------------------
# 15.7 render.yaml structure
# ---------------------------------------------------------------------------

class TestRenderYaml(TestCase):
    def _load_render_yaml(self):
        import pathlib
        render_path = pathlib.Path(__file__).resolve().parents[2] / 'render.yaml'
        with open(render_path) as f:
            return yaml.safe_load(f)

    def test_service_present(self):
        data = self._load_render_yaml()
        self.assertIn('services', data)
        self.assertGreater(len(data['services']), 0)

    def test_service_type_web(self):
        data = self._load_render_yaml()
        svc = data['services'][0]
        self.assertEqual(svc.get('type'), 'web')

    def test_build_command_present(self):
        data = self._load_render_yaml()
        svc = data['services'][0]
        self.assertIn('buildCommand', svc)
        self.assertIn('pip install', svc['buildCommand'])

    def test_start_command_present(self):
        data = self._load_render_yaml()
        svc = data['services'][0]
        self.assertIn('startCommand', svc)
        self.assertIn('gunicorn', svc['startCommand'])

    def test_required_env_vars_present(self):
        data = self._load_render_yaml()
        svc = data['services'][0]
        env_keys = {e['key'] for e in svc.get('envVars', [])}
        for required in ('SECRET_KEY', 'DEBUG', 'ALLOWED_HOSTS'):
            self.assertIn(required, env_keys)


# ---------------------------------------------------------------------------
# 15.8 End-to-end smoke test
# ---------------------------------------------------------------------------

SMOKE_FASTA = b""">seq1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>seq2
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
>seq3
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG
"""


class TestEndToEndSmoke(TestCase):
    def setUp(self):
        self.client = Client()

    def test_upload_returns_200(self):
        f = io.BytesIO(SMOKE_FASTA)
        f.name = 'smoke.fasta'
        res = self.client.post('/upload/', {'file': f})
        self.assertEqual(res.status_code, 200)
        data = res.json()
        self.assertIn('sequence_count', data)
        self.assertEqual(data['sequence_count'], 3)

    def test_generate_returns_200_with_required_fields(self):
        # Upload first
        f = io.BytesIO(SMOKE_FASTA)
        f.name = 'smoke.fasta'
        self.client.post('/upload/', {'file': f})

        res = self.client.post(
            '/generate/',
            data=json.dumps({'layout': 'rectangular', 'method': 'nj', 'theme': 'default'}),
            content_type='application/json',
        )
        self.assertEqual(res.status_code, 200)
        data = res.json()
        for field in ('tree_svg', 'sequence_count', 'alignment_length', 'tree_method', 'closest_relatives'):
            self.assertIn(field, data)
            self.assertIsNotNone(data[field])

    def test_download_newick_returns_200_with_correct_content_type(self):
        # Upload + generate first
        f = io.BytesIO(SMOKE_FASTA)
        f.name = 'smoke.fasta'
        self.client.post('/upload/', {'file': f})
        self.client.post(
            '/generate/',
            data=json.dumps({'layout': 'rectangular', 'method': 'nj', 'theme': 'default'}),
            content_type='application/json',
        )

        res = self.client.get('/download/newick/')
        self.assertEqual(res.status_code, 200)
        self.assertIn('text/plain', res.get('Content-Type', ''))
        self.assertIn('attachment', res.get('Content-Disposition', ''))
