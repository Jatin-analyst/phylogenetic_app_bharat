import io
import json
import os
from io import StringIO

import matplotlib.pyplot as plt
from django.shortcuts import render
from django.http import HttpResponse, JsonResponse, FileResponse, HttpResponseNotAllowed
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from modules.parser import SequenceParser, ParseError
from modules.validator import SequenceValidator
from modules.aligner import SequenceAligner, AlignmentError
from modules.tree_builder import TreeBuilder, TreeBuildError
from modules.visualizer import TreeVisualizer
from modules.exporter import OutputExporter
from modules.data_models import Settings


# ---------------------------------------------------------------------------
# Task 6 — Session serialization helpers
# ---------------------------------------------------------------------------

def serialize_results(tree, alignment, settings, closest_relatives) -> dict:
    """Serialize pipeline results to JSON-safe dict for session storage."""
    exporter = OutputExporter()
    newick = exporter.export_newick(tree)
    alignment_list = [
        {"id": rec.id, "sequence": str(rec.seq), "description": rec.description}
        for rec in alignment
    ]
    settings_dict = {
        "layout": settings.layout,
        "theme": settings.theme,
        "tree_method": settings.tree_method,
        "dpi": settings.dpi,
    }
    relatives_dict = {
        seq_id: {
            "closest_relative": rel.closest_relative,
            "distance": rel.distance,
            "explanation": rel.explanation,
        }
        for seq_id, rel in closest_relatives.items()
    }
    return {
        "newick": newick,
        "alignment": alignment_list,
        "settings": settings_dict,
        "closest_relatives": relatives_dict,
    }


def deserialize_results(data: dict):
    """Reconstruct Bio.Phylo tree and MultipleSeqAlignment from session data."""
    tree = Phylo.read(StringIO(data["newick"]), "newick")
    records = [
        SeqRecord(Seq(item["sequence"]), id=item["id"], description=item.get("description", ""))
        for item in data["alignment"]
    ]
    alignment = MultipleSeqAlignment(records)
    settings = Settings(**data["settings"])
    return tree, alignment, settings


# ---------------------------------------------------------------------------
# Task 7 — figure_to_svg adapter
# ---------------------------------------------------------------------------

def figure_to_svg(fig) -> str:
    """Convert a matplotlib Figure to an SVG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format='svg', bbox_inches='tight')
    buf.seek(0)
    svg = buf.getvalue().decode('utf-8')
    plt.close(fig)
    return svg


# ---------------------------------------------------------------------------
# Task 8 — index_view
# ---------------------------------------------------------------------------

def index_view(request):
    """Render the main page. Populate context from session if results exist."""
    context = {
        "tree_svg": None,
        "sequence_count": None,
        "alignment_length": None,
        "closest_relatives": None,
        "tree_method": None,
    }
    results = request.session.get("phylo_results")
    if results:
        context["tree_svg"] = results.get("tree_svg")
        context["sequence_count"] = results.get("sequence_count")
        context["alignment_length"] = results.get("alignment_length")
        context["closest_relatives"] = results.get("closest_relatives")
        context["tree_method"] = results.get("tree_method")
    return render(request, "phylo_app/index.html", context)


# ---------------------------------------------------------------------------
# Task 9 — upload_view
# ---------------------------------------------------------------------------

def upload_view(request):
    """Handle file upload. Validates and stores file content in session."""
    if request.method != "POST":
        return HttpResponseNotAllowed(["POST"])

    uploaded_file = request.FILES.get("file")
    if not uploaded_file:
        return JsonResponse({"error": "No file provided."}, status=400)

    name = uploaded_file.name.lower()
    if not (name.endswith(".fasta") or name.endswith(".fa") or name.endswith(".txt")):
        return JsonResponse(
            {"error": "Unsupported file type. Please upload .fasta, .fa, or .txt files."},
            status=400,
        )

    if uploaded_file.size == 0:
        return JsonResponse({"error": "Uploaded file is empty."}, status=400)
    if uploaded_file.size > 10 * 1024 * 1024:
        return JsonResponse({"error": "File too large. Maximum size is 10 MB."}, status=400)

    try:
        content = uploaded_file.read().decode("utf-8")
    except UnicodeDecodeError:
        return JsonResponse({"error": "File could not be decoded as UTF-8."}, status=400)

    try:
        parser = SequenceParser()
        sequences = parser.parse_file(content)
    except ParseError as e:
        return JsonResponse({"error": str(e)}, status=400)

    request.session["file_content"] = content
    request.session.pop("phylo_results", None)

    return JsonResponse({"filename": uploaded_file.name, "sequence_count": len(sequences)})


# ---------------------------------------------------------------------------
# Task 10 — generate_view
# ---------------------------------------------------------------------------

def generate_view(request):
    """Run the full pipeline and return results as JSON."""
    if request.method != "POST":
        return HttpResponseNotAllowed(["POST"])

    file_content = request.session.get("file_content")
    if not file_content:
        return JsonResponse(
            {"error": "No file uploaded. Please upload a FASTA or TXT file first."},
            status=400,
        )

    try:
        body = json.loads(request.body)
    except Exception:
        body = request.POST

    layout = body.get("layout", "rectangular")
    method = body.get("method", "nj")
    theme = body.get("theme", "default")

    if layout not in ("rectangular", "circular"):
        layout = "rectangular"
    if method not in ("nj", "upgma"):
        method = "nj"
    if theme not in ("default", "dark", "colorful", "minimal"):
        theme = "default"

    settings = Settings(layout=layout, theme=theme, tree_method=method, dpi=300)

    try:
        parser = SequenceParser()
        sequences = parser.parse_file(file_content)

        validator = SequenceValidator(min_length=10)
        validation_result = validator.validate_sequences(sequences)
        valid_sequences = validation_result.valid_sequences

        if len(valid_sequences) < 3:
            return JsonResponse(
                {
                    "error": (
                        f"At least 3 valid sequences are required. "
                        f"Found {len(valid_sequences)} valid sequence(s)."
                    )
                },
                status=400,
            )

        aligner = SequenceAligner()
        alignment = aligner.align_sequences(valid_sequences)

        builder = TreeBuilder()
        tree = builder.build_tree(alignment, method=method)

        visualizer = TreeVisualizer()
        fig = visualizer.draw_tree(tree, layout=layout, theme=theme)
        tree_svg = figure_to_svg(fig)

        exporter = OutputExporter()
        relatives = exporter.calculate_closest_relatives(tree)

        serialized = serialize_results(tree, alignment, settings, relatives)
        serialized["tree_svg"] = tree_svg
        serialized["sequence_count"] = len(valid_sequences)
        serialized["alignment_length"] = alignment.get_alignment_length()
        serialized["tree_method"] = method
        request.session["phylo_results"] = serialized
        request.session.modified = True

        return JsonResponse({
            "tree_svg": tree_svg,
            "sequence_count": len(valid_sequences),
            "alignment_length": alignment.get_alignment_length(),
            "tree_method": method,
            "closest_relatives": serialized["closest_relatives"],
        })

    except AlignmentError as e:
        return JsonResponse({"error": str(e)}, status=400)
    except TreeBuildError as e:
        return JsonResponse({"error": str(e)}, status=400)
    except Exception:
        return JsonResponse({"error": "An unexpected error occurred. Please try again."}, status=500)


# ---------------------------------------------------------------------------
# Task 11 — download_view
# ---------------------------------------------------------------------------

def download_view(request, format):
    """Return a downloadable file for the given format."""
    results = request.session.get("phylo_results")
    if not results:
        return JsonResponse(
            {"error": "No results available. Please generate a tree first."},
            status=400,
        )

    try:
        tree, alignment, settings = deserialize_results(results)
        exporter = OutputExporter()

        if format == "png":
            visualizer = TreeVisualizer()
            fig = visualizer.draw_tree(tree, layout=settings.layout, theme=settings.theme)
            png_bytes = visualizer.get_png_bytes(fig, dpi=settings.dpi)
            response = HttpResponse(png_bytes, content_type="image/png")
            response["Content-Disposition"] = 'attachment; filename="phylogenetic_tree.png"'
            return response

        elif format == "newick":
            newick_str = exporter.export_newick(tree)
            response = HttpResponse(newick_str, content_type="text/plain")
            response["Content-Disposition"] = 'attachment; filename="tree.nwk"'
            return response

        elif format == "fasta":
            fasta_str = exporter.export_aligned_fasta(alignment)
            response = HttpResponse(fasta_str, content_type="text/plain")
            response["Content-Disposition"] = 'attachment; filename="aligned.fasta"'
            return response

        elif format == "json":
            json_str = json.dumps(results["closest_relatives"], indent=2)
            response = HttpResponse(json_str, content_type="application/json")
            response["Content-Disposition"] = 'attachment; filename="closest_relatives.json"'
            return response

        else:
            return JsonResponse({"error": "Unknown download format."}, status=404)

    except Exception as e:
        return JsonResponse({"error": f"Download failed: {str(e)}"}, status=500)
