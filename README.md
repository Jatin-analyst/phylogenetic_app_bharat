# 🧬 Phylogenetics App

A beginner-friendly web application for phylogenetic analysis built with Streamlit. Build phylogenetic trees from your sequence data with ease!

## ✨ Features

- 📁 **Multiple Input Formats**: FASTA, multi-FASTA, TXT (raw or named sequences)
- 🧹 **Automated Sequence Cleaning**: Removes invalid characters and validates sequences
- 🔄 **Multi-Sequence Alignment**: Progressive alignment using Biopython
- 🌳 **Phylogenetic Tree Construction**: Neighbor-Joining (NJ) and UPGMA methods
- 🎨 **Multiple Visualization Styles**: Rectangular and circular tree layouts
- 🎨 **Color Themes**: Default, dark, colorful, and minimal themes
- 📊 **High-Resolution Export**: PNG images at 300+ DPI
- 📄 **Newick Format Export**: Standard phylogenetic tree format
- 🤝 **Closest Relatives Analysis**: JSON export with distance metrics and explanations
- ⚡ **Performance Optimized**: Caching for faster repeated operations

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone <repository-url>
cd phylogenetics-app

# Install dependencies
pip install -r requirements.txt
```

### Running Locally

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`

### Using the App

1. **Upload Your Sequences**: Choose a FASTA or TXT file
2. **Configure Settings**: Select tree layout, color theme, and building method in the sidebar
3. **Build Tree**: Click the "Build Tree" button
4. **View Results**: Explore the phylogenetic tree and closest relatives
5. **Download**: Get PNG, Newick, aligned FASTA, and JSON outputs

## 📁 Sample Data

Sample input files are provided in the `sample_data/` directory:
- `sample.fasta`: Protein sequences (hemoglobin from different species)
- `sample_dna.txt`: DNA sequences in named TXT format

## 🧪 Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run specific test modules
pytest tests/test_parser.py -v
pytest tests/test_integration.py -v

# Run with coverage
pytest tests/ --cov=modules --cov-report=html
```

**Test Coverage**: 50 tests covering all modules with property-based testing using Hypothesis.

## 📂 Project Structure

```
phylogenetics-app/
├── app.py                      # Main Streamlit application
├── modules/                    # Core modules
│   ├── __init__.py
│   ├── data_models.py         # Data classes (Sequence, Settings, etc.)
│   ├── parser.py              # Sequence parser (FASTA, TXT)
│   ├── validator.py           # Sequence validator and cleaner
│   ├── aligner.py             # Multi-sequence alignment
│   ├── tree_builder.py        # Phylogenetic tree construction
│   ├── visualizer.py          # Tree visualization with matplotlib
│   └── exporter.py            # Output exporter (Newick, FASTA, JSON)
├── tests/                      # Comprehensive test suite
│   ├── test_parser.py         # Parser tests (9 tests)
│   ├── test_validator.py      # Validator tests (9 tests)
│   ├── test_aligner.py        # Aligner tests (6 tests)
│   ├── test_tree_builder.py   # Tree builder tests (7 tests)
│   ├── test_exporter.py       # Exporter tests (8 tests)
│   ├── test_visualizer.py     # Visualizer tests (6 tests)
│   └── test_integration.py    # Integration tests (5 tests)
├── sample_data/                # Sample input files
│   ├── sample.fasta
│   └── sample_dna.txt
├── requirements.txt            # Python dependencies
└── README.md                   # This file
```

## 🔧 Configuration

### Supported File Formats

**FASTA Format:**
```
>Sequence_1
ACGTACGTACGT
>Sequence_2
ACGTACGGACGT
```

**TXT Format (Named):**
```
Seq_1:ACGTACGTACGT
Seq_2:ACGTACGGACGT
```

**TXT Format (Raw):**
```
ACGTACGTACGT
ACGTACGGACGT
```

### Tree Building Methods

- **Neighbor-Joining (NJ)**: Generally more accurate, doesn't assume a molecular clock
- **UPGMA**: Simpler method, assumes a constant rate of evolution (molecular clock)

### Visualization Options

- **Layouts**: Rectangular (traditional), Circular (radial)
- **Themes**: Default (black on white), Dark (white on dark), Colorful (gradient colors), Minimal (gray tones)

## 🌐 Deployment

### Streamlit Cloud

1. Push your code to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your repository
4. Deploy!

The app is designed to run on Streamlit Cloud without any additional configuration.

### Requirements

- Python 3.8+
- All dependencies are pip-installable
- No system-level dependencies required
- Memory optimized for Streamlit Cloud (~1GB limit)

## 🧬 Technical Details

### Algorithms

- **Alignment**: Progressive pairwise alignment using Biopython
- **Distance Calculation**: Identity-based distance matrix
- **Tree Construction**: Biopython's Phylo module (NJ and UPGMA)
- **Visualization**: Matplotlib with custom styling

### Performance

- **Optimal**: 5-20 sequences, 500-5000 characters each
- **Maximum**: 100 sequences, 10,000 characters each
- **Caching**: Streamlit's @st.cache_data for expensive operations

### Testing

- **Property-Based Testing**: Using Hypothesis for comprehensive coverage
- **Integration Tests**: End-to-end workflow validation
- **45 Unit Tests**: Covering all modules
- **5 Integration Tests**: Complete workflow scenarios

## 📝 License

MIT License - feel free to use this project for any purpose.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📧 Support

For issues or questions, please open an issue on GitHub.

---

Built with ❤️ using Streamlit, Biopython, and Python
