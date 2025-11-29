# Design Document

## Overview

The Phylogenetics App is a Streamlit-based web application that provides an end-to-end pipeline for phylogenetic analysis. The system follows a linear workflow: file upload → sequence parsing → cleaning/validation → alignment → tree construction → visualization → output generation. The architecture emphasizes simplicity, using lightweight Python libraries compatible with Streamlit Cloud deployment constraints.

The application uses Biopython for sequence handling and alignment, and either Biopython's Phylo module or ETE3 for tree construction and visualization. All processing occurs server-side with results cached to improve performance for repeated operations.

## Architecture

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     Streamlit Frontend                       │
│  (File Upload, Controls, Visualization Display, Downloads)  │
└────────────────────┬────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│                   Application Layer                          │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐     │
│  │   Parser     │  │   Aligner    │  │ Tree Builder │     │
│  │   Module     │  │   Module     │  │    Module    │     │
│  └──────────────┘  └──────────────┘  └──────────────┘     │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐     │
│  │  Validator   │  │ Visualizer   │  │   Exporter   │     │
│  │   Module     │  │   Module     │  │    Module    │     │
│  └──────────────┘  └──────────────┘  └──────────────┘     │
└─────────────────────────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────┐
│                    Data Layer                                │
│         (Sequence Objects, Tree Objects, Cache)             │
└─────────────────────────────────────────────────────────────┘
```

### Technology Stack

- **Frontend Framework**: Streamlit (latest stable version)
- **Sequence Processing**: Biopython (Bio.SeqIO, Bio.Seq, Bio.AlignIO)
- **Alignment**: Biopython's pairwise2 or MUSCLE wrapper (if available)
- **Tree Construction**: Biopython Phylo with distance-based methods (UPGMA or Neighbor-Joining)
- **Visualization**: Matplotlib with Biopython Phylo.draw or Toytree (lightweight alternative)
- **File Handling**: Python standard library (io, tempfile)
- **Data Export**: JSON (standard library), Biopython Newick writer

### Deployment Considerations

- All dependencies must be pip-installable and lightweight
- No system-level dependencies (no MUSCLE binary, no FastTree binary)
- Use pure Python implementations where possible
- Streamlit Cloud has ~1GB memory limit - keep data structures efficient
- Session state management for caching intermediate results

## Components and Interfaces

### 1. Parser Module

**Responsibility**: Parse uploaded files and extract sequences with identifiers

**Interface**:
```python
class SequenceParser:
    def parse_file(file_content: str, file_type: str) -> List[Sequence]
    def detect_format(file_content: str) -> str
    def parse_fasta(content: str) -> List[Sequence]
    def parse_multi_fasta(content: str) -> List[Sequence]
    def parse_txt_raw(content: str) -> List[Sequence]
    def parse_txt_named(content: str) -> List[Sequence]
```

**Key Behaviors**:
- Auto-detect format when possible
- Generate default IDs for unnamed sequences (e.g., "Seq_1", "Seq_2")
- Return structured Sequence objects with id and sequence string
- Raise ParseError with descriptive messages for invalid formats

### 2. Validator Module

**Responsibility**: Clean and validate sequences for phylogenetic analysis

**Interface**:
```python
class SequenceValidator:
    def clean_sequence(seq: str) -> str
    def detect_sequence_type(sequences: List[str]) -> str  # "nucleotide" or "protein"
    def validate_sequences(sequences: List[Sequence]) -> ValidationResult
    def remove_invalid_characters(seq: str, seq_type: str) -> str
    def check_minimum_length(seq: str, min_length: int = 10) -> bool
```

**Key Behaviors**:
- Remove whitespace, numbers, and invalid characters
- Detect nucleotide (ACGT) vs protein sequences
- Flag sequences that are too short (< 10 characters)
- Return ValidationResult with cleaned sequences and summary of actions
- Preserve original IDs through cleaning process

### 3. Aligner Module

**Responsibility**: Perform multi-sequence alignment

**Interface**:
```python
class SequenceAligner:
    def align_sequences(sequences: List[Sequence]) -> Alignment
    def get_alignment_method() -> str  # Returns method name for display
```

**Key Behaviors**:
- Use Biopython's pairwise alignment for small datasets (< 10 sequences)
- For larger datasets, use progressive alignment approach
- Return Bio.Align.MultipleSeqAlignment object
- Handle alignment failures gracefully with informative errors
- Preserve sequence IDs in alignment

**Implementation Note**: Since MUSCLE/ClustalW binaries may not be available on Streamlit Cloud, implement a simple progressive alignment using pairwise alignments, or use Biopython's built-in methods.

### 4. Tree Builder Module

**Responsibility**: Construct phylogenetic tree from aligned sequences

**Interface**:
```python
class TreeBuilder:
    def build_tree(alignment: Alignment, method: str = "nj") -> Tree
    def calculate_distance_matrix(alignment: Alignment) -> DistanceMatrix
    def neighbor_joining(distance_matrix: DistanceMatrix) -> Tree
    def upgma(distance_matrix: DistanceMatrix) -> Tree
```

**Key Behaviors**:
- Calculate pairwise distances from alignment
- Use Neighbor-Joining (NJ) as default method (faster, more accurate for most cases)
- Provide UPGMA as alternative (simpler, assumes molecular clock)
- Return Biopython Tree object with branch lengths
- Handle edge cases (identical sequences, very divergent sequences)

### 5. Visualizer Module

**Responsibility**: Generate tree visualizations in multiple styles

**Interface**:
```python
class TreeVisualizer:
    def draw_tree(tree: Tree, layout: str, theme: str) -> Figure
    def draw_rectangular(tree: Tree, theme: str) -> Figure
    def draw_circular(tree: Tree, theme: str) -> Figure
    def apply_color_theme(fig: Figure, theme: str) -> Figure
    def save_high_res_png(fig: Figure, filename: str, dpi: int = 300)
```

**Key Behaviors**:
- Support "rectangular" and "circular" layouts
- Provide color themes: "default", "dark", "colorful", "minimal"
- Generate matplotlib Figure objects for display
- Ensure labels are readable and non-overlapping
- Export PNG at 300 DPI minimum
- Add hover information using Streamlit's native capabilities

**Color Themes**:
- Default: Black branches, black labels on white background
- Dark: White branches, white labels on dark background
- Colorful: Gradient colors based on branch depth
- Minimal: Gray branches, minimal styling

### 6. Exporter Module

**Responsibility**: Generate downloadable output files

**Interface**:
```python
class OutputExporter:
    def export_newick(tree: Tree) -> str
    def export_aligned_fasta(alignment: Alignment) -> str
    def export_closest_relatives(tree: Tree) -> str  # JSON format
    def calculate_closest_relatives(tree: Tree) -> Dict[str, RelativeInfo]
```

**Key Behaviors**:
- Write Newick format with branch lengths
- Write aligned sequences in FASTA format
- Generate JSON with closest relative information for each sequence
- Include distance metrics in closest relatives output
- Format JSON for human readability

**Closest Relatives JSON Structure**:
```json
{
  "Seq_1": {
    "closest_relative": "Seq_3",
    "distance": 0.042,
    "explanation": "Seq_1 is most closely related to Seq_3 with a genetic distance of 0.042"
  }
}
```

### 7. Streamlit UI Module

**Responsibility**: Provide user interface and orchestrate workflow

**Interface**:
```python
def main():
    # Streamlit app entry point
    
def render_upload_section() -> UploadedFile
def render_settings_section() -> Settings
def render_build_button() -> bool
def render_progress(step: str)
def render_results(tree: Tree, alignment: Alignment, outputs: Dict)
def render_downloads(outputs: Dict)
```

**Key Behaviors**:
- Single-page layout with clear sections
- File uploader at top
- Settings panel (layout, theme) in sidebar
- "Build Tree" button prominently displayed
- Progress indicators during processing
- Results section with interactive tree display
- Download buttons for all outputs
- Error messages displayed with st.error()
- Use st.cache_data for expensive operations

## Data Models

### Sequence
```python
@dataclass
class Sequence:
    id: str
    sequence: str
    description: str = ""
    
    def __len__(self) -> int:
        return len(self.sequence)
```

### ValidationResult
```python
@dataclass
class ValidationResult:
    valid_sequences: List[Sequence]
    invalid_sequences: List[Sequence]
    cleaning_summary: Dict[str, int]  # e.g., {"removed_spaces": 5, "removed_numbers": 2}
    sequence_type: str  # "nucleotide" or "protein"
```

### Settings
```python
@dataclass
class Settings:
    layout: str = "rectangular"  # "rectangular" or "circular"
    theme: str = "default"  # "default", "dark", "colorful", "minimal"
    tree_method: str = "nj"  # "nj" or "upgma"
    dpi: int = 300
```

### RelativeInfo
```python
@dataclass
class RelativeInfo:
    closest_relative: str
    distance: float
    explanation: str
```

### Tree
Using Biopython's Bio.Phylo.BaseTree.Tree object directly

### Alignment
Using Biopython's Bio.Align.MultipleSeqAlignment object directly


## Correctness Properties

*A property is a characteristic or behavior that should hold true across all valid executions of a system—essentially, a formal statement about what the system should do. Properties serve as the bridge between human-readable specifications and machine-verifiable correctness guarantees.*

### Property 1: Parsing preserves sequence count and identifiers

*For any* valid input file (FASTA, multi-FASTA, or TXT format), parsing should extract all sequences and assign each a unique identifier, with the number of output sequences matching the number of input sequences.

**Validates: Requirements 1.1, 1.2, 1.3, 1.4**

### Property 2: Invalid input rejection

*For any* invalid or malformed input file, the parser should raise an appropriate error and not return partial or corrupted sequence data.

**Validates: Requirements 1.5**

### Property 3: Sequence cleaning preserves valid characters

*For any* sequence containing a mix of valid and invalid characters, cleaning should remove only invalid characters while preserving all valid nucleotide or protein characters in their original order.

**Validates: Requirements 2.1**

### Property 4: Sequence type detection consistency

*For any* set of sequences of the same type (all nucleotide or all protein), the validator should correctly identify the sequence type consistently across all sequences.

**Validates: Requirements 2.2**

### Property 5: Short sequence filtering

*For any* set of sequences, all sequences shorter than the minimum length threshold should be excluded from the validation result, and all sequences meeting the threshold should be included.

**Validates: Requirements 2.3**

### Property 6: Alignment preserves sequences and identifiers

*For any* set of valid sequences, alignment should produce an output containing the same number of sequences with the same identifiers as the input, with each sequence having the same length (due to gap insertion).

**Validates: Requirements 3.1, 3.4**

### Property 7: Aligned output is valid FASTA

*For any* alignment result, exporting to FASTA format should produce a valid FASTA file that can be parsed back into the same number of sequences.

**Validates: Requirements 3.2**

### Property 8: Tree construction produces valid tree

*For any* valid alignment, tree construction should produce a tree object with the same number of leaf nodes as input sequences, where each leaf is labeled with a sequence identifier from the input.

**Validates: Requirements 4.1**

### Property 9: Newick round-trip consistency

*For any* constructed phylogenetic tree, exporting to Newick format and parsing back should preserve the tree topology and branch lengths within numerical precision.

**Validates: Requirements 4.2**

### Property 10: Branch lengths are non-negative

*For any* constructed phylogenetic tree, all branch lengths should be non-negative real numbers representing evolutionary distances.

**Validates: Requirements 4.4**

### Property 11: Layout changes preserve tree structure

*For any* phylogenetic tree, changing the visualization layout (rectangular to circular or vice versa) should not modify the underlying tree structure, topology, or branch lengths.

**Validates: Requirements 5.4**

### Property 12: PNG generation with minimum resolution

*For any* tree visualization, the generated PNG file should have a resolution of at least 300 DPI and be a valid PNG format that can be opened by standard image viewers.

**Validates: Requirements 6.1, 6.2**

### Property 13: Closest relatives completeness

*For any* phylogenetic tree with N sequences, the closest relatives JSON should contain exactly N entries, where each sequence has an identified closest relative and a non-negative distance value.

**Validates: Requirements 7.1, 7.2, 7.5**

## Error Handling

### Error Categories

1. **Input Errors**
   - Invalid file format
   - Empty files
   - Corrupted sequence data
   - Unsupported characters

2. **Validation Errors**
   - All sequences too short
   - Mixed sequence types (nucleotide and protein)
   - Insufficient sequences for tree building (< 3 sequences)

3. **Processing Errors**
   - Alignment failure (sequences too divergent)
   - Tree construction failure (numerical instability)
   - Memory limits exceeded

4. **Export Errors**
   - File write failures
   - Format conversion errors

### Error Handling Strategy

- All errors should be caught and displayed using Streamlit's `st.error()` with beginner-friendly messages
- Error messages should include:
  - What went wrong (in plain language)
  - Why it might have happened
  - Suggested next steps
- Processing should stop at the first error to avoid cascading failures
- Partial results should not be displayed if processing fails
- Use Python's exception hierarchy: custom exceptions inherit from appropriate base classes

### Example Error Messages

```python
# Input error
"Unable to parse file: The uploaded file doesn't appear to be in FASTA or TXT format. 
Please check that your file contains sequence data in one of the supported formats."

# Validation error
"Not enough sequences: Phylogenetic trees require at least 3 sequences. 
Your file contains only 2 valid sequences after cleaning."

# Processing error
"Alignment failed: The sequences are too different to align reliably. 
Try using sequences that are more closely related."
```

## Testing Strategy

### Unit Testing

Unit tests will verify specific examples and edge cases:

- **Parser tests**: Test each supported format with known examples
- **Validator tests**: Test cleaning with specific invalid characters
- **Edge cases**: Empty files, single sequence, identical sequences
- **Error conditions**: Malformed input, invalid characters, insufficient data

### Property-Based Testing

Property-based tests will verify universal properties across randomly generated inputs using the **Hypothesis** library for Python:

- **Configuration**: Each property test should run a minimum of 100 iterations
- **Tagging**: Each test must include a comment referencing the design property
- **Format**: `# Feature: phylogenetics-app, Property {N}: {property description}`
- **Coverage**: Each correctness property listed above must have exactly one corresponding property-based test

**Property Test Examples**:

1. **Parsing property**: Generate random valid FASTA files with varying numbers of sequences and verify all are extracted
2. **Cleaning property**: Generate sequences with random invalid characters and verify only valid characters remain
3. **Alignment property**: Generate random sequences and verify alignment output has correct count and IDs
4. **Tree property**: Generate random alignments and verify tree has correct number of leaves

**Hypothesis Strategies**:

- Use `st.text()` with character whitelists for generating valid sequences
- Use `st.lists()` for generating multiple sequences
- Use `st.integers()` for sequence lengths
- Implement custom strategies for FASTA format generation
- Use `assume()` to filter out invalid test cases

### Integration Testing

Integration tests will verify the complete workflow:

- Upload → Parse → Clean → Align → Build Tree → Export
- Test with real biological sequence examples
- Verify all output files are generated correctly
- Test with different settings combinations (layouts, themes, methods)

### Testing Workflow

1. Implement core functionality first
2. Write property-based tests immediately after implementing each component
3. Run tests frequently during development
4. Fix any failing tests before moving to next component
5. Use unit tests to debug specific failures found by property tests

## Performance Considerations

### Computational Complexity

- **Parsing**: O(n) where n is file size
- **Alignment**: O(m²·l) where m is number of sequences, l is sequence length
- **Tree construction**: O(m³) for Neighbor-Joining
- **Visualization**: O(m) for tree rendering

### Optimization Strategies

1. **Caching**: Use `@st.cache_data` for:
   - Parsed sequences (keyed by file content hash)
   - Alignment results (keyed by sequence set)
   - Tree construction (keyed by alignment)

2. **Lazy Loading**: Only compute visualizations when requested

3. **Limits**: Impose reasonable limits for Streamlit Cloud:
   - Maximum 100 sequences
   - Maximum 10,000 characters per sequence
   - Display warning if approaching limits

4. **Memory Management**:
   - Clear large objects from session state when no longer needed
   - Use generators where possible for file processing
   - Avoid storing multiple copies of sequence data

### Scalability Constraints

The application is designed for educational and small-scale research use:
- Optimal: 5-20 sequences, 500-5000 characters each
- Maximum: 100 sequences, 10,000 characters each
- Not suitable for: Genome-scale analysis, thousands of sequences

## Security Considerations

1. **File Upload**: Limit file size to 10MB
2. **Input Validation**: Sanitize all user input before processing
3. **Resource Limits**: Prevent infinite loops or excessive memory use
4. **No Code Execution**: Never execute user-provided code or commands
5. **Session Isolation**: Ensure user data doesn't leak between sessions (handled by Streamlit)

## Future Enhancements

Potential features for future versions:
- Support for additional alignment algorithms
- Maximum likelihood tree methods
- Bootstrap support values
- Tree editing capabilities
- Sequence annotation display
- Export to additional formats (SVG, PDF)
- Batch processing of multiple files
- Comparison of multiple trees
