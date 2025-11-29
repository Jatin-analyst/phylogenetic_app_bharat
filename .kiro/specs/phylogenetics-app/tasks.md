# Implementation Plan

- [x] 1. Set up project structure and dependencies



  - Create project directory structure (modules for parser, validator, aligner, tree_builder, visualizer, exporter)
  - Create requirements.txt with Biopython, Streamlit, Matplotlib, Hypothesis for testing
  - Create main Streamlit app file (app.py)
  - Set up basic Streamlit app skeleton with title and file uploader
  - _Requirements: 10.3, 10.5_

- [x] 2. Implement sequence parser module



  - Create SequenceParser class with methods for each format (FASTA, multi-FASTA, TXT raw, TXT named)
  - Implement format auto-detection logic
  - Implement default ID generation for unnamed sequences
  - Add error handling for invalid formats with descriptive messages
  - _Requirements: 1.1, 1.2, 1.3, 1.4, 1.5_

- [x] 2.1 Write property test for parser




  - **Property 1: Parsing preserves sequence count and identifiers**
  - **Validates: Requirements 1.1, 1.2, 1.3, 1.4**


- [x] 2.2 Write property test for invalid input handling



  - **Property 2: Invalid input rejection**
  - **Validates: Requirements 1.5**

- [x] 3. Implement sequence validator module



  - Create SequenceValidator class with cleaning and validation methods
  - Implement invalid character removal (preserve only valid nucleotide/protein characters)
  - Implement sequence type detection (nucleotide vs protein)
  - Implement minimum length checking and filtering
  - Create ValidationResult dataclass to return cleaned sequences and summary
  - _Requirements: 2.1, 2.2, 2.3, 2.5_

- [x] 3.1 Write property test for sequence cleaning



  - **Property 3: Sequence cleaning preserves valid characters**

  - **Validates: Requirements 2.1**

- [x] 3.2 Write property test for sequence type detection

  - **Property 4: Sequence type detection consistency**
  - **Validates: Requirements 2.2**

- [x] 3.3 Write property test for short sequence filtering


  - **Property 5: Short sequence filtering**
  - **Validates: Requirements 2.3**

- [x] 4. Implement sequence aligner module



  - Create SequenceAligner class using Biopython alignment tools
  - Implement progressive alignment for multiple sequences using pairwise alignments
  - Handle alignment of both nucleotide and protein sequences
  - Add error handling for alignment failures
  - Return Bio.Align.MultipleSeqAlignment object
  - _Requirements: 3.1, 3.4_

- [x] 4.1 Write property test for alignment preservation


  - **Property 6: Alignment preserves sequences and identifiers**
  - **Validates: Requirements 3.1, 3.4**

- [x] 4.2 Write property test for aligned FASTA export


  - **Property 7: Aligned output is valid FASTA**
  - **Validates: Requirements 3.2**

- [x] 5. Implement tree builder module



  - Create TreeBuilder class with distance matrix calculation
  - Implement Neighbor-Joining tree construction method
  - Implement UPGMA tree construction as alternative method
  - Use Biopython's Phylo module for tree objects
  - Add error handling for tree construction failures
  - _Requirements: 4.1, 4.4_

- [x] 5.1 Write property test for tree construction


  - **Property 8: Tree construction produces valid tree**
  - **Validates: Requirements 4.1**






- [ ] 5.2 Write property test for branch lengths
  - **Property 10: Branch lengths are non-negative**
  - **Validates: Requirements 4.4**

- [ ] 6. Implement output exporter module
  - Create OutputExporter class with export methods
  - Implement Newick format export using Biopython
  - Implement aligned FASTA export


  - Implement closest relatives calculation from tree structure
  - Generate JSON output with closest relative information and distances





  - Format JSON for human readability with explanatory text
  - _Requirements: 3.2, 4.2, 7.1, 7.2, 7.5_

- [ ] 6.1 Write property test for Newick round-trip
  - **Property 9: Newick round-trip consistency**
  - **Validates: Requirements 4.2**

- [x] 6.2 Write property test for closest relatives JSON


  - **Property 13: Closest relatives completeness**
  - **Validates: Requirements 7.1, 7.2, 7.5**






- [ ] 7. Implement tree visualizer module
  - Create TreeVisualizer class with Matplotlib-based rendering
  - Implement rectangular tree layout using Biopython's Phylo.draw
  - Implement circular tree layout
  - Implement color theme system (default, dark, colorful, minimal)


  - Add method to save high-resolution PNG (300 DPI minimum)
  - Ensure labels are readable and properly positioned
  - _Requirements: 5.1, 5.2, 5.3, 5.4, 5.5, 6.1, 6.2, 6.3_

- [ ] 7.1 Write property test for layout changes
  - **Property 11: Layout changes preserve tree structure**


  - **Validates: Requirements 5.4**

- [ ] 7.2 Write property test for PNG generation
  - **Property 12: PNG generation with minimum resolution**
  - **Validates: Requirements 6.1, 6.2**




- [ ] 8. Integrate modules into Streamlit UI
  - Create main app layout with file upload section
  - Add settings sidebar (layout selection, color theme, tree method)
  - Implement "Build Tree" button with workflow orchestration
  - Add progress indicators for each processing step (parsing, cleaning, alignment, tree building)


  - Display validation summary after cleaning

  - Wire together: upload → parse → validate → align → build tree → visualize
  - _Requirements: 9.1, 9.2, 9.3, 9.4_

- [ ] 9. Implement results display and downloads
  - Display interactive tree visualization in main area
  - Show closest relatives information in readable format
  - Add download buttons for all outputs (PNG, Newick, aligned FASTA, JSON)
  - Implement proper file naming for downloads
  - Display summary statistics (number of sequences, alignment length, tree method used)
  - _Requirements: 7.3, 7.4, 8.1, 8.2, 9.4_

- [ ] 10. Add error handling and user feedback
  - Implement error catching for all processing steps
  - Display beginner-friendly error messages using st.error()
  - Add helpful suggestions for common errors
  - Prevent processing with insufficient sequences (< 3)
  - Handle edge cases (empty files, single sequence, all invalid sequences)
  - _Requirements: 1.5, 2.5, 3.3, 4.3, 9.5_

- [ ] 11. Implement caching and performance optimizations
  - Add @st.cache_data decorators to expensive operations (parsing, alignment, tree building)
  - Implement session state management for intermediate results
  - Add input size limits (max 100 sequences, max 10KB per sequence)
  - Display warnings when approaching resource limits
  - Clear cached data appropriately when new files are uploaded
  - _Requirements: 10.2_

- [ ] 12. Checkpoint - Ensure all tests pass
  - Ensure all tests pass, ask the user if questions arise.

- [x] 13. Create integration tests

  - Test complete workflow with sample biological sequences
  - Test with different file formats (FASTA, multi-FASTA, TXT)
  - Test with different settings combinations
  - Verify all output files are generated correctly
  - Test error scenarios end-to-end

- [x] 14. Prepare for deployment


  - Create requirements.txt with pinned versions
  - Add README.md with usage instructions
  - Create sample input files for testing
  - Add configuration for Streamlit Cloud (config.toml if needed)
  - Test locally with `streamlit run app.py`
  - _Requirements: 10.1, 10.3, 10.4_

- [x] 15. Final checkpoint - Ensure all tests pass



  - Ensure all tests pass, ask the user if questions arise.
