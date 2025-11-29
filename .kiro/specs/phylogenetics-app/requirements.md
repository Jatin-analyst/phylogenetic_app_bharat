# Requirements Document

## Introduction

The Phylogenetics App is a beginner-friendly web application that enables users to analyze genetic sequences and visualize evolutionary relationships through phylogenetic trees. The system accepts multiple sequence input formats, performs automated alignment and tree construction, and provides high-quality visualizations with explanatory outputs. The application is designed to be accessible to non-biology users while providing professional-grade outputs for researchers.

## Glossary

- **Phylogenetic Tree**: A branching diagram showing the evolutionary relationships between biological sequences
- **FASTA**: A text-based format for representing nucleotide or protein sequences
- **Multi-FASTA**: A FASTA file containing multiple sequences
- **Multi-Sequence Alignment (MSA)**: The process of aligning three or more biological sequences to identify regions of similarity
- **Newick Format**: A standard text format for representing phylogenetic trees using nested parentheses
- **Streamlit App**: The web application interface built using the Streamlit framework
- **Sequence Cleaning**: The process of removing invalid characters and formatting sequences properly
- **Closest-Known-Relative**: The most similar sequence to a given sequence in the phylogenetic tree

## Requirements

### Requirement 1

**User Story:** As a user, I want to upload sequence files in multiple formats, so that I can work with data from different sources without manual conversion.

#### Acceptance Criteria

1. WHEN a user uploads a FASTA file, THEN the Streamlit App SHALL parse and validate the sequence data
2. WHEN a user uploads a multi-FASTA file, THEN the Streamlit App SHALL extract all sequences and their identifiers
3. WHEN a user uploads a TXT file with raw sequences, THEN the Streamlit App SHALL process the sequences and assign default identifiers
4. WHEN a user uploads a TXT file with named sequences, THEN the Streamlit App SHALL parse both sequence names and sequence data
5. WHEN a user uploads an invalid file format, THEN the Streamlit App SHALL display a clear error message and prevent further processing

### Requirement 2

**User Story:** As a user, I want my sequences to be automatically cleaned and validated, so that I can proceed with analysis without worrying about data quality issues.

#### Acceptance Criteria

1. WHEN the Streamlit App receives sequence data, THEN the Streamlit App SHALL remove invalid characters from sequences
2. WHEN the Streamlit App detects sequences with inconsistent types, THEN the Streamlit App SHALL identify whether sequences are nucleotide or protein
3. WHEN the Streamlit App finds sequences that are too short for analysis, THEN the Streamlit App SHALL notify the user and exclude those sequences
4. WHEN sequence cleaning is complete, THEN the Streamlit App SHALL display a summary of cleaning actions performed
5. WHEN all sequences are invalid, THEN the Streamlit App SHALL prevent tree building and display an explanatory error message

### Requirement 3

**User Story:** As a user, I want sequences to be automatically aligned, so that I can build accurate phylogenetic trees without manual alignment work.

#### Acceptance Criteria

1. WHEN a user initiates tree building with multiple sequences, THEN the Streamlit App SHALL perform multi-sequence alignment on the cleaned sequences
2. WHEN alignment is complete, THEN the Streamlit App SHALL generate an aligned FASTA file for download
3. WHEN alignment fails, THEN the Streamlit App SHALL display an error message with guidance on resolving the issue
4. WHEN alignment completes successfully, THEN the Streamlit App SHALL preserve sequence identifiers in the aligned output

### Requirement 4

**User Story:** As a user, I want a phylogenetic tree to be automatically built from my sequences, so that I can visualize evolutionary relationships without manual tree construction.

#### Acceptance Criteria

1. WHEN aligned sequences are available, THEN the Streamlit App SHALL construct a phylogenetic tree using a distance-based or likelihood method
2. WHEN tree construction is complete, THEN the Streamlit App SHALL generate a Newick format file for download
3. WHEN tree construction fails, THEN the Streamlit App SHALL display an error message explaining the failure
4. WHEN the tree is built, THEN the Streamlit App SHALL calculate branch lengths representing evolutionary distances

### Requirement 5

**User Story:** As a user, I want to visualize the phylogenetic tree in multiple styles, so that I can choose the presentation that best suits my needs.

#### Acceptance Criteria

1. WHEN a phylogenetic tree is available, THEN the Streamlit App SHALL display the tree in the user-selected layout style
2. WHERE the user selects rectangular layout, THEN the Streamlit App SHALL render the tree with horizontal branches
3. WHERE the user selects circular layout, THEN the Streamlit App SHALL render the tree in a radial arrangement
4. WHEN the user changes the layout style, THEN the Streamlit App SHALL update the visualization without rebuilding the tree
5. WHEN rendering the tree, THEN the Streamlit App SHALL apply the user-selected color theme to branches and labels

### Requirement 6

**User Story:** As a user, I want to download high-resolution tree images, so that I can use them in presentations and publications.

#### Acceptance Criteria

1. WHEN a tree visualization is displayed, THEN the Streamlit App SHALL generate a high-resolution PNG file
2. WHEN the user downloads the PNG file, THEN the Streamlit App SHALL provide an image with minimum 300 DPI resolution
3. WHEN generating the PNG, THEN the Streamlit App SHALL include all labels and branch information clearly
4. WHEN the user changes visualization settings, THEN the Streamlit App SHALL regenerate the PNG with updated settings

### Requirement 7

**User Story:** As a non-biology user, I want explanations of closest relatives for each sequence, so that I can understand the relationships without specialized knowledge.

#### Acceptance Criteria

1. WHEN a phylogenetic tree is built, THEN the Streamlit App SHALL identify the closest relative for each sequence
2. WHEN closest relatives are identified, THEN the Streamlit App SHALL generate a JSON file containing relationship information
3. WHEN displaying closest relatives, THEN the Streamlit App SHALL present the information in plain language
4. WHEN the user views the closest relatives explanation, THEN the Streamlit App SHALL show which sequences are most similar to each other
5. WHEN generating the JSON file, THEN the Streamlit App SHALL include distance metrics for each relationship

### Requirement 8

**User Story:** As a user, I want to interact with the tree visualization in the web interface, so that I can explore relationships dynamically.

#### Acceptance Criteria

1. WHEN the Streamlit App displays a tree, THEN the Streamlit App SHALL provide an interactive view within the browser
2. WHEN the user hovers over tree nodes, THEN the Streamlit App SHALL display sequence names and distance information
3. WHEN the user applies color themes, THEN the Streamlit App SHALL update the visualization colors immediately
4. WHEN the tree is displayed, THEN the Streamlit App SHALL ensure all labels are readable and non-overlapping

### Requirement 9

**User Story:** As a user, I want a simple interface with minimal steps, so that I can complete my analysis quickly without confusion.

#### Acceptance Criteria

1. WHEN a user opens the Streamlit App, THEN the Streamlit App SHALL display a clear file upload interface
2. WHEN a user uploads a file, THEN the Streamlit App SHALL provide a single "Build Tree" button to initiate analysis
3. WHEN analysis is in progress, THEN the Streamlit App SHALL display progress indicators for each step
4. WHEN analysis is complete, THEN the Streamlit App SHALL present all outputs in a clearly organized layout
5. WHEN errors occur, THEN the Streamlit App SHALL display beginner-friendly error messages with suggested solutions

### Requirement 10

**User Story:** As a user, I want the app to run on Streamlit Cloud without local installations, so that I can access it from anywhere without setup.

#### Acceptance Criteria

1. WHEN the Streamlit App is deployed to Streamlit Cloud, THEN the Streamlit App SHALL run without requiring local software installation
2. WHEN the Streamlit App performs computations, THEN the Streamlit App SHALL complete within Streamlit Cloud resource limits
3. WHEN the Streamlit App is accessed, THEN the Streamlit App SHALL load all dependencies from the requirements file
4. WHEN multiple users access the app, THEN the Streamlit App SHALL handle concurrent sessions independently
5. WHEN the Streamlit App uses external libraries, THEN the Streamlit App SHALL use only lightweight dependencies compatible with Streamlit Cloud
