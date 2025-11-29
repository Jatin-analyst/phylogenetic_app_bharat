"""Phylogenetics App - Main Streamlit Application"""
import streamlit as st
import hashlib
from modules.parser import SequenceParser, ParseError
from modules.validator import SequenceValidator
from modules.aligner import SequenceAligner, AlignmentError
from modules.tree_builder import TreeBuilder, TreeBuildError
from modules.visualizer import TreeVisualizer
from modules.exporter import OutputExporter
from modules.data_models import Settings


@st.cache_data
def parse_sequences(file_content: str):
    """Parse sequences with caching."""
    parser = SequenceParser()
    return parser.parse_file(file_content)


@st.cache_data
def validate_sequences(sequences_data, min_length: int = 10):
    """Validate sequences with caching."""
    from modules.data_models import Sequence
    # Reconstruct Sequence objects from data
    sequences = [Sequence(**seq) for seq in sequences_data]
    validator = SequenceValidator(min_length=min_length)
    return validator.validate_sequences(sequences)


@st.cache_data
def align_sequences(sequences_data):
    """Align sequences with caching."""
    from modules.data_models import Sequence
    sequences = [Sequence(**seq) for seq in sequences_data]
    aligner = SequenceAligner()
    return aligner.align_sequences(sequences)


@st.cache_data
def build_tree(alignment_data, method: str = "nj"):
    """Build tree with caching."""
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    
    # Reconstruct alignment
    records = [SeqRecord(Seq(rec['seq']), id=rec['id']) for rec in alignment_data]
    alignment = MultipleSeqAlignment(records)
    
    builder = TreeBuilder()
    return builder.build_tree(alignment, method=method)


def main():
    """Main application entry point."""
    st.set_page_config(
        page_title="Phylogenetics App",
        page_icon="🧬",
        layout="wide"
    )
    
    st.title("🧬 Phylogenetics App")
    st.markdown("""
    Build phylogenetic trees from your sequence data with ease!
    
    **Supported formats:** FASTA, multi-FASTA, TXT (raw or named sequences)
    """)
    
    # Sidebar settings
    st.sidebar.header("⚙️ Settings")
    
    layout = st.sidebar.selectbox(
        "Tree Layout",
        ["rectangular", "circular"],
        help="Choose how to display the phylogenetic tree"
    )
    
    theme = st.sidebar.selectbox(
        "Color Theme",
        ["default", "dark", "colorful", "minimal"],
        help="Choose the color scheme for the tree visualization"
    )
    
    tree_method = st.sidebar.selectbox(
        "Tree Building Method",
        ["nj", "upgma"],
        format_func=lambda x: "Neighbor-Joining (NJ)" if x == "nj" else "UPGMA",
        help="NJ is generally more accurate; UPGMA assumes a molecular clock"
    )
    
    settings = Settings(
        layout=layout,
        theme=theme,
        tree_method=tree_method,
        dpi=300
    )
    
    # File upload section
    st.header("1. Upload Your Sequences")
    uploaded_file = st.file_uploader(
        "Choose a file",
        type=["fasta", "fa", "txt"],
        help="Upload a FASTA or TXT file containing your sequences"
    )
    
    if uploaded_file is not None:
        st.success(f"✅ File uploaded: {uploaded_file.name}")
        
        # Check file size
        file_size = uploaded_file.size
        if file_size > 10 * 1024 * 1024:  # 10MB limit
            st.error("❌ File too large! Please upload a file smaller than 10MB.")
            return
        
        if file_size == 0:
            st.error("❌ File is empty! Please upload a file with sequence data.")
            return
        
        # Read file content
        try:
            file_content = uploaded_file.read().decode('utf-8')
        except UnicodeDecodeError:
            st.error("❌ Unable to read file! Please ensure it's a text file (FASTA or TXT).")
            return
        
        if not file_content.strip():
            st.error("❌ File contains no data! Please upload a file with sequence data.")
            return
        
        # Build Tree button
        st.header("2. Build Phylogenetic Tree")
        
        if st.button("🌳 Build Tree", type="primary", use_container_width=True):
            try:
                # Initialize progress
                progress_bar = st.progress(0)
                status_text = st.empty()
                
                # Step 1: Parse sequences
                status_text.text("📄 Parsing sequences...")
                progress_bar.progress(10)
                
                sequences = parse_sequences(file_content)
                
                # Check sequence count limits
                if len(sequences) > 100:
                    st.warning(f"⚠️ Large dataset detected ({len(sequences)} sequences). Processing may take longer.")
                elif len(sequences) > 50:
                    st.info(f"ℹ️ Processing {len(sequences)} sequences...")
                
                st.success(f"✅ Parsed {len(sequences)} sequences")
                progress_bar.progress(20)
                
                # Check sequence lengths
                max_seq_len = max(len(seq.sequence) for seq in sequences)
                if max_seq_len > 10000:
                    st.warning(f"⚠️ Very long sequences detected (max: {max_seq_len} characters). Processing may take longer.")
                
                # Step 2: Validate and clean sequences
                status_text.text("🧹 Cleaning and validating sequences...")
                progress_bar.progress(30)
                
                validator = SequenceValidator(min_length=10)
                validation_result = validator.validate_sequences(sequences)
                
                # Display validation summary
                if validation_result.cleaning_summary:
                    with st.expander("🔍 Validation Summary"):
                        st.write(f"**Sequence Type:** {validation_result.sequence_type}")
                        st.write(f"**Valid Sequences:** {len(validation_result.valid_sequences)}")
                        st.write(f"**Invalid Sequences:** {len(validation_result.invalid_sequences)}")
                        
                        if validation_result.cleaning_summary:
                            st.write("**Cleaning Actions:**")
                            for action, count in validation_result.cleaning_summary.items():
                                if count > 0:
                                    st.write(f"  - {action.replace('_', ' ').title()}: {count}")
                
                # Check if we have enough valid sequences
                if len(validation_result.valid_sequences) < 3:
                    st.error(
                        f"❌ Not enough sequences: Phylogenetic trees require at least 3 sequences. "
                        f"Your file contains only {len(validation_result.valid_sequences)} valid sequence(s) after cleaning."
                    )
                    return
                
                st.success(f"✅ Validated {len(validation_result.valid_sequences)} sequences")
                progress_bar.progress(40)
                
                # Step 3: Align sequences
                status_text.text("🔄 Aligning sequences...")
                progress_bar.progress(50)
                
                aligner = SequenceAligner()
                alignment = aligner.align_sequences(validation_result.valid_sequences)
                
                st.success(f"✅ Aligned {len(alignment)} sequences")
                progress_bar.progress(60)
                
                # Step 4: Build tree
                status_text.text("🌳 Building phylogenetic tree...")
                progress_bar.progress(70)
                
                builder = TreeBuilder()
                tree = builder.build_tree(alignment, method=settings.tree_method)
                
                st.success(f"✅ Tree built using {settings.tree_method.upper()} method")
                progress_bar.progress(80)
                
                # Step 5: Visualize tree
                status_text.text("🎨 Generating visualization...")
                progress_bar.progress(90)
                
                visualizer = TreeVisualizer()
                fig = visualizer.draw_tree(tree, layout=settings.layout, theme=settings.theme)
                
                progress_bar.progress(100)
                status_text.text("✅ Complete!")
                
                # Store results in session state
                st.session_state['tree'] = tree
                st.session_state['alignment'] = alignment
                st.session_state['fig'] = fig
                st.session_state['settings'] = settings
                st.session_state['sequences'] = validation_result.valid_sequences
                
            except ParseError as e:
                st.error(f"❌ **Parsing Error**")
                st.error(str(e))
                st.info("💡 **Suggestion:** Make sure your file is in FASTA or TXT format with valid sequence data.")
            except AlignmentError as e:
                st.error(f"❌ **Alignment Error**")
                st.error(str(e))
                st.info("💡 **Suggestion:** Try using sequences that are more closely related or check for very short sequences.")
            except TreeBuildError as e:
                st.error(f"❌ **Tree Building Error**")
                st.error(str(e))
                st.info("💡 **Suggestion:** Ensure you have at least 3 valid sequences and they are not identical.")
            except Exception as e:
                st.error(f"❌ **Unexpected Error**")
                st.error(str(e))
                st.info("💡 **Suggestion:** Please check your input file format and try again. If the problem persists, contact support.")
                with st.expander("🔍 Technical Details"):
                    st.exception(e)
        
        # Display results if available
        if 'tree' in st.session_state and 'fig' in st.session_state:
            st.header("3. Results")
            
            # Display tree visualization
            st.subheader("📊 Phylogenetic Tree")
            st.pyplot(st.session_state['fig'])
            
            # Display summary statistics
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Sequences", len(st.session_state['sequences']))
            with col2:
                st.metric("Alignment Length", len(st.session_state['alignment'][0]))
            with col3:
                st.metric("Tree Method", st.session_state['settings'].tree_method.upper())
            
            # Closest Relatives section
            st.subheader("🤝 Closest Relatives")
            
            exporter = OutputExporter()
            relatives_dict = exporter.calculate_closest_relatives(st.session_state['tree'])
            
            # Display in a nice format
            for seq_id, rel_info in relatives_dict.items():
                with st.expander(f"**{seq_id}**"):
                    st.write(rel_info.explanation)
                    st.metric("Genetic Distance", f"{rel_info.distance:.4f}")
            
            # Downloads section
            st.header("4. Downloads")
            st.markdown("Download your results in various formats:")
            
            col1, col2, col3, col4 = st.columns(4)
            
            # PNG Download
            with col1:
                visualizer = TreeVisualizer()
                png_bytes = visualizer.get_png_bytes(
                    st.session_state['fig'], 
                    dpi=st.session_state['settings'].dpi
                )
                st.download_button(
                    label="📷 Download PNG",
                    data=png_bytes,
                    file_name="phylogenetic_tree.png",
                    mime="image/png",
                    use_container_width=True
                )
            
            # Newick Download
            with col2:
                newick_str = exporter.export_newick(st.session_state['tree'])
                st.download_button(
                    label="🌳 Download Newick",
                    data=newick_str,
                    file_name="tree.nwk",
                    mime="text/plain",
                    use_container_width=True
                )
            
            # Aligned FASTA Download
            with col3:
                aligned_fasta = exporter.export_aligned_fasta(st.session_state['alignment'])
                st.download_button(
                    label="📄 Download Aligned FASTA",
                    data=aligned_fasta,
                    file_name="aligned.fasta",
                    mime="text/plain",
                    use_container_width=True
                )
            
            # Closest Relatives JSON Download
            with col4:
                relatives_json = exporter.export_closest_relatives(st.session_state['tree'])
                st.download_button(
                    label="📊 Download Relatives JSON",
                    data=relatives_json,
                    file_name="closest_relatives.json",
                    mime="application/json",
                    use_container_width=True
                )


if __name__ == "__main__":
    main()
