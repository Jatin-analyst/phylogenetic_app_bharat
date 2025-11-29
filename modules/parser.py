"""Sequence parser module."""
from typing import List
from modules.data_models import Sequence


class ParseError(Exception):
    """Exception raised for parsing errors."""
    pass


class SequenceParser:
    """Parser for multiple sequence file formats."""
    
    def parse_file(self, file_content: str, file_type: str = None) -> List[Sequence]:
        """
        Parse a file and extract sequences.
        
        Args:
            file_content: The content of the file as a string
            file_type: Optional file type hint ('fasta', 'txt', or None for auto-detect)
            
        Returns:
            List of Sequence objects
            
        Raises:
            ParseError: If the file format is invalid or cannot be parsed
        """
        if not file_content or not file_content.strip():
            raise ParseError("File is empty or contains no valid data.")
        
        # Auto-detect format if not specified
        if file_type is None:
            file_type = self.detect_format(file_content)
        
        # Parse based on detected format
        if file_type == 'fasta':
            return self.parse_fasta(file_content)
        elif file_type == 'txt':
            # Try to determine if it's named or raw
            if self._is_named_txt(file_content):
                return self.parse_txt_named(file_content)
            else:
                return self.parse_txt_raw(file_content)
        else:
            raise ParseError(
                "Unable to parse file: The uploaded file doesn't appear to be in FASTA or TXT format. "
                "Please check that your file contains sequence data in one of the supported formats."
            )
    
    def detect_format(self, file_content: str) -> str:
        """
        Auto-detect the file format.
        
        Args:
            file_content: The content of the file
            
        Returns:
            'fasta' or 'txt'
        """
        lines = file_content.strip().split('\n')
        
        # Check for FASTA format (starts with >)
        for line in lines:
            line = line.strip()
            if line:  # Skip empty lines
                if line.startswith('>'):
                    return 'fasta'
                break
        
        # Default to TXT format
        return 'txt'
    
    def parse_fasta(self, content: str) -> List[Sequence]:
        """
        Parse FASTA format (single or multi-FASTA).
        
        Args:
            content: FASTA file content
            
        Returns:
            List of Sequence objects
        """
        sequences = []
        current_id = None
        current_seq = []
        current_desc = ""
        
        lines = content.strip().split('\n')
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id is not None:
                    seq_data = ''.join(current_seq)
                    if seq_data:  # Only add if sequence is not empty
                        sequences.append(Sequence(
                            id=current_id,
                            sequence=seq_data,
                            description=current_desc
                        ))
                
                # Parse header line
                header = line[1:].strip()
                if ' ' in header:
                    current_id = header.split()[0]
                    current_desc = header[len(current_id):].strip()
                else:
                    current_id = header
                    current_desc = ""
                current_seq = []
            else:
                # Sequence line
                current_seq.append(line)
        
        # Save last sequence
        if current_id is not None:
            seq_data = ''.join(current_seq)
            if seq_data:  # Only add if sequence is not empty
                sequences.append(Sequence(
                    id=current_id,
                    sequence=seq_data,
                    description=current_desc
                ))
        
        if not sequences:
            raise ParseError("No valid sequences found in FASTA file.")
        
        return sequences
    
    def parse_multi_fasta(self, content: str) -> List[Sequence]:
        """
        Parse multi-FASTA format (same as parse_fasta).
        
        Args:
            content: Multi-FASTA file content
            
        Returns:
            List of Sequence objects
        """
        return self.parse_fasta(content)
    
    def parse_txt_raw(self, content: str) -> List[Sequence]:
        """
        Parse TXT file with raw sequences (one per line).
        
        Args:
            content: TXT file content with raw sequences
            
        Returns:
            List of Sequence objects with auto-generated IDs
        """
        sequences = []
        lines = content.strip().split('\n')
        
        seq_num = 1
        for line in lines:
            line = line.strip()
            if line:  # Skip empty lines
                sequences.append(Sequence(
                    id=f"Seq_{seq_num}",
                    sequence=line,
                    description=""
                ))
                seq_num += 1
        
        if not sequences:
            raise ParseError("No valid sequences found in TXT file.")
        
        return sequences
    
    def parse_txt_named(self, content: str) -> List[Sequence]:
        """
        Parse TXT file with named sequences (name:sequence or name sequence format).
        
        Args:
            content: TXT file content with named sequences
            
        Returns:
            List of Sequence objects
        """
        sequences = []
        lines = content.strip().split('\n')
        
        seq_num = 1
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            # Try colon separator first
            if ':' in line:
                parts = line.split(':', 1)
                seq_id = parts[0].strip()
                seq_data = parts[1].strip()
            # Try space/tab separator
            elif ' ' in line or '\t' in line:
                parts = line.split(None, 1)
                seq_id = parts[0].strip()
                seq_data = parts[1].strip() if len(parts) > 1 else ""
            else:
                # No separator found, treat as raw sequence
                seq_id = f"Seq_{seq_num}"
                seq_data = line
            
            if seq_data:  # Only add if we have sequence data
                sequences.append(Sequence(
                    id=seq_id,
                    sequence=seq_data,
                    description=""
                ))
                seq_num += 1
        
        if not sequences:
            raise ParseError("No valid sequences found in TXT file.")
        
        return sequences
    
    def _is_named_txt(self, content: str) -> bool:
        """
        Determine if TXT content has named sequences.
        
        Args:
            content: TXT file content
            
        Returns:
            True if content appears to have named sequences
        """
        lines = content.strip().split('\n')
        
        # Check first few non-empty lines
        check_count = 0
        named_count = 0
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            # Check if line has name separator (: or space/tab)
            if ':' in line or (' ' in line and not line.startswith(' ')):
                named_count += 1
            
            check_count += 1
            if check_count >= 3:  # Check first 3 lines
                break
        
        # If majority of checked lines have separators, consider it named
        return named_count > check_count / 2 if check_count > 0 else False
