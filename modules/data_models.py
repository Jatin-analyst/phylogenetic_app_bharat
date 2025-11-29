"""Data models for the phylogenetics app."""
from dataclasses import dataclass, field
from typing import List, Dict


@dataclass
class Sequence:
    """Represents a biological sequence with identifier."""
    id: str
    sequence: str
    description: str = ""
    
    def __len__(self) -> int:
        return len(self.sequence)


@dataclass
class ValidationResult:
    """Result of sequence validation and cleaning."""
    valid_sequences: List[Sequence]
    invalid_sequences: List[Sequence]
    cleaning_summary: Dict[str, int] = field(default_factory=dict)
    sequence_type: str = ""  # "nucleotide" or "protein"


@dataclass
class Settings:
    """User settings for tree visualization."""
    layout: str = "rectangular"  # "rectangular" or "circular"
    theme: str = "default"  # "default", "dark", "colorful", "minimal"
    tree_method: str = "nj"  # "nj" or "upgma"
    dpi: int = 300


@dataclass
class RelativeInfo:
    """Information about closest relative for a sequence."""
    closest_relative: str
    distance: float
    explanation: str
