"""Tree builder module."""
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment


class TreeBuildError(Exception):
    """Exception raised for tree building errors."""
    pass


class TreeBuilder:
    """Builder for phylogenetic trees from aligned sequences."""
    
    def __init__(self):
        """Initialize the tree builder."""
        self.distance_calculator = DistanceCalculator('identity')
        self.tree_constructor = DistanceTreeConstructor()
    
    def build_tree(self, alignment: MultipleSeqAlignment, method: str = "nj"):
        """
        Build a phylogenetic tree from aligned sequences.
        
        Args:
            alignment: Bio.Align.MultipleSeqAlignment object
            method: Tree building method ('nj' for Neighbor-Joining or 'upgma' for UPGMA)
            
        Returns:
            Bio.Phylo.BaseTree.Tree object
            
        Raises:
            TreeBuildError: If tree construction fails
        """
        if not alignment or len(alignment) == 0:
            raise TreeBuildError("No sequences provided for tree building")
        
        if len(alignment) < 3:
            raise TreeBuildError(
                "Not enough sequences: Phylogenetic trees require at least 3 sequences. "
                f"Your alignment contains only {len(alignment)} sequence(s)."
            )
        
        try:
            # Calculate distance matrix
            distance_matrix = self.calculate_distance_matrix(alignment)
            
            # Build tree using specified method
            if method.lower() == "nj":
                tree = self.neighbor_joining(distance_matrix)
            elif method.lower() == "upgma":
                tree = self.upgma(distance_matrix)
            else:
                raise TreeBuildError(f"Unknown tree building method: {method}")
            
            return tree
            
        except Exception as e:
            if isinstance(e, TreeBuildError):
                raise
            raise TreeBuildError(
                f"Tree construction failed: {str(e)}. "
                "This may occur if sequences are too similar or too divergent."
            )
    
    def calculate_distance_matrix(self, alignment: MultipleSeqAlignment):
        """
        Calculate pairwise distance matrix from alignment.
        
        Args:
            alignment: Bio.Align.MultipleSeqAlignment object
            
        Returns:
            Bio.Phylo.TreeConstruction.DistanceMatrix object
        """
        try:
            return self.distance_calculator.get_distance(alignment)
        except Exception as e:
            raise TreeBuildError(f"Failed to calculate distance matrix: {str(e)}")
    
    def neighbor_joining(self, distance_matrix):
        """
        Build tree using Neighbor-Joining method.
        
        Args:
            distance_matrix: Distance matrix
            
        Returns:
            Bio.Phylo.BaseTree.Tree object
        """
        try:
            tree = self.tree_constructor.nj(distance_matrix)
            return tree
        except Exception as e:
            raise TreeBuildError(f"Neighbor-Joining failed: {str(e)}")
    
    def upgma(self, distance_matrix):
        """
        Build tree using UPGMA method.
        
        Args:
            distance_matrix: Distance matrix
            
        Returns:
            Bio.Phylo.BaseTree.Tree object
        """
        try:
            tree = self.tree_constructor.upgma(distance_matrix)
            return tree
        except Exception as e:
            raise TreeBuildError(f"UPGMA failed: {str(e)}")

