"""Output exporter module."""
import json
from typing import Dict
from io import StringIO
from Bio import Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment
from modules.data_models import RelativeInfo


class OutputExporter:
    """Exporter for various output formats."""
    
    def export_newick(self, tree) -> str:
        """
        Export tree to Newick format.
        
        Args:
            tree: Bio.Phylo.BaseTree.Tree object
            
        Returns:
            Newick format string
        """
        output = StringIO()
        Phylo.write(tree, output, 'newick')
        return output.getvalue()
    
    def export_aligned_fasta(self, alignment: MultipleSeqAlignment) -> str:
        """
        Export alignment to FASTA format.
        
        Args:
            alignment: Bio.Align.MultipleSeqAlignment object
            
        Returns:
            FASTA format string
        """
        output = StringIO()
        AlignIO.write(alignment, output, 'fasta')
        return output.getvalue()
    
    def export_closest_relatives(self, tree) -> str:
        """
        Export closest relatives information as JSON.
        
        Args:
            tree: Bio.Phylo.BaseTree.Tree object
            
        Returns:
            JSON format string
        """
        relatives_dict = self.calculate_closest_relatives(tree)
        
        # Convert to JSON-serializable format
        json_dict = {}
        for seq_id, rel_info in relatives_dict.items():
            json_dict[seq_id] = {
                'closest_relative': rel_info.closest_relative,
                'distance': rel_info.distance,
                'explanation': rel_info.explanation
            }
        
        return json.dumps(json_dict, indent=2)
    
    def calculate_closest_relatives(self, tree) -> Dict[str, RelativeInfo]:
        """
        Calculate closest relative for each sequence in the tree.
        
        Args:
            tree: Bio.Phylo.BaseTree.Tree object
            
        Returns:
            Dictionary mapping sequence IDs to RelativeInfo objects
        """
        relatives = {}
        terminals = tree.get_terminals()
        
        # For each terminal (leaf), find its closest relative
        for terminal in terminals:
            seq_id = terminal.name
            min_distance = float('inf')
            closest_relative = None
            
            # Calculate distance to all other terminals
            for other_terminal in terminals:
                if other_terminal.name != seq_id:
                    # Calculate distance between the two terminals
                    distance = tree.distance(terminal, other_terminal)
                    
                    if distance < min_distance:
                        min_distance = distance
                        closest_relative = other_terminal.name
            
            # Create explanation
            if closest_relative:
                explanation = (
                    f"{seq_id} is most closely related to {closest_relative} "
                    f"with a genetic distance of {min_distance:.4f}"
                )
            else:
                explanation = f"{seq_id} has no close relatives in this tree"
                min_distance = 0.0
                closest_relative = "None"
            
            relatives[seq_id] = RelativeInfo(
                closest_relative=closest_relative,
                distance=min_distance,
                explanation=explanation
            )
        
        return relatives

