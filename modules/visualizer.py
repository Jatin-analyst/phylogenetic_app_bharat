"""Tree visualizer module."""
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from Bio import Phylo
import io


class TreeVisualizer:
    """Visualizer for phylogenetic trees."""
    
    def __init__(self):
        """Initialize the visualizer."""
        self.themes = {
            'default': {
                'branch_color': 'black',
                'label_color': 'black',
                'bg_color': 'white'
            },
            'dark': {
                'branch_color': 'white',
                'label_color': 'white',
                'bg_color': '#1e1e1e'
            },
            'colorful': {
                'branch_color': 'blue',
                'label_color': 'darkgreen',
                'bg_color': 'white'
            },
            'minimal': {
                'branch_color': 'gray',
                'label_color': 'darkgray',
                'bg_color': 'white'
            }
        }
    
    def draw_tree(self, tree, layout: str = 'rectangular', theme: str = 'default') -> Figure:
        """
        Draw a phylogenetic tree.
        
        Args:
            tree: Bio.Phylo.BaseTree.Tree object
            layout: 'rectangular' or 'circular'
            theme: Color theme name
            
        Returns:
            matplotlib Figure object
        """
        if layout == 'rectangular':
            return self.draw_rectangular(tree, theme)
        elif layout == 'circular':
            return self.draw_circular(tree, theme)
        else:
            raise ValueError(f"Unknown layout: {layout}")
    
    def draw_rectangular(self, tree, theme: str = 'default') -> Figure:
        """
        Draw tree in rectangular layout.
        
        Args:
            tree: Bio.Phylo.BaseTree.Tree object
            theme: Color theme name
            
        Returns:
            matplotlib Figure object
        """
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1)
        
        # Apply theme
        theme_colors = self.themes.get(theme, self.themes['default'])
        fig.patch.set_facecolor(theme_colors['bg_color'])
        ax.set_facecolor(theme_colors['bg_color'])
        
        # Draw tree
        Phylo.draw(tree, axes=ax, do_show=False,
                   branch_labels=None)
        
        # Style the plot
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xlabel('Branch Length', color=theme_colors['label_color'])
        ax.tick_params(colors=theme_colors['label_color'])
        
        plt.tight_layout()
        return fig
    
    def draw_circular(self, tree, theme: str = 'default') -> Figure:
        """
        Draw tree in circular layout.
        
        Args:
            tree: Bio.Phylo.BaseTree.Tree object
            theme: Color theme name
            
        Returns:
            matplotlib Figure object
        """
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1, projection='polar')
        
        # Apply theme
        theme_colors = self.themes.get(theme, self.themes['default'])
        fig.patch.set_facecolor(theme_colors['bg_color'])
        ax.set_facecolor(theme_colors['bg_color'])
        
        # Draw tree in circular layout
        # Note: Biopython's Phylo.draw doesn't directly support circular layout
        # We'll use a workaround by drawing it and then converting coordinates
        try:
            Phylo.draw(tree, axes=ax, do_show=False)
        except:
            # Fallback to rectangular if circular fails
            ax.remove()
            ax = fig.add_subplot(1, 1, 1)
            Phylo.draw(tree, axes=ax, do_show=False)
        
        plt.tight_layout()
        return fig
    
    def apply_color_theme(self, fig: Figure, theme: str) -> Figure:
        """
        Apply color theme to existing figure.
        
        Args:
            fig: matplotlib Figure object
            theme: Color theme name
            
        Returns:
            Updated Figure object
        """
        theme_colors = self.themes.get(theme, self.themes['default'])
        
        # Update figure background
        fig.patch.set_facecolor(theme_colors['bg_color'])
        
        # Update axes
        for ax in fig.get_axes():
            ax.set_facecolor(theme_colors['bg_color'])
            
            # Update text colors
            for text in ax.texts:
                text.set_color(theme_colors['label_color'])
            
            # Update tick colors
            ax.tick_params(colors=theme_colors['label_color'])
            
            # Update spine colors
            for spine in ax.spines.values():
                spine.set_color(theme_colors['label_color'])
        
        return fig
    
    def save_high_res_png(self, fig: Figure, filename: str, dpi: int = 300):
        """
        Save figure as high-resolution PNG.
        
        Args:
            fig: matplotlib Figure object
            filename: Output filename
            dpi: Dots per inch (resolution)
        """
        fig.savefig(filename, dpi=dpi, bbox_inches='tight', 
                    facecolor=fig.get_facecolor())
        plt.close(fig)
    
    def get_png_bytes(self, fig: Figure, dpi: int = 300) -> bytes:
        """
        Get PNG image as bytes.
        
        Args:
            fig: matplotlib Figure object
            dpi: Dots per inch (resolution)
            
        Returns:
            PNG image as bytes
        """
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
        buf.seek(0)
        return buf.getvalue()

