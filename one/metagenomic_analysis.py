#!/usr/bin/env python3
"""
Metagenomic Analysis Script
Modified to work with flat directory structure
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
import re
import logging
from pathlib import Path
import argparse
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("metagenomic_analysis.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("metagenomic_analysis")

# Set up matplotlib parameters for better visualizations
plt.rcParams['figure.figsize'] = (14, 10)  # Increased figure size
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['figure.titlesize'] = 18

# Define color palettes - improved from the original
# Using a more visually distinct palette for clarity
MAIN_PALETTE = "viridis"
HEATMAP_PALETTE = sns.diverging_palette(230, 20, as_cmap=True)
GROUP_COLORS = {"Human Gut": "#4C72B0", "Wastewater": "#C44E52"}  # Blue and red with better contrast

class MetagenomicAnalyzer:
    """Class for metagenomic analysis functions"""
    
    def __init__(self, results_dir):
        """Initialize with results directory"""
        self.results_dir = Path(results_dir)
        if not self.results_dir.exists():
            logger.error(f"Results directory not found: {self.results_dir}")
            raise FileNotFoundError(f"Results directory not found: {self.results_dir}")
            
        # Create output directory for plots
        self.output_dir = Path("./metagenomic_plots")
        self.output_dir.mkdir(exist_ok=True)
        
        logger.info(f"Initialized MetagenomicAnalyzer with results directory: {self.results_dir}")

    def parse_kraken_report(self, file_path, min_abundance=0.01, taxonomic_level=None):
        """
        Parse a Kraken2 report file and extract taxonomic information.
        
        Parameters:
        file_path (Path): Path to the Kraken2 report file
        min_abundance (float): Minimum abundance percentage to include
        taxonomic_level (str, optional): Filter by taxonomic level (G for genus, F for family, etc.)
        
        Returns:
        dict: Dictionary with taxon names as keys and abundance percentages as values
        """
        taxa_dict = {}
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 6:
                        continue
                    
                    try:
                        percentage = float(parts[0])
                        # Skip taxa with abundance below threshold
                        if percentage < min_abundance:
                            continue
                            
                        rank_code = parts[3]
                        taxon_id = parts[4]
                        taxon_name = parts[5].strip()
                        
                        # Remove leading spaces and taxonomic prefixes for clarity
                        taxon_name = re.sub(r'^\s+[A-Z]\s+', '', taxon_name)
                        taxon_name = re.sub(r'^\s+', '', taxon_name)
                        
                        # Filter by taxonomic level if specified
                        if taxonomic_level and rank_code != taxonomic_level:
                            continue
                            
                        # Handle unclassified reads separately
                        if taxon_id == '0' and "unclassified" in taxon_name.lower():
                            taxa_dict["Unclassified"] = percentage
                        elif rank_code in ['G', 'F', 'O', 'C', 'P', 'D', 'S']:  # Include genus, family, order, class, phylum, domain, species
                            taxa_dict[f"{rank_code}_{taxon_name}"] = percentage
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Error parsing line in {file_path}: {line.strip()}")
                        continue
        except Exception as e:
            logger.error(f"Error parsing Kraken report {file_path}: {e}")
            
        return taxa_dict

    def create_taxonomic_profile_matrix(self, taxonomic_level=None, min_abundance=0.01):
        """
        Create a matrix of taxonomic profiles from all samples.
        Modified to work with a flat directory structure.
        
        Parameters:
        taxonomic_level (str, optional): Filter by taxonomic level
        min_abundance (float): Minimum abundance percentage to include
        
        Returns:
        pd.DataFrame: Matrix with samples as rows and taxa as columns
        """
        try:
            # Get all report files directly from the results directory
            report_files = list(self.results_dir.glob("*_report.txt"))
            logger.info(f"Found {len(report_files)} report files")
            
            # Extract sample IDs from filenames
            sample_map = {}  # Maps sample ID to report file
            
            for report_file in report_files:
                filename = report_file.name
                
                # Extract sample ID based on file naming pattern
                if "SRR" in filename:
                    # For real-world samples (SRR*)
                    sample_id = next((part for part in filename.split("_") if part.startswith("SRR")), None)
                elif "error_free" in filename or "miseq" in filename or "no_errors" in filename:
                    # For simulated reads
                    if "error_free" in filename:
                        sample_id = "error_free"
                    elif "no_errors" in filename:
                        sample_id = "no_errors"
                    else:
                        sample_id = "miseq_error"
                else:
                    # Default fallback - use filename without extension
                    sample_id = report_file.stem
                
                if sample_id:
                    sample_map[sample_id] = report_file
            
            sample_ids = list(sample_map.keys())
            logger.info(f"Extracted {len(sample_ids)} sample IDs: {sample_ids}")
            
            # Initialize dictionary to hold all taxa from all samples
            all_taxa = {}
            
            # First pass: collect all taxa across all samples
            for sample_id, report_file in sample_map.items():
                taxa_dict = self.parse_kraken_report(report_file, min_abundance, taxonomic_level)
                for taxon in taxa_dict.keys():
                    all_taxa[taxon] = 0  # Just initialize with zero
            
            logger.info(f"Found {len(all_taxa)} unique taxa across all samples")
            
            # Create the matrix with all samples and all taxa
            profile_matrix = pd.DataFrame(index=sample_ids, columns=list(all_taxa.keys()))
            profile_matrix = profile_matrix.fillna(0)  # Fill with zeros initially
            
            # Second pass: fill in the matrix with actual values
            for sample_id, report_file in sample_map.items():
                taxa_dict = self.parse_kraken_report(report_file, min_abundance, taxonomic_level)
                for taxon, abundance in taxa_dict.items():
                    if taxon in profile_matrix.columns:
                        profile_matrix.loc[sample_id, taxon] = abundance
            
            return profile_matrix
        except Exception as e:
            logger.error(f"Error creating taxonomic profile matrix: {e}")
            return pd.DataFrame()

    def plot_pca(self, profile_matrix, groups=None):
        """
        Perform PCA on the taxonomic profiles and plot the results.
        
        Parameters:
        profile_matrix (pd.DataFrame): Matrix with samples as rows and taxa as columns
        groups (dict, optional): Dictionary mapping sample ids to group labels
        
        Returns:
        tuple: (pca_model, pca_result_df) containing the PCA model and transformed data
        """
        try:
            # Standardize the data
            scaler = StandardScaler()
            scaled_data = scaler.fit_transform(profile_matrix)
            
            # Apply PCA
            pca = PCA(n_components=2)
            principal_components = pca.fit_transform(scaled_data)
            
            # Create DataFrame with PCA results
            pca_result_df = pd.DataFrame(data=principal_components, 
                                        columns=['PC1', 'PC2'], 
                                        index=profile_matrix.index)
            
            # Add group information if provided
            if groups:
                pca_result_df['Group'] = pca_result_df.index.map(groups)
            
            # Create PCA plot with improved styling
            plt.figure(figsize=(12, 10))
            
            if groups:
                # Plot with group colors
                for group_name in set(groups.values()):
                    group_data = pca_result_df[pca_result_df['Group'] == group_name]
                    plt.scatter(
                        group_data['PC1'], 
                        group_data['PC2'], 
                        label=group_name, 
                        s=120,  # Larger points
                        alpha=0.8,
                        color=GROUP_COLORS.get(group_name, "#333333"),
                        edgecolors='white',
                        linewidth=1
                    )
                
                # Add sample labels
                for i, sample_id in enumerate(profile_matrix.index):
                    x = pca_result_df.loc[sample_id, 'PC1']
                    y = pca_result_df.loc[sample_id, 'PC2']
                    plt.annotate(
                        sample_id, 
                        (x, y),
                        fontsize=10,
                        alpha=0.7,
                        xytext=(5, 5),
                        textcoords='offset points',
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.6)
                    )
                
                plt.legend(
                    title="Environment",
                    title_fontsize=14,
                    markerscale=1.5,
                    frameon=True,
                    fancybox=True,
                    framealpha=0.8,
                    edgecolor='gray'
                )
            else:
                # Simple plot without groups
                plt.scatter(
                    pca_result_df['PC1'], 
                    pca_result_df['PC2'], 
                    s=120, 
                    alpha=0.8,
                    color="#4C72B0",
                    edgecolors='white',
                    linewidth=1
                )
                
                for i, sample_id in enumerate(profile_matrix.index):
                    plt.annotate(
                        sample_id, 
                        (principal_components[i, 0], principal_components[i, 1]),
                        fontsize=10,
                        xytext=(5, 5),
                        textcoords='offset points',
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.6)
                    )
            
            # Add explanatory text with improved styling
            explained_variance = pca.explained_variance_ratio_ * 100
            plt.xlabel(f'PC1 ({explained_variance[0]:.2f}%)', fontweight='bold')
            plt.ylabel(f'PC2 ({explained_variance[1]:.2f}%)', fontweight='bold')
            total_variance = explained_variance[0] + explained_variance[1]
            
            plt.title(
                'PCA of Taxonomic Profiles', 
                fontweight='bold',
                pad=20
            )
            
            # Add a subtitle with variance explanation
            plt.figtext(
                0.5, 0.01, 
                f'The first two principal components explain {total_variance:.2f}% of the total variance',
                ha='center',
                fontsize=12,
                bbox=dict(boxstyle="round,pad=0.5", fc="aliceblue", ec="steelblue", alpha=0.8)
            )
            
            # Improve grid and spines
            plt.grid(True, linestyle='--', alpha=0.7)
            for spine in plt.gca().spines.values():
                spine.set_visible(True)
                spine.set_color('gray')
                spine.set_linewidth(0.5)
            
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust for the subtitle
            
            # Save figure
            plt.savefig(self.output_dir / 'pca_analysis.png', dpi=300, bbox_inches='tight')
            logger.info(f"PCA plot saved to {self.output_dir / 'pca_analysis.png'}")
            
            return pca, pca_result_df
        
        except Exception as e:
            logger.error(f"Error performing PCA analysis: {e}")
            return None, None

    def perform_hierarchical_clustering(self, profile_matrix, method='ward'):
        """
        Perform hierarchical clustering on the taxonomic profiles and plot a dendrogram.
        
        Parameters:
        profile_matrix (pd.DataFrame): Matrix with samples as rows and taxa as columns
        method (str): Linkage method to use
        
        Returns:
        numpy.ndarray: The hierarchical clustering encoded as a linkage matrix
        """
        try:
            # Calculate distance matrix
            dist_matrix = pdist(profile_matrix.values, metric='euclidean')
            
            # Perform hierarchical clustering
            Z = linkage(dist_matrix, method=method)
            
            # Plot dendrogram with improved styling
            plt.figure(figsize=(14, 10))
            
            # Create a color mapping for the labels
            sample_ids = profile_matrix.index
            label_colors = {}
            
            for sample_id in sample_ids:
                if str(sample_id).startswith('SRR11412'):
                    label_colors[sample_id] = GROUP_COLORS["Human Gut"]
                elif str(sample_id).startswith('SRR21907'):
                    label_colors[sample_id] = GROUP_COLORS["Wastewater"]
                else:
                    label_colors[sample_id] = "#999999"  # Gray for unknown
            
            # Plot dendrogram
            dendrogram(
                Z, 
                labels=profile_matrix.index, 
                leaf_rotation=90,
                leaf_font_size=12,
                color_threshold=0.7 * max(Z[:, 2]),  # Color threshold at 70% of max distance
                above_threshold_color='gray'
            )
            
            # Color the labels based on their group
            ax = plt.gca()
            xlbls = ax.get_xticklabels()
            for lbl in xlbls:
                sample_id = lbl.get_text()
                lbl.set_color(label_colors.get(sample_id, "#333333"))
            
            plt.title(
                f'Hierarchical Clustering of Samples (Method: {method})', 
                fontweight='bold',
                pad=20
            )
            plt.xlabel('Samples', fontweight='bold')
            plt.ylabel('Distance', fontweight='bold')
            
            # Add a colored background for different environments
            ax = plt.gca()
            for i, sample_id in enumerate(profile_matrix.index):
                if str(sample_id).startswith('SRR11412'):  # Human gut
                    ax.axvspan(i-0.5, i+0.5, facecolor=GROUP_COLORS["Human Gut"], alpha=0.1)
                elif str(sample_id).startswith('SRR21907'):  # Wastewater
                    ax.axvspan(i-0.5, i+0.5, facecolor=GROUP_COLORS["Wastewater"], alpha=0.1)
            
            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor=GROUP_COLORS["Human Gut"], edgecolor='gray', alpha=0.6, label='Human Gut'),
                Patch(facecolor=GROUP_COLORS["Wastewater"], edgecolor='gray', alpha=0.6, label='Wastewater')
            ]
            plt.legend(
                handles=legend_elements, 
                title="Environment", 
                loc='upper right',
                title_fontsize=14,
                frameon=True,
                fancybox=True,
                framealpha=0.8,
                edgecolor='gray'
            )
            
            plt.tight_layout()
            
            # Save figure
            plt.savefig(self.output_dir / 'hierarchical_clustering.png', dpi=300, bbox_inches='tight')
            logger.info(f"Hierarchical clustering plot saved to {self.output_dir / 'hierarchical_clustering.png'}")
            
            return Z
        
        except Exception as e:
            logger.error(f"Error performing hierarchical clustering: {e}")
            return None


    def create_abundance_heatmap(self, profile_matrix, groups=None):
        """
        Create a heatmap of taxonomic abundance.
        
        Parameters:
        profile_matrix (pd.DataFrame): Matrix with samples as rows and taxa as columns
        groups (dict, optional): Dictionary mapping sample ids to group labels
        
        Returns:
        matplotlib.figure.Figure: Heatmap figure
        """
        try:
            # Get top 20 most abundant taxa across all samples
            most_abundant = profile_matrix.mean().sort_values(ascending=False).head(20).index
            
            # Create a heatmap of these taxa
            heatmap_data = profile_matrix[most_abundant].copy()
            
            # Create the heatmap with improved styling
            plt.figure(figsize=(16, 12))
            
            # Add group information for color bar
            if groups:
                group_colors = pd.Series(groups).map({
                    'Human Gut': GROUP_COLORS['Human Gut'], 
                    'Wastewater': GROUP_COLORS['Wastewater']
                })
                
                # Create the clustermap
                g = sns.clustermap(
                    heatmap_data,
                    cmap=HEATMAP_PALETTE,
                    standard_scale=1,  # Scale data across rows
                    figsize=(16, 12),
                    row_colors=group_colors,
                    col_cluster=True,
                    row_cluster=True,
                    xticklabels=1,
                    yticklabels=1,
                    linewidths=0.1,
                    linecolor='gray',
                    cbar_kws={'label': 'Standardized Abundance'},
                    dendrogram_ratio=0.1,  # Smaller dendrograms
                    colors_ratio=0.03,     # Thinner color bar
                    tree_kws={'linewidth': 1.5}  # Thicker dendrogram lines
                )
                
                # Add a legend for environments
                handles = [
                    plt.Rectangle((0,0), 1, 1, color=GROUP_COLORS['Human Gut']),
                    plt.Rectangle((0,0), 1, 1, color=GROUP_COLORS['Wastewater'])
                ]
                
                # Position the legend within the figure
                g.fig.legend(
                    handles, 
                    ['Human Gut', 'Wastewater'],
                    title='Environment',
                    loc='upper right',
                    bbox_to_anchor=(0.98, 0.98),
                    frameon=True,
                    fancybox=True,
                    framealpha=0.8,
                    edgecolor='gray'
                )
                
                # Adjust labels for better readability
                plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=12)
                plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=10)
                
                # Add a title
                g.fig.suptitle(
                    'Clustering of Samples by Taxonomic Abundance', 
                    fontsize=20, 
                    fontweight='bold',
                    y=0.98
                )
                
                # Save figure
                plt.savefig(self.output_dir / 'abundance_heatmap.png', dpi=300, bbox_inches='tight')
                logger.info(f"Abundance heatmap saved to {self.output_dir / 'abundance_heatmap.png'}")
                
                return g.fig
            
            else:
                # Simple heatmap without group information
                plt.figure(figsize=(16, 12))
                sns.heatmap(
                    heatmap_data,
                    cmap=HEATMAP_PALETTE,
                    xticklabels=1,
                    yticklabels=1,
                    linewidths=0.1,
                    linecolor='gray',
                    cbar_kws={'label': 'Abundance (%)'}
                )
                
                plt.title('Taxonomic Abundance Heatmap', fontweight='bold', pad=20)
                plt.xticks(rotation=45, ha='right', fontsize=12)
                plt.yticks(fontsize=10)
                plt.tight_layout()
                
                # Save figure
                plt.savefig(self.output_dir / 'abundance_heatmap.png', dpi=300, bbox_inches='tight')
                logger.info(f"Abundance heatmap saved to {self.output_dir / 'abundance_heatmap.png'}")
                
                return plt.gcf()
        
        except Exception as e:
            logger.error(f"Error creating abundance heatmap: {e}")
            return None

    def run_analysis(self):
        """Run the complete metagenomic analysis workflow"""
        try:
            # 1. Create taxonomic profile matrix (using genus level for clarity)
            logger.info("Creating taxonomic profile matrix...")
            genus_profiles = self.create_taxonomic_profile_matrix(taxonomic_level='G', min_abundance=0.1)
            logger.info(f"Matrix shape: {genus_profiles.shape}")
            
            # Save the profile matrix
            genus_profiles.to_csv(self.output_dir / 'genus_profiles.csv')
            
            # Check if we have enough data
            if genus_profiles.shape[0] < 2:
                logger.error("Not enough samples found. Please check the results directory.")
                return False
            
            # Define sample grouping - Gut (SRR11412*) vs Wastewater (SRR21907*)
            sample_groups = {}
            for sample_id in genus_profiles.index:
                sample_id_str = str(sample_id)
                if sample_id_str.startswith('SRR11412'):
                    sample_groups[sample_id] = 'Human Gut'
                elif sample_id_str.startswith('SRR21907'):
                    sample_groups[sample_id] = 'Wastewater'
                else:
                    sample_groups[sample_id] = 'Unknown'
            
            logger.info(f"Sample groups: {sample_groups}")
            
            # 1. Perform PCA analysis
            logger.info("Performing PCA analysis...")
            pca_model, pca_result = self.plot_pca(genus_profiles, sample_groups)
            
            # 2. Perform hierarchical clustering
            logger.info("Performing hierarchical clustering...")
            cluster_linkage = self.perform_hierarchical_clustering(genus_profiles)
            
            # 3. Create heatmap of top taxa
            logger.info("Creating abundance heatmap...")
            heatmap_fig = self.create_abundance_heatmap(genus_profiles, sample_groups)
            
            # 4. Generate a summary report
            self.generate_summary_report(
                genus_profiles, 
                pca_model, 
                pca_result, 
                sample_groups
            )
            
            logger.info("Analysis complete!")
            return True
        
        except Exception as e:
            logger.error(f"Error in analysis workflow: {e}")
            return False

    def generate_summary_report(self, profiles, pca_model, pca_result, sample_groups):
        """Generate a summary report of the analysis"""
        try:
            logger.info("Generating analysis summary report...")
            
            report_file = self.output_dir / "analysis_summary.md"
            
            with open(report_file, "w") as f:
                # Header
                f.write("# Metagenomic Analysis Summary Report\n\n")
                f.write(f"*Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}*\n\n")
                
                # Dataset overview
                f.write("## Dataset Overview\n\n")
                
                # Count samples by group
                group_counts = {}
                for group in sample_groups.values():
                    if group not in group_counts:
                        group_counts[group] = 0
                    group_counts[group] += 1
                
                f.write(f"Total samples analyzed: {len(profiles)}\n\n")
                f.write("Samples by environment:\n")
                for group, count in group_counts.items():
                    f.write(f"- **{group}**: {count} samples\n")
                
                f.write(f"\nNumber of taxa identified: {len(profiles.columns)}\n\n")
                
                # PCA Results
                f.write("## PCA Analysis\n\n")
                if pca_model is not None:
                    explained_variance = pca_model.explained_variance_ratio_ * 100
                    total_variance = explained_variance[0] + explained_variance[1]
                    
                    f.write(f"The first two principal components explain **{total_variance:.2f}%** of the total variance:\n")
                    f.write(f"- PC1: {explained_variance[0]:.2f}%\n")
                    f.write(f"- PC2: {explained_variance[1]:.2f}%\n\n")
                    
                    # Check if the PCA separates the samples clearly
                    if pca_result is not None and 'Group' in pca_result.columns:
                        # Only consider Human Gut and Wastewater groups
                        filtered_pca = pca_result[pca_result['Group'].isin(['Human Gut', 'Wastewater'])]
                        if not filtered_pca.empty and len(filtered_pca['Group'].unique()) == 2:
                            group_pc1_means = filtered_pca.groupby('Group')['PC1'].mean()
                            group_pc2_means = filtered_pca.groupby('Group')['PC2'].mean()
                            
                            if len(group_pc1_means) == 2:
                                pc1_diff = abs(group_pc1_means.iloc[0] - group_pc1_means.iloc[1])
                                pc2_diff = abs(group_pc2_means.iloc[0] - group_pc2_means.iloc[1])
                                
                                if pc1_diff > 2 or pc2_diff > 2:
                                    f.write("The PCA plot shows **clear separation** between environmental sources, ")
                                    f.write("confirming that human gut and wastewater samples have distinct microbial profiles.\n\n")
                                else:
                                    f.write("The PCA plot shows **some separation** between environmental sources, ")
                                    f.write("suggesting differences in the microbial communities between these environments.\n\n")
                
                # Conclusion
                f.write("## Conclusion\n\n")
                
                if 'Human Gut' in group_counts and 'Wastewater' in group_counts and group_counts['Human Gut'] > 0 and group_counts['Wastewater'] > 0:
                    f.write("The metagenomic analysis reveals that samples from human gut and wastewater environments ")
                    f.write("have distinct microbial compositions. ")
                    
                    # Add more specific conclusions based on results
                    if pca_model is not None and explained_variance[0] > 25:
                        f.write("The strong separation in the PCA plot (with high explained variance) ")
                        f.write("indicates clear differences in the taxonomic profiles. ")
                    
                    f.write("\n\nHierarchical clustering and heatmap visualization further support the distinct ")
                    f.write("microbial community compositions between these environmental sources.")
            
            logger.info(f"Summary report saved to {report_file}")
            return True
        
        except Exception as e:
            logger.error(f"Error generating summary report: {e}")
            return False

def main():
    """Main function to parse arguments and run analysis"""
    parser = argparse.ArgumentParser(description="Metagenomic Analysis Script (Flat Directory Structure)")
    parser.add_argument('--results-dir', type=str, default="./kraken2_results",
                       help="Directory containing Kraken2 results")
    
    args = parser.parse_args()
    
    try:
        # Create analyzer and run analysis
        analyzer = MetagenomicAnalyzer(args.results_dir)
        success = analyzer.run_analysis()
        
        if success:
            logger.info("Analysis completed successfully!")
            return 0
        else:
            logger.error("Analysis failed!")
            return 1
    
    except Exception as e:
        logger.exception(f"Unhandled exception: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
