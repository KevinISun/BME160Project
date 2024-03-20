import os
import sys
import seaborn as sns
import pysam
from matplotlib import (pyplot as plt, lines)
import argparse

def estimate_genome_length(bam_file):
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Initialize variables
    max_alignment_position = 0

    # Iterate over aligned reads
    for read in bam.fetch():
        # Get the end position of the read alignment
        end_position = read.reference_end
        if end_position is not None and end_position > max_alignment_position:
            max_alignment_position = end_position

    # Close BAM file
    bam.close()

    # Estimate genome length
    genome_length = max_alignment_position + 1  # Add 1 to account for 0-based indexing

    return genome_length

def gen_depth_file(input_bam, output_depth):
    # Open BAM file
    bam = pysam.AlignmentFile(input_bam, "rb")

    pysam.depth("-a", input_bam, "-o", output_depth)

    # Close BAM file
    bam.close()

def calculate_depth(bam_file):
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Initialize variables
    depth = {}
    for pileupcolumn in bam.pileup():
        depth[pileupcolumn.pos] = pileupcolumn.nsegments

    # Close BAM file
    bam.close()

    # Sort depth dictionary by position
    depth_sorted = sorted(depth.items())

    return [dp for pos, dp in depth_sorted] 

def calculate_genome_coverage(depth_array):
    # Calculate genome coverage
    genome_coverage = sum(depth_array) / len(depth_array)
    return genome_coverage

def parse_depth(depth_input, genome_size):
    """Parse depth file.
 
    Args:
        depth_input (str): Path to depth file.
        genome_size (int): Genome size.
 
    Returns:
        list: List with depth.
 
    """
    depth = [0] * genome_size
    references = set()
 
    with open(depth_input) as depth_object:
        for row in depth_object:
            genome_id, position, depth_count = row.split()
 
            references.add(genome_id)
 
            if len(references) > 1:
                raise Exception(' This script only handles one genome - contig.')
 
            depth[int(position)] = int(depth_count)
 
    return depth

def plot_depth(depth_report, output_name, plot_title, genome_size, normalize=False, depth_cut_off=20):
    """Plot genome Depth across genome.
 
    Args:
        depth_report (str): Path to samtool's depth file.
        output_name (str): Path to output PNG image.
        plot_title (str): Plot title.
        genome_size (int): Genome size.
        normalize (bool): If `True`, normalizes the depth by the largest depth (default = `False`).
        depth_cut_off (int): Plot a line to represent a targeted depth (default = 20).
 
    """
    data = parse_depth(depth_report, genome_size)

    y_label = "Normalized Depth" if normalize else "Depth"
    data = [xx / max(data) for xx in data] if normalize else data
 
    sns.set(color_codes=True)
    plt.title(plot_title)
    ax = plt.subplot(111)
 
    sns_plot = sns.lineplot(x=range(len(data)), y=data)
    sns_plot.set(xlabel='Genome Position (bp)', ylabel=y_label)
 
    if not normalize:
        ax.add_line(lines.Line2D([0, genome_size + 1], [depth_cut_off], color="r"))
 
    plt.savefig(output_name, bbox_inches='tight', dpi=400)
    plt.close()

def generate_depth_graph(bam_file, output_png):
    # Estimate genome length
    genome_length = estimate_genome_length(bam_file)
    
    # Generate depth file
    gen_depth_file(bam_file, 'temp_depth_file.depth')
    
    # Calculate depth
    depth = calculate_depth(bam_file)
    
    # Plot depth graph
    plot_depth('temp_depth_file.depth', output_png, f"Genome Depth for {bam_file}" , genome_length, False, 20)
    
    # Remove temporary depth file
    os.remove('temp_depth_file.depth')

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate depth graph from a BAM file")
    parser.add_argument("-i", "--input", help="Input BAM file", required=True)
    parser.add_argument("-o", "--output", help="Output PNG file", required=True)
    args = parser.parse_args()
    
    # Call generate_depth_graph function with command-line arguments
    generate_depth_graph(args.input, args.output)


if __name__ == "__main__":
    main()