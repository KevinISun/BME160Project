import os
import pysam
from matplotlib import pyplot as plt
import argparse

def estimate_genome_length(bam_file):
    """
    Estimate the length of the genome from a BAM file.

    Args:
        bam_file (str): Path to the input BAM file.

    Returns:
        int: Estimated length of the genome.
    """
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
    """
    Generate a depth file from a BAM file.

    Args:
        input_bam (str): Path to the input BAM file.
        output_depth (str): Path to the output depth file.
    """
    # Open BAM file
    bam = pysam.AlignmentFile(input_bam, "rb")

    # Generate depth file
    pysam.depth("-a", input_bam, "-o", output_depth)

    # Close BAM file
    bam.close()

def calculate_depth(bam_file):
    """
    Calculate the depth of coverage from a BAM file.

    Args:
        bam_file (str): Path to the input BAM file.

    Returns:
        list: List containing depth of coverage at each position.
    """
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
    """
    Calculate the average genome coverage from a depth array.

    Args:
        depth_array (list): List containing depth of coverage at each position.

    Returns:
        float: Average genome coverage.
    """
    # Calculate genome coverage
    genome_coverage = sum(depth_array) / len(depth_array)
    return genome_coverage

def parse_depth(depth_input, genome_size):
    """
    Parse depth file and extract depth values.

    Args:
        depth_input (str): Path to depth file.
        genome_size (int): Size of the genome.

    Returns:
        list: List containing depth of coverage at each position.
    """
    depth = [0] * genome_size
    references = set()

    with open(depth_input) as depth_object:
        for row in depth_object:
            genome_id, position, depth_count = row.split()
            references.add(genome_id)
            depth[int(position)] = int(depth_count)

    return depth

import os
from matplotlib import pyplot as plt, lines

def plot_depth(depth_file_path, output_file_name, plot_title, genome_size, normalize=False, depth_cut_off=20):
    """
    Plot genome depth across the genome using Matplotlib.

    Args:
        depth_file_path (str): Path to depth file.
        output_file_name (str): Path to output PNG image.
        plot_title (str): Plot title.
        genome_size (int): Size of the genome.
        normalize (bool): If True, normalize the depth by the largest depth (default=False).
        depth_cut_off (int): Plot a line to represent a targeted depth (default=20).
    """
    data = parse_depth(depth_file_path, genome_size)
    y_label = "Depth"

    if normalize:
        plot_title += " (Normalized)"
        y_label = "Normalized Depth"
        data = [xx / max(data) for xx in data]

    plt.figure()
    plt.title(plot_title)
    
    plt.plot(range(len(data)), data)

    plt.xlabel('Genome Position (bp)')
    plt.ylabel(y_label)

    if not normalize:
        plt.axhline(y=depth_cut_off, color="r")

    plt.savefig(output_file_name, dpi=400)
    plt.close()


def generate_depth_graph(bam_file, output_png, normalize=False, depth_cut_off=20):
    """
    Generate depth graph from a BAM file.

    Args:
        bam_file (str): Path to the input BAM file.
        output_png (str): Path to the output PNG file.
    """
    # Estimate genome length
    genome_length = estimate_genome_length(bam_file)
    
    # Generate depth file
    gen_depth_file(bam_file, 'temp_depth_file.depth')
    
    # Calculate depth
    depth = calculate_depth(bam_file)
    
    # Plot depth graph
    plot_depth('temp_depth_file.depth', output_png, f"Genome Depth for {bam_file}" , genome_length, normalize, depth_cut_off)
    
    # Remove temporary depth file
    os.remove('temp_depth_file.depth')

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate depth graph from a BAM file")
    parser.add_argument("-i", "--input", help="Input BAM file", required=True)
    parser.add_argument("-o", "--output", help="Output PNG file", required=True)
    parser.add_argument("-n", "--normalize", action="store_true", help="Normalize the depth")
    parser.add_argument("-c", "--depth_cutoff", type=int, help="Depth cutoff for the red line")
    args = parser.parse_args()
    
    # Call generate_depth_graph function with command-line arguments
    generate_depth_graph(args.input, args.output, args.normalize, args.depth_cutoff)

if __name__ == "__main__":
    main()
