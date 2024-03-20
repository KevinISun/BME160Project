
import os
import sys
import seaborn as sns
import pysam
import matplotlib.pyplot as plt

def calculate_depth(bam_file):
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Initialize variables
    depth = {}
    cnt = 0
    # Iterate over each read in the BAM file
    for read in bam.fetch():
        # Get the start and end positions of the read alignment
        start_pos = read.reference_start
        end_pos = read.reference_end

        loop_count = 0
        # Increment the depth count for each position covered by the read
        for pos in range(start_pos, end_pos):
            depth[pos] = depth.get(pos, 0) + 1
            loop_count += 1
            if loop_count % 1000000 == 0:
                print('lc', loop_count)
        cnt += 1
        if cnt % 1000000 == 0:
            print(cnt)

    # Close BAM file
    bam.close()

    return depth

# def calculate_depth(bam_file):
#     # Open BAM file
#     bam = pysam.AlignmentFile(bam_file, "rb")

#     # Initialize variables
#     depth = {}
#     for pileupcolumn in bam.pileup():
#         depth[pileupcolumn.pos] = pileupcolumn.nsegments

#     # Close BAM file
#     bam.close()

#     # Create depth array covering all positions in the genome
#     max_pos = max(depth.keys())
#     depth_array = [depth[pos] if pos in depth else 0 for pos in range(max_pos + 1)]

#     return depth_array

def calculate_genome_coverage(depth_array):
    # Calculate genome coverage
    genome_coverage = sum(depth_array) / len(depth_array)
    return genome_coverage

# def plot_depth(depth):
#     # Plotting
#     positions = list(depth.keys())
#     depths = list(depth.values())
    
#     plt.figure(figsize=(10, 5))
#     plt.plot(positions, depths, color='blue')
#     plt.xlabel('Position')
#     plt.ylabel('Depth')
#     plt.title('Depth Distribution')
#     plt.grid(True)
#     plt.show()
def plot_depth_new(data, output_name, plot_title, genome_size, normalize=False, depth_cut_off=20):
    """Plot genome Depth across genome.
 
    Args:
        depth_report (str): Path to samtool's depth file.
        output_name (str): Path to output PNG image.
        plot_title (str): Plot title.
        genome_size (int): Genome size.
        normalize (bool): If `True`, normalizes the depth by the largest depth (default = `False`).
        depth_cut_off (int): Plot a line to represent a targeted depth (default = 20).
 
    """
 
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
 
    print("Done :)")

def main():
    bam_file = "alignment_1.bam"
    depth = calculate_depth(bam_file)
    plot_depth_new(depth, 'S.png', 'SSRgenome_1', 30000, False, 20)

if __name__ == "__main__":
    main()