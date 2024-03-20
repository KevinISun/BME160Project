#BAM Depth Plotter
##Overview

BAM Depth Plotter is a Python script that generates depth plots from BAM files, visualizing the depth of coverage across the genome. It utilizes the pysam library to parse BAM files and matplotlib to create the plots.
##Features

    Estimate genome length from a BAM file.
    Generate depth file from a BAM file.
    Calculate depth of coverage from a BAM file.
    Parse depth file and extract depth values.
    Plot genome depth across the genome.
    Normalize depth by the largest depth.
    Specify a depth cutoff for visualization.
    Command-line interface for easy execution.

##Usage
Run the script with the following command:
python main.py -i input.bam -o output.png

Replace input.bam with the path to your input BAM file and output.png with the desired output PNG file.
Optional Arguments

    -i, --input: Path to the input BAM file (required).
    -o, --output: Path to the output PNG file (required).
    -n, --normalize: Normalize the depth by the largest depth (default=False).
    -c, --cutoff: Plot a line to represent a targeted depth (default=20).

Example with optional arguments:

python main.py -i input.bam -o output.png -n -c 30

