__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

from snakemake.shell import shell
from os import path
import os

ref_genome = snakemake.params.get("ref_genome", "")

shell("bioawk -c fastx '{{ print $name, length($seq) }}' < {ref_genome} > {snakemake.output.reference_size_file}")
shell("bedtools makewindows -g {snakemake.output.reference_size_file} -w 1000 > {snakemake.output.reference_window_file}")

