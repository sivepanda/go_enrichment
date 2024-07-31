# go_enrichment
## Introduction
This project aims to create p-values for the abundance of each GO Categories in each chromosome queried by doing a T-Test of number of hits in a particular chromosome vs those that are NOT on the chromosome. This helps us identify if there is a significant correlation between a particular GO Category and a chromosome.

# Usage
Usage is fairly simple, as the program is designed to be a python library. To run, all you need to do is install the library with `pip install .` in the main directory in the project. This will install dependencies as well. Then, you can run `go_enrichment {output file} {number of hits}`. Be sure to replace the output file and number of hits parameters with your desired values.


*This project was created at the Wren Lab at Oklahoma Medical Research Foundation by Siven Panda.*
