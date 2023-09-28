# Graphite V2
### _Longest exact matches, LMEMS, in a collection of sequences_

_Graphite_ starts with two graph files (from Cuttlefish) and a set of query identifiers. It then builds a suffix array of the queries along with other datastructures to speed up matching. Then each sequence (i.e "reference") is read from the graph file and mapped onto the Suffix array. Each mapping is an identical sequence between the query and ref, also called Maximum Exact Matches (MEMs). Each time a MEM is found its length is compared to previously discovered MEMs. When the current MEM covers a region extending beyond a previous, or multiple previous MEMs _Graphite_ replaces it. 

#### Usage
Clone the repo 
`git clone https://github.com/rickbeeloo/GraphiteV2` 
This already includes `libasais` but can then only be used from within the `GraphiteV2` folder. 

To call `Graphite`, use:

`julia main.jl -g test_data/phage_graph.cf_seq -s test_data/phage_graph.cf_seg -k 31 -o test_data/output.txt -q test_data/query.txt`

- `-g`,  the graph path file, when constructed using Cuttlefish the `.cf_seq` file.
- `-s`, the node sequence file, when constructed using Cuttlefish the `.cf_seq` file 
- `-k`, the k-mer size used to build the graph (i.e. 31 nucleotides)
- `-q`, a file with query identifiers. These should match those in `-g` 
- `-o`, the output file path


#### Output 
The output file has four columns:
1) The query identifier (as in `-q`)
2) The start position in nucleotides on the query sequence 
3) The end position in nucleotides on the query sequence 
4) The original LMEM length*

*The LMEM length can differ from the size of the output match (as calculted from end-start+1). This is because the LMEM can have been replaced by another, longer, MEM, yet we keep track of the original length and output this to get an idea of how the matches were shortened