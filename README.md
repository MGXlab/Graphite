# Graphite V2
### _Longest exact matches, LMEMS, in a collection of sequences_

_Graphite_ starts with two graph files (from Cuttlefish) and a set of query identifiers. It then builds a suffix array of the queries along with other datastructures to speed up matching. Then each sequence (i.e "reference") is read from the graph file and mapped onto the Suffix array. Each mapping is an identical sequence between the queries and ref, also called Maximum Exact Matches (MEMs). Each time a MEM is found its length is compared to previously discovered MEMs to only retain the Longest MEM (LMEM). 


#### Usage
Clone the repo 
`git clone https://github.com/rickbeeloo/GraphiteV2` 
This already includes [libasais](https://github.com/IlyaGrebnov/libsais) but can then only be used from within the `GraphiteV2` folder.

To use `graphite`:

`julia main.jl -g test_data/phage_graph.cf_seq -s test_data/phage_graph.cf_seg -k 31 -o test_data/output.txt -q test_data/query.txt`

- `-g`,  the graph path file, when constructed using Cuttlefish the `.cf_seq` file.
- `-s`, the node sequence file, when constructed using Cuttlefish the `.cf_seq` file 
- `-k`, the k-mer size used to build the graph (i.e. 31 nucleotides)
- `-q`, a file with query identifiers. These should match those in `-g` 
- `-o`, the output file path


#### Output 
The output file looks like:
```
Reference:1_Sequence:MT952848.1	2	1	112	112
Reference:1_Sequence:MT952848.1	1	150	184	35
Reference:1_Sequence:MT952848.1	1	192	337	146
Reference:1_Sequence:MT952848.1	1	339	398	60
Reference:1_Sequence:MT952848.1	1	432	491	60
```

Which has five columns:
1) The query identifier (as in `-q`)
2) The reference index (1-based)
2) The start position in nucleotides on the query sequence 
3) The end position in nucleotides on the query sequence 
4) The original MEM length*

So for example, `Reference:1_Sequence:MT952848.1` matches the first 112 nucleotides with the second reference in the file (`test_data/phage_graph.cf_seq`)

*The original MEM size corresponding to the LMEM, to get the LMEM length use end-start+1