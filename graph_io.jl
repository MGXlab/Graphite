
const MAX_NODE_ID = 2^30 - 1 

function parse_numbers(path::AbstractString, new_id_map::Dict{Int64, Int32})
    numbers = Vector{Int32}()
    for node in split(path, " ")

        # Use Int32 reduced representation
        if isempty(new_id_map)
            number = parse(Int32, view(node, 1:length(node)-1))
            number = node[end] == '-' ? flipnode(number) : noflip(number)
        
        # Use regular IDs (always used unless some node ids exceed int32)
        else
            int64_id = parse(Int64, view(node, 1:length(node)-1))
            number = new_id_map[int64_id]
            number = node[end] == '-' ? flipnode(number) : noflip(number) # Int32 again
        end 

        push!(numbers, number)
        
    end 
    return numbers
end

function read_queries(f::String, query_ids::OrderedSet{String}, node_id_remapping::Dict{Int64, Int32})
    queries = Vector{Vector{Int32}}()
    query_ids_file_order = OrderedSet{String}() # Might find in different order in the file than query_ids
    for (i, line) in enumerate(eachline(f))
        identifier, path = split(line, "\t")
        if identifier in query_ids
            path_numbers = parse_numbers(path, node_id_remapping)
            push!(query_ids_file_order, identifier)
            push!(queries, path_numbers)
            length(query_ids_file_order) == length(query_ids) && break
        end
    end
    return queries, query_ids_file_order
end

function read_node_sizes(f::String)
    size_map = Dict{Int32, Int32}()
    for line in eachline(f)
        node_id, seq = split(line)
        seq_size = Int32(length(seq))
        node_id = parse(Int32, node_id)
        size_map[node_id] = seq_size
    end
    return size_map
end

function read_node_sizes_int32_fit_try(seg_file::String)
    size_map = Dict{Int32, Int32}()
    new_id_map = Dict{Int64, Int32}()
    new_id::Int32 = 0
    old_max_node = 0
    for (i, line) in enumerate(eachline(seg_file))
        new_id < MAX_NODE_ID || error("Could not reduce to fit within int32")
        node_id, seq = split(line)
        seq_size = Int32(length(seq))  # This should still fit though
        old_node_id = parse(Int64, node_id)
        old_max_node = old_node_id > old_max_node ? old_node_id : old_max_node
        size_map[new_id] = seq_size
        new_id_map[old_node_id] = new_id
        new_id += 1
        if mod(i, 1_000_000) == 0 # Print per million nodes
            println("Converted ", i)
        end
    end
    println("Old max node ID: ", old_max_node, " new max ID: ", new_id - 1)
    return new_id_map, size_map
end

function read_ids_from_file(f::String)
    ids = OrderedSet{String}()
    for line in eachline(f)
        push!(ids, strip(line))
    end
    return ids
end

function processGFA(gfa::String, query_file::String, node_id_remapping::Dict{Int64, Int32})
    # Read the query ids from the file 
    query_ids = read_ids_from_file(query_file)
    queries, query_ids = read_queries(gfa, query_ids, node_id_remapping)
    return queries, query_ids
end

function writeResults(ca::Vector{Int32}, color::Color, query_ids::OrderedSet{String}, out_file::String, size_map::Dict{Int32, Int32})
    h = open(out_file, "w+")

    info_vector = color.info_vector

    prev_ori =  info_vector[1]
    aln_start = 1
    genome_loc = 0
    q_count = 1

   
    for i in 1:length(info_vector)-1
  
        if ca[i] > 0 
            node_id = get_original(ca[i])
            node_size = size_map[node_id]
            genome_loc += node_size - color.k_size + 1
        end

        if ca[i+1] < 0 || prev_ori.id != info_vector[i+1].id || prev_ori.pos != info_vector[i+1].pos
            if info_vector[i].len > 0 
                println(h, query_ids[q_count], "\t", info_vector[i].id, "\t", aln_start, "\t", genome_loc+color.k_size-1, "\t", info_vector[i].len)
            end
            prev_ori = info_vector[i+1]
            aln_start = copy(genome_loc) + 1 
        end

        if ca[i+1] < 0
            q_count +=1 
            prev_ori = info_vector[i+1]
            genome_loc = 0
        end 

    end
        

end

