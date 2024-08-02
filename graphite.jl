using Dates
using Base.Threads
const MASK = Int32(1<<30) 

function flipnode(n::Int32)
    n < 2^30 - 1 || error("Encountered node ID that is too large (>= 2^30 - 1)")
    n ⊻ MASK
end

isflipped(n::Int32) = ifelse(n & MASK != 0, true, false)
noflip(n::Int32) = n & (MASK - 1)
convert_node(n::Int32) = isflipped(n) ? noflip(n) : flipnode(n)
get_original(n::Int32) = isflipped(n) ? noflip(n) : n

const NUM_LOCKS = 1024  # You can adjust this number based on your needs
const color_locks = [ReentrantLock() for _ in 1:NUM_LOCKS]


mutable struct ColorInfo
    len::Atomic{UInt32}
    id::Atomic{Int32}
    pos::Atomic{Int32}
    lck::Base.Threads.SpinLock
end

ColorInfo() = ColorInfo(Atomic{UInt32}(0), Atomic{Int32}(-1), Atomic{Int32}(-1), Base.Threads.SpinLock())

struct Color 
    info_vector::Vector{ColorInfo}
    size_map::Dict{Int32, UInt32}
    k_size::Int32
end

function atomic_update_color_info!(color_info::ColorInfo, new_len::UInt32, new_id::Int32, new_pos::Int32)
    while true
        old_len = color_info.len[]
        if new_len <= old_len
            return
        end
        if atomic_cas!(color_info.len, old_len, new_len) == old_len
            atomic_cas!(color_info.id, color_info.id[], new_id)
            atomic_cas!(color_info.pos, color_info.pos[], new_pos)
            return
        end
    end
end

function update_color!(color::Color, ref_id::Int32, match_start::Int32, match_size::Int32, ca::Vector{Int32})    
    match_end = match_start + match_size - Int32(1)

    at_start_info = color.info_vector[match_start]
    at_end_info = color.info_vector[match_end]
    
    if at_start_info.id[] > 0 && at_end_info.id[] == at_start_info.id[] && at_start_info.pos[] == at_end_info.pos[]
        return 
    
    elseif match_start == match_end && at_start_info.len[] > 0
        return 
    end
    
    match_size_nt = 0
    
    for i in match_start:match_end
        match_size_nt += color.size_map[get_original(ca[i])]
    end

    # Lock range within as well 
    match_size_nt = match_size_nt - ((match_size-1) * (color.k_size-1))

    for i in match_start:match_end
        lock(color.info_vector[i].lck)
    end

    for i in match_start:match_end
        atomic_update_color_info!(color.info_vector[i], UInt32(match_size_nt), ref_id, match_start)
    end

    for i in reverse(match_start:match_end) 
        unlock(color.info_vector[i].lck)
    end

    
end

function reverse_complement_ref!(ref::Vector{Int32})
    reverse!(ref)
    @inbounds for i in eachindex(ref)
        ref[i] = convert_node(ref[i])
    end
end


function decide_on_seed(insert_point::Int32, ca::Vector{Int32}, sa::Vector{Int32}, ref::AbstractVector{Int32}, ref_start::Int32)
    left_of_ip = insert_point > 1 ? check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point-Int32(1),  Int32(0)) : 0 
    left_of_ip > 0 && return insert_point-Int32(1), Int32(left_of_ip)

    right_of_ip = insert_point <= length(sa) ? check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point, Int32(0)) : 0
    right_of_ip > 0 && return insert_point, Int32(right_of_ip)

    return Int32(0), Int32(0)
end

function matches_till(ref::AbstractVector{Int32}, ref_start::Int32, ca::Vector{Int32}, q_start::Int32)
    (ref_start > length(ref) || q_start > length(ca)) && return 0
    smallest_n = min(length(ref)-ref_start+1, length(ca)-q_start+1)
    for i in 1:smallest_n
        if ref[ref_start + i - 1] != ca[q_start+i-1]
            return Int32(i - 1)
        end 
    end 
    return Int32(smallest_n)
end

function check_this_point(ca::Vector{Int32}, sa::Vector{Int32}, ref::AbstractVector{Int32}, ref_start::Int32, point::Int32, skip::Int32)
    ca_suffix_index = sa[point]
    ca_start = ca_suffix_index + skip
    ref_start = ref_start + skip
    match_size::Int32 = matches_till(ref, ref_start, ca, ca_start) + skip
    return match_size
end

function extend_from_point!(ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, lcp::Vector{Int32}, point::Int32, forward::Bool, ref_start::Int32, match_size::Int32, ref_id::Int32, color::Color)
    move_dir = forward ? 1 : -1
    lcp_dir  = forward ? 0 :  1
    i = point += move_dir
    while i > 1 && i <= length(sa) && lcp[i + lcp_dir] > 0
        start_check_from = Int32(min(lcp[i + lcp_dir], match_size))
        match_size::Int32 = check_this_point(ca, sa, ref, ref_start, Int32(i), start_check_from)
        update_color!(color, ref_id, sa[i], Int32(match_size), ca)
        i += move_dir        
    end
end

function do_the_work(ref_id::Int32, color::Color, ca::Vector{Int32}, sa::Vector{Int32}, s_idx::Int32, e_idx::Int32, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
    ref_start = s_idx

    while ref_start <= e_idx  
        insert_point = locate_insert_point(sa, ca, view(ref, ref_start:length(ref)))
        max_match_index, max_match_size = decide_on_seed(insert_point, ca, sa, ref, ref_start)

        if max_match_size > 0
            while ref_start <= e_idx
                max_match_size = check_this_point(ca, sa, ref, ref_start, max_match_index, Int32(max_match_size-1))
                
                max_match_size == 0 && break 

                update_color!(color, ref_id, sa[max_match_index], max_match_size, ca)
                            
                extend_from_point!(ca, sa, ref, lcp, max_match_index, false, ref_start, Int32(max_match_size), ref_id, color)
                extend_from_point!(ca, sa, ref, lcp, max_match_index, true, ref_start, Int32(max_match_size), ref_id, color)

                max_match_index = inv_perm_sa[sa[max_match_index]+1]
                ref_start += Int32(1)
            end 
        else 
            ref_start += Int32(1)
        end
    end
end

function run(gfa::String, seq_file::String, query_file::String, k_size::Int32, out_file::String; blacklist::String = "") 
    blacklist_ids = !isempty(blacklist) ? read_ids_from_file(blacklist) : OrderedSet{String}()

    println("[STEP 1] Reading node sizes...")
    node_id_remapping =  Dict{Int64, Int32}()
    size_map = Dict{Int32, Int32}()
    try
        size_map = read_node_sizes(seq_file)
    catch 
        println("Probably overflow of node ids, trying to fit in int32s..")
        node_id_remapping, size_map = read_node_sizes_int32_fit_try(seq_file)
    end

    println("[STEP 2] Reading query paths...")
    queries, query_ids = processGFA(gfa, query_file, node_id_remapping)
    
    println("[STEP 3] Construction SA of queries")
    ca, sa = create_k_suffix_array(queries, Int32(0))

    println("[STEP 4] Construction ISA of queries")
    inv_sa_perm = inverse_perm_sa(sa)

    println("[STEP 5] Construction LCP from SA")
    lcp = build_lcp(sa, ca)

    info_vector = [ColorInfo() for _ in 1:length(ca)]
    color = Color(info_vector, size_map, k_size)

    println("[STEP 6] Starting graph file iteration..")
    
    nthreads = Threads.nthreads()
    
    for (ref_id, line) in enumerate(eachline(gfa))
        identifier, path = split(line, "\t")        
        if !(identifier in query_ids) && !(identifier in blacklist_ids)
            path_numbers = parse_numbers(path, node_id_remapping)
            chunk_size = div(length(path_numbers) + nthreads - 1, nthreads)
            
            # Forward pass
            Threads.@threads for i in 1:nthreads
                start_idx = Int32((i - 1) * chunk_size + 1)
                end_idx = Int32(min(i * chunk_size, length(path_numbers)))
                do_the_work(Int32(ref_id), color, ca, sa, start_idx, end_idx, path_numbers, inv_sa_perm, lcp)
            end

            # Reverse complement pass
            reverse_complement_ref!(path_numbers)

            Threads.@threads for i in 1:nthreads
                start_idx = Int32((i - 1) * chunk_size + 1)
                end_idx = Int32(min(i * chunk_size, length(path_numbers)))
                do_the_work(Int32(ref_id), color, ca, sa, start_idx, end_idx, path_numbers, inv_sa_perm, lcp)
            end

            println("$ref_id - $identifier done!")
        end 
    end

    println("[STEP 7] Writing results")
    writeResults(ca, color, query_ids, out_file, size_map)
end