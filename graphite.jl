using Dates
const MASK = Int32(1<<30) 

function flipnode(n::Int32)
    n < 2^30 - 1 || error("Encounterd node ID that is too large (>= 2^30 - 1)")
    n ⊻ MASK
end

# color = Color( [ColorInfo(Base.Threads.SpinLock(), 10, 10, 10), ColorInfo(Base.Threads.SpinLock(), 10, 10, 10)], Dict{Int32, UInt32}(), 3  )

isflipped(n::Int32) = ifelse(n & MASK != 0, true, false)
noflip(n::Int32) = n & (MASK - 1)
convert_node(n::Int32) = isflipped(n) ? noflip(n) : flipnode(n)
get_original(n::Int32) = isflipped(n) ? noflip(n) : n

mutable struct ColorInfo
    lck::Base.Threads.SpinLock
    len::UInt32 
    id::Int32
    pos::Int32
end

struct Color 
    info_vector::Vector{ColorInfo}
    size_map::Dict{Int32, UInt32}
    k_size::Int32
end

function update_color!(color::Color, ref_id::Int32, match_start::Int32, match_size::Int32, ca::Vector{Int32})    
    match_end = match_start+match_size-Int32(1)
    
    # Get info from pervious matches
    at_start_info = color.info_vector[match_start]
    at_end_info   = color.info_vector[match_end]

    if at_start_info.id > 0 && at_end_info.id == at_end_info.id && at_start_info.pos == at_end_info.pos
        return 

    # Don't have to bother about single node matches as they can't be longer anyway
    elseif match_start == match_end && at_start_info.len > 0
        return 
    else

        match_size_nt = 0
        for i in match_start:match_end
            match_size_nt += color.size_map[get_original(ca[i])]
        end

        # We have to consider the overlap between k-mers as well 
        match_size_nt = match_size_nt - ((match_size-1) * (color.k_size-1)) 


        for i in match_start:match_end
            color_info = color.info_vector[i]
            if color_info.len < match_size_nt 
                lock(color.info_vector[i].lck)
                color.info_vector[i].len = match_size_nt
                color.info_vector[i].id = ref_id
                color.info_vector[i].pos = match_start
                unlock(color.info_vector[i].lck)
            end
        end

        return 
    end
end

# Note, in place version
function reverse_complement_ref!(ref::Vector{Int32})
    reverse!(ref)
    @inbounds for i in eachindex(ref)
        ref[i] = convert_node(ref[i])
    end
end

function add_reversed_ref(ref::Vector{Int32})
    concat_ref = Vector{Int}(undef, length(ref) * 2)

    # Copy old ref 
    @inbounds @simd for i in eachindex(ref)
        concat_ref[i] = ref
    end

    # Copy reversed 
    ref_l = length(ref)
    @inbounds @simd for i in reverse(0:length(ref)) # reverse iter (from end of ref to beginning)
        concat_ref[ref_l + i] = flipnode(ref[i])
    end

    return concat_ref
end


# Check left and right from insert point to see if we have a match
function decide_on_seed(insert_point::Int32, ca::Vector{Int32}, sa::Vector{Int32}, ref::AbstractVector{Int32}, ref_start::Int32)
    # Check left for a match
    left_of_ip = insert_point > 1 ? check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point-Int32(1),  Int32(0)) : 0 
    left_of_ip > 0 && return insert_point-Int32(1), Int32(left_of_ip)

    # Check right for a match, no need to check if it's outside the bounds of the SA <= length(sa)
    right_of_ip = insert_point <= length(sa) ? check_this_point(ca, sa, view(ref, ref_start:length(ref)), Int32(1), insert_point, Int32(0)) : 0
    right_of_ip > 0 && return insert_point, Int32(right_of_ip)

    # Neither actually matches our Ref, return 0,0 to move on to the next one
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
    # Given a point in the suffix array, compare the suffix to the Ref 
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
        # We can skip the LCP part when extending, note though we also have to 
        # check the previous match size so min(lcp valu, prev match size)
        start_check_from = Int32(min(lcp[i + lcp_dir], match_size))
        # Check the size of this match starting from +1 of the LCP value)
        match_size::Int32 = check_this_point(ca, sa, ref, ref_start, Int32(i), start_check_from )
        update_color!(color, ref_id, sa[i], Int32(match_size), ca)
        i += move_dir        
    end
end

function do_the_work(ref_id::Int32,  color::Color, ca::Vector{Int32}, sa::Vector{Int32}, s_idx::Int32, e_idx::Int32, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
    
    ref_start = s_idx

    while ref_start <= e_idx  

         # Do binary search to locate the insert point
         insert_point = locate_insert_point(sa, ca, view(ref, ref_start:length(ref)))
         max_match_index, max_match_size = decide_on_seed(insert_point, ca, sa, ref, ref_start)
 
         ref_suffix = view(ref, ref_start:ref_start+max_match_size-Int32(1))
 
         # If we have a match keep using the linked list to extend 
         if max_match_size > 0
             
             while ref_start <= e_idx
                 
                 # Check the match size at this point 
                 max_match_size = check_this_point(ca, sa, ref, ref_start, max_match_index, Int32(max_match_size-1)) # skip k-1
                 
                 # If we don't have any match we don't have to check the flanks
                 max_match_size == 0 && break 
 
                 ref_suffix = view(ref, ref_start:max_match_size)
 
                 # Else update color
                 update_color!(color, ref_id, sa[max_match_index], max_match_size, ca)
                               
                 # Check up and down in suffix array for other matches (LCP)
                 extend_from_point!(ca, sa, ref, lcp, max_match_index, false, ref_start, Int32(max_match_size), ref_id, color)
                 extend_from_point!(ca, sa, ref, lcp, max_match_index, true, ref_start, Int32(max_match_size), ref_id, color)
 
                 # Move to next location in suffix array for the second around (ISA)
                 max_match_index = inv_perm_sa[sa[max_match_index]+1]
    #             println("ref start increased")
                 ref_start += Int32(1)
                # next!(p)
             end 
         else 
             # No match at current point, move + 1 to do a binary search again
             #println("Will do binary search again: ", ref_start)
             ref_start += Int32(1)
             #next!(p)
         end
     end
end
 


function align_forward_and_reverse(ref_id::Int32, color::Color, ca::Vector{Int32}, sa::Vector{Int32}, start_idx::Int32, end_idx::Int32, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
    
    # First do the forward align 
    do_the_work(ref_id, color, ca, sa,  start_idx, end_idx, ref, inv_perm_sa, lcp)
    
    # Flip the nodes and reverse to do the reverse alignment 
    reverse_complement_ref!(ref)
    do_the_work(ref_id, color, ca, sa, start_idx, end_idx, ref, inv_perm_sa, lcp)

end


function run(gfa::String, seq_file::String, query_file::String, k_size::Int32, out_file::String; blacklist::String = "") 
    
    blacklist_ids = !isempty(blacklist) ? read_ids_from_file(blacklist) : OrderedSet{String}()

    # Read node size
    println("[STEP 1] Reading node sizes...")
    node_id_remapping =  Dict{Int64, Int32}()
    size_map = Dict{Int32, Int32}()
    try
        size_map = read_node_sizes(seq_file)
    catch 
        # We can try to fit the node ids in int32, else abort
        println("Probably overflow of node ids, trying to fit in int32s..")
        node_id_remapping, size_map = read_node_sizes_int32_fit_try(seq_file)
    end

    # Read node paths
    println("[STEP 2] Reading query paths...")
    queries, query_ids = processGFA(gfa, query_file, node_id_remapping)
    
    println("[STEP 3] Construction SA of queries")
    ca, sa = create_k_suffix_array(queries, Int32(0))

  #  println("CA: ", ca)

    println("[STEP 4] Construction ISA of queries")
    inv_sa_perm = inverse_perm_sa(sa)

    println("[STEP 5] Construction LCP from SA")
    lcp = build_lcp(sa, ca)

    info_vector = [ColorInfo(Base.Threads.SpinLock(), 0, -1, -1) for i in 1:length(ca)]
    color = Color(info_vector, size_map, k_size)

    println("[STEP 6] Starting graph file iteration..")
    
    # lets push suffixes to a channel
    nthreads = Threads.nthreads()
    
    for (ref_id, line) in enumerate(eachline(gfa))
        identifier, path = split(line, "\t")        
        if !(identifier in query_ids) && !(identifier in blacklist_ids)
            path_numbers = parse_numbers(path, node_id_remapping)
            nthreads = Threads.nthreads()
            chunk_size = div(length(path_numbers)  + nthreads - 1, nthreads)
            
            Threads.@threads for i in 1:nthreads
                start_idx = Int32((i - 1) * chunk_size + 1)
                end_idx = Int32(min(i * chunk_size, length(path_numbers)))
               # println("S: ", start_idx, " E:", end_idx, " total: ", length(path_numbers))
                align_forward_and_reverse(Int32(ref_id), color, ca, sa, start_idx, end_idx, path_numbers, inv_sa_perm, lcp)
            end

            

            println("$ref_id - $identifier done!")
            
        end 
    end


    println("[STEP 7] Writing results")
    writeResults(ca, color, query_ids, out_file, size_map)

end


# Profile.clear(); @profile run("/net/phage/linuxhome/mgx/people/rickb/Campy/complete_genomees/graph/complete_campy_graph.cf_seq", "/net/phage/linuxhome/mgx/people/rickb/Campy/complete_genomees/graph/complete_campy_graph.cf_seg", "/net/phage/linuxhome/mgx/people/rickb/Campy/complete_genomees/graph/coli_query2.txt"  , Int32(31), "t.txt")