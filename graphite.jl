const MASK = Int32(1<<30) 

function flipnode(n::Int32)
    n < 2^30 - 1 || error("Encounterd node ID that is too large (>= 2^30 - 1)")
    n ⊻ MASK
end
isflipped(n::Int32) = ifelse(n & MASK != 0, true, false)
noflip(n::Int32) = n & (MASK - 1)
convert_node(n::Int32) = n < 0 ? flipnode(abs(n)) : noflip(n) 
get_original(n::Int32) = isflipped(n) ? flipnode(n) : flipnode(n)

struct Origin 
    id::Int32 
    pos::Int32 
end

struct Color 
    len::Vector{UInt32}
    origin::Vector{Origin}
    size_map::Dict{Int32, UInt32}
    k_size::Int32
end

function update_color!(color::Color, ref_id::Int32, match_start::Int32, match_size::Int32, ca::Vector{Int32})
    
    match_end = match_start+match_size-1
    
    # Get info from pervious matches
    at_start = color.origin[match_start]
    at_end   = color.origin[match_end]

    if at_start.id > 0 && at_start.id == at_end.id && at_start.pos == at_end.pos
        return
    # Don't have to bother about single node matches as they can't be longer anyway
    elseif match_start ==  match_end && color.len[match_start] > 0
        return
    else
        match_size_nt = sum(get.(Ref(color.size_map), view(ca, match_start:match_end), 0))
        # We have to consider to overlap between k-mers as well 
        match_size_nt = match_size_nt - ((match_size-1) * (color.k_size-1)) 
        for i in match_start:match_end
            if color.len[i] < match_size_nt 
                color.len[i]  = match_size_nt
                color.origin[i] = Origin(ref_id, match_start)
            end
        end
    end
end

function reverse_complement_ref!(ref::Vector{Int32})
    reverse!(ref)
    @inbounds for i in eachindex(ref)
        ref[i] = flipnode(ref[i])
    end
end

function convert_nodes!(in_vector::Vector{Int32})
    for i in eachindex(in_vector)
        in_vector[i] = convert_node(in_vector[i])
    end
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
    match_size = matches_till(ref, ref_start, ca, ca_start) + skip
    return match_size
end

function extend_from_point!(ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, lcp::Vector{Int32}, point::Int32, forward::Bool, ref_start::Int32, match_size::Int32, ref_id::Int32, color::Color)
    move_dir = forward ? 1 : -1
    lcp_dir  = forward ? 0 :  1
   # println("Extending this point")
    i = point += move_dir
    while i > 1 && i <= length(sa) && lcp[i + lcp_dir] > 0
        # We can skip the LCP part when extending, note though we also have to 
        # check the previous match size so min(lcp valu, prev match size)
        start_check_from = Int32(min(lcp[i + lcp_dir], match_size))
        # Check the size of this match starting from +1 of the LCP value)
        match_size = check_this_point(ca, sa, ref, ref_start, Int32(i), start_check_from )
        update_color!(color, ref_id, sa[i], Int32(match_size), ca)
        i += move_dir        
    end
end

function align(ref_id::Int32, color::Color, ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
   # println(ref_id, " at: ", ref[307])
    max_match_index = Int32(0)
    max_match_size = Int32(0)
    ref_start = Int32(1)
    while ref_start <= length(ref)

        # Do binary search to locate the insert point
        insert_point = locate_insert_point(sa, ca, view(ref, ref_start:length(ref)))
        max_match_index, max_match_size = decide_on_seed(insert_point, ca, sa, ref, ref_start)

        # If we have a match keep using the linked list to extend 
        if max_match_size > 0 
            while ref_start <= length(ref)
                       # Check the match size at this point 
                max_match_size = check_this_point(ca, sa, ref, ref_start, max_match_index, Int32(max_match_size-1)) # skip k-1
                
                # If we don't have any match we don't have to check the flanks
                max_match_size == 0 && break 
                update_color!(color, ref_id, sa[max_match_index], Int32(max_match_size), ca)
                              
                # Check up and down in suffix array for other matches
                extend_from_point!(ca, sa, ref, lcp, max_match_index, false, ref_start, Int32(max_match_size), ref_id, color)
                extend_from_point!(ca, sa, ref, lcp, max_match_index, true, ref_start, Int32(max_match_size), ref_id, color)

                # # Move to next location in suffix array for the second around
                max_match_index = inv_perm_sa[sa[max_match_index]+1]
                ref_start += Int32(1)
            end 
        else 
            # No match at current point, move +1 to do a binary search again
            ref_start += Int32(1)
        end
    end
end

function align_forward_and_reverse(ref_id::Int32, color::Color, ca::Vector{Int32}, sa::Vector{Int32}, ref::Vector{Int32}, inv_perm_sa::Vector{Int32}, lcp::Vector{Int32})
    # First do the forward align 
    align(ref_id, color, ca, sa, ref, inv_perm_sa,lcp)
    # Flip the nodes and reverse to do the reverse alignment 
    reverse_complement_ref!(ref)
    align(ref_id, color, ca, sa, ref, inv_perm_sa,lcp)
end

function run(gfa::String, seq_file::String, query_file::String, k_size::Int32, out_file::String; blacklist::String = "") 
    
    blacklist_ids = !isempty(blacklist) ? read_ids_from_file(blacklist) : OrderedSet{String}()

    queries, query_ids = processGFA(gfa, query_file)
       
    ca, sa = create_k_suffix_array(queries, Int32(0))
    inv_sa_perm = inverse_perm_sa(sa)
    lcp = build_lcp(sa, ca)

  
    size_map = read_node_sizes(seq_file)
    len = zeros(Int32, length(ca))
    ori =  [Origin(-1,-1) for i in 1:length(ca)] # Vector{Origin}(undef, length(ca)) 
    color = Color(len, ori, size_map, k_size)

    for (ref_id, line) in enumerate(eachline(gfa))
       identifier, path = split(line, "\t")
       if !(identifier in query_ids) && !(identifier in blacklist_ids)
           path_numbers = parse_numbers(path)
           align_forward_and_reverse(Int32(ref_id), color, ca, sa, path_numbers, inv_sa_perm, lcp)
       end
    end

    writeResults(ca, color, query_ids, out_file, size_map)

end

# run("test", "test", "test", Int32(31), "test")
