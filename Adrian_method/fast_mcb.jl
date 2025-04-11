module FastMCB

using LinearAlgebra
using Graphs
using SimpleWeightedGraphs # For weighted graph and shortest paths
using SparseArrays        # For efficient GF(2) vectors
using DataStructures      # For PriorityQueue in Dijkstra if needed, although Graphs.jl handles it

# --- Helper Functions ---

"""
Calculate Euclidean distance between two points (vectors).
Assumes points have the same dimension.
"""
function euclidean_distance(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real})
    return norm(p1 - p2)
end

"""
Compute XOR sum (addition in GF(2)) of two sparse vectors, modifying v1 in place.
Assumes vectors are sparse, contain only 0s and 1s, and have the same length.
"""
function sparse_xor!(v1::SparseVector{Tv, Ti}, v2::SparseVector{Tv, Ti}) where {Tv<:Integer, Ti<:Integer}
    # Ensure v1's storage can accommodate potential new non-zeros from v2
    # This part is tricky with SparseVector, might be easier to create a new vector
    # Let's try creating a new vector for simplicity and correctness first.
    
    nz1 = v1.nzind
    nz2 = v2.nzind
    
    # Combine indices and values, handling XOR logic
    new_indices = Ti[]
    new_values = Tv[]
    
    idx1 = 1
    idx2 = 1
    len1 = length(nz1)
    len2 = length(nz2)
    
    while idx1 <= len1 || idx2 <= len2
        if idx1 <= len1 && (idx2 > len2 || nz1[idx1] < nz2[idx2])
            # Only in v1
            push!(new_indices, nz1[idx1])
            push!(new_values, v1.nzval[idx1]) # Assumes value is 1
            idx1 += 1
        elseif idx2 <= len2 && (idx1 > len1 || nz2[idx2] < nz1[idx1])
            # Only in v2
            push!(new_indices, nz2[idx2])
            push!(new_values, v2.nzval[idx2]) # Assumes value is 1
            idx2 += 1
        else # Indices match (nz1[idx1] == nz2[idx2])
            # Both have a 1 at this index, XOR is 0, so skip
            idx1 += 1
            idx2 += 1
        end
    end
    
    # Construct the new sparse vector
    if isempty(new_indices)
        return spzeros(Tv, v1.n)
    else
        # Need to ensure no explicit zeros are stored if Tv allows them
        filter!(x -> x != 0, new_values) # Should not be needed if input is GF(2) {0,1}
        return sparsevec(new_indices, new_values, v1.n)
    end
end


"""
Compute dot product of two sparse vectors over GF(2).
"""
function dot_gf2(v1::SparseVector{<:Integer, Ti}, v2::SparseVector{<:Integer, Ti}) where Ti
    common_indices = intersect(v1.nzind, v2.nzind)
    # In GF(2), only 1*1 contributes 1 to the sum.
    # So we just need the count of common indices where both are non-zero (assumed 1).
    return length(common_indices) % 2
end

"""
Reconstructs the path (as a list of edge indices) from Dijkstra parents array.
"""
function get_path_edges(parents::Vector{Ti}, target::Ti, source::Ti, graph_info) where Ti
    edge_indices = Ti[]
    curr = target
    edge_map = graph_info.edge_to_index

    while curr != source && parents[curr] != 0 # Check for 0 parent (unreachable)
        prev = parents[curr]
        # Find the edge index corresponding to (prev, curr) or (curr, prev)
        edge_tuple_fwd = (prev, curr)
        edge_tuple_rev = (curr, prev)

        edge_idx = 0
        if haskey(edge_map, edge_tuple_fwd)
            edge_idx = edge_map[edge_tuple_fwd]
        elseif haskey(edge_map, edge_tuple_rev)
             edge_idx = edge_map[edge_tuple_rev]
        else
             # This shouldn't happen if path exists in the original graph structure used for mapping
             @warn "Edge ($prev, $curr) not found in edge_to_index map during path reconstruction."
             return nothing # Indicate error or incomplete path
        end
        push!(edge_indices, edge_idx)
        curr = prev
    end

    if curr != source # Path didn't reach the source
         @warn "Path reconstruction failed to reach source node $source from target $target."
         return nothing
    end
    
    return reverse(edge_indices) # Reverse to get path from source to target
end

"""
Converts a list of edge indices into a sparse GF(2) vector.
"""
function cycle_to_vector(edge_indices::Vector{Ti}, num_total_edges::Ti) where Ti
    if isempty(edge_indices)
        return spzeros(Int, num_total_edges)
    end
    # Ensure indices are unique and within bounds
    unique_indices = unique(edge_indices)
     if any(idx -> idx < 1 || idx > num_total_edges, unique_indices)
         error("Edge index out of bounds [1, $num_total_edges].")
     end
    vals = ones(Int, length(unique_indices))
    return sparsevec(unique_indices, vals, num_total_edges)
end

"""
Converts a list of edge indices representing a cycle into a sequence of vertices.
Assumes the edges form a single simple cycle.
"""
function edge_indices_to_vertex_sequence(
    cycle_edge_indices::Vector{Ti}, 
    graph_info
    )::Union{Vector{Ti}, Nothing} where Ti

    if isempty(cycle_edge_indices)
        return Ti[]
    end

    # Build adjacency list for the cycle edges only
    adj = Dict{Ti, Vector{Ti}}()
    vertices_in_cycle = Set{Ti}()
    num_total_vertices = graph_info.num_vertices # Use the count from the graph potentially including isolated ones

    for edge_idx in cycle_edge_indices
        if edge_idx < 1 || edge_idx > graph_info.num_edges
             @warn "Invalid edge index $edge_idx found in cycle during vertex sequence conversion."
             return nothing
        end
        u, v = graph_info.edge_list[edge_idx]
        
        # Add to adjacency list (undirected)
        push!(get!(adj, u, Ti[]), v)
        push!(get!(adj, v, Ti[]), u)
        # Add vertices to set
        push!(vertices_in_cycle, u)
        push!(vertices_in_cycle, v)
    end

    # Basic validation: In a simple cycle, every vertex must have degree 2.
    for v in keys(adj)
        if length(adj[v]) != 2
            @warn "Cycle edges do not form a simple cycle (vertex $v has degree $(length(adj[v]))). Cannot reliably extract vertex sequence."
            # This might happen if the basis contains sums of cycles.
            # How to handle this? Return nothing, or try a more complex traversal?
            # For now, return nothing.
            return nothing 
        end
    end
    
    if isempty(vertices_in_cycle)
        return Ti[] # No vertices means empty cycle
    end

    # Start traversal from an arbitrary vertex
    start_node = first(vertices_in_cycle)
    path = Ti[start_node]
    prev_node = -1 # Sentinel value
    curr_node = start_node

    while length(path) <= length(vertices_in_cycle) # Prevent infinite loops
        neighbors = adj[curr_node]
        next_node = -1
        
        # Find the neighbor that isn't the previous node
        if neighbors[1] != prev_node
            next_node = neighbors[1]
        elseif length(neighbors) > 1 && neighbors[2] != prev_node # Should always have 2 neighbors here based on check above
            next_node = neighbors[2]
        else
             # This case implies an issue, possibly degree 1 or error in logic
             @warn "Traversal error: Could not find next node from $curr_node (prev=$prev_node, neighbors=$neighbors)."
             return nothing
        end

        if next_node == start_node
            # Completed the cycle
            if length(path) == length(vertices_in_cycle)
                 return path
             else
                 # Reached start too early - indicates issue (maybe not a single connected cycle?)
                 @warn "Traversal completed cycle prematurely. Expected $(length(vertices_in_cycle)) vertices, found $(length(path))."
                 return nothing
             end
        end

        push!(path, next_node)
        prev_node = curr_node
        curr_node = next_node
    end
    
    # If loop finished without returning, something went wrong (e.g., didn't close)
    @warn "Traversal failed to close the cycle within expected length."
    return nothing
end

# --- Core Algorithm Components ---

"""
Computes the minimum weight cycle C such that dot_gf2(C, S) == 1.
Implementation based on finding shortest paths between endpoints of edges e
where S[e]=1, using only edges f where S[f]=0.
"""
function compute_min_weight_cycle(
    S_vector::SparseVector{Int, Ti},
    graph_info
    )::Tuple{SparseVector{Int, Ti}, Float64} where Ti

    min_cycle_vec = spzeros(Int, graph_info.num_edges)
    min_weight = Inf

    # 1. Identify edges 'e' where S_vector[edge_index(e)] == 1.
    edges_in_S_indices = S_vector.nzind # Indices where S is 1

    if isempty(edges_in_S_indices)
        # If S is the zero vector, orthogonality cannot be 1.
        # Paper likely handles this case; maybe returns inf weight?
        @warn "S_vector is zero in compute_min_weight_cycle. Cannot satisfy orthogonality=1."
        return min_cycle_vec, min_weight
    end

    # 2. Create temporary graph G' with edges 'f' where S[f] == 0.
    g_prime = SimpleWeightedGraph(graph_info.num_vertices)
    for j in 1:graph_info.num_edges
        if S_vector[j] == 0 # Check if edge j is NOT in S
            u, v = graph_info.edge_list[j]
            weight = graph_info.weights[j]
            add_edge!(g_prime, u, v, weight)
        end
    end
    
    if ne(g_prime) == 0 && !isempty(edges_in_S_indices)
        # If G' is empty but S is not, no paths can be formed.
        @warn "Graph G' (edges not in S) is empty in compute_min_weight_cycle."
         return min_cycle_vec, min_weight
    end

    # 3. For each edge e=(u,v) in S, find shortest path P(u,v) in G'.
    for edge_idx in edges_in_S_indices
        u, v = graph_info.edge_list[edge_idx]
        edge_weight = graph_info.weights[edge_idx]

        # Check if endpoints u, v exist in g_prime
        if u > nv(g_prime) || v > nv(g_prime) || !has_vertex(g_prime, u) || !has_vertex(g_prime, v)
            # If endpoints are not in g_prime (isolated by removing S edges), no path exists
            continue
        end

        # Find shortest path in G'
        try
            dijkstra_result = dijkstra_shortest_paths(g_prime, u)
            
            path_dist = dijkstra_result.dists[v]

            if isfinite(path_dist)
                cycle_weight = path_dist + edge_weight
                if cycle_weight < min_weight
                    min_weight = cycle_weight
                    # Reconstruct path and create cycle vector
                    path_edge_indices = get_path_edges(dijkstra_result.parents, v, u, graph_info)
                     if path_edge_indices !== nothing
                        cycle_edge_indices = vcat(path_edge_indices, [edge_idx])
                        # Convert to sparse vector, handling potential duplicate edges in cycle via XOR logic implicitly
                        current_cycle_vec = spzeros(Int, graph_info.num_edges)
                        for idx in cycle_edge_indices
                            current_cycle_vec[idx] = 1 - current_cycle_vec[idx] # XOR toggle
                        end
                        # Remove explicit zeros if any crept in (shouldn't with 0/1 toggling)
                        dropzeros!(current_cycle_vec)
                        
                        min_cycle_vec = current_cycle_vec
                    else
                         @warn "Path reconstruction failed for edge $edge_idx between $u and $v, though distance was finite."
                         # Reset min_weight if the best path couldn't be reconstructed? Or just skip this candidate?
                         # Let's skip this candidate for now.
                         continue
                    end
                end
            end
        catch e
             if isa(e, KeyError) || isa(e, BoundsError)
                 # Handle cases where a vertex might not be in the graph or paths are impossible
                 # These might occur if the graph becomes disconnected significantly
                 # Or if vertex indexing is inconsistent
                 @warn "Error during Dijkstra for edge $edge_idx ($u, $v): $e. Skipping this candidate."
                 continue
             else
                 rethrow(e)
             end
         end
    end
    
     if !isfinite(min_weight)
        @warn "compute_min_weight_cycle did not find any valid cycle."
    end

    return min_cycle_vec, min_weight
end


"""
Updates the second list of S-vectors (S2_original_list) to be orthogonal
to the newly found cycles (cycles_from_call1), using the updated S1 vectors
(S1_out) perhaps as part of the basis transformation (details depend on paper).
This version implements the Gaussian elimination style update.
"""
function update(
    S1_out::Vector{SparseVector{Int, Ti}}, # Updated S vectors from first half (role depends on paper)
    S2_original_list::Vector{SparseVector{Int, Ti}},
    cycles_from_call1::Vector{SparseVector{Int, Ti}},
    graph_info # Not strictly needed for this GF(2) update logic
    )::Vector{SparseVector{Int, Ti}} where Ti

    T2_list = deepcopy(S2_original_list) # Start with copies

    # For each S_j in the second half, make it orthogonal to all C_p found so far.
    for p in 1:length(cycles_from_call1)
        Cp = cycles_from_call1[p]
        if nnz(Cp) == 0 continue end # Skip empty cycles

        for j in 1:length(T2_list)
            Tj = T2_list[j]
            if dot_gf2(Cp, Tj) == 1
                # If not orthogonal, add Cp to Tj (XOR)
                T2_list[j] = sparse_xor!(Tj, Cp) # Update Tj in the list
            end
        end
    end

    return T2_list
end

# --- Recursive Core Function ---

"""
Recursive procedure to extend the cycle basis. Modifies C_basis in place.
"""
function extend_cb!(
    C_basis::Vector{SparseVector{Int, Ti}},
    S_list::Vector{SparseVector{Int, Ti}},
    k::Int,
    graph_info
    )::Vector{SparseVector{Int, Ti}} where Ti

    if k == 0
        return SparseVector{Int, Ti}[] # Return empty list for S
    end

    if k == 1
        if isempty(S_list)
             error("S_list is empty when k=1")
        end
        S_current = S_list[1]

        min_cycle_vector, weight = compute_min_weight_cycle(S_current, graph_info)

        if nnz(min_cycle_vector) > 0 # Check if a non-empty cycle was found
            push!(C_basis, min_cycle_vector)
            # Paper needs to clarify what happens to S_current. Assume it's 'consumed' / transformed.
            # Returning the original S_current for now as a placeholder for the "updated" S from this branch.
             return [S_current] 
        else
            @warn "Failed to find min weight cycle for k=1 base case. S = $(findnz(S_current)[1])"
            # How to handle failure? Paper might specify. Return original S?
            return [S_current]
        end

    else # k > 1, use recursion
        k1 = k ÷ 2 # Floor division
        k2 = k - k1

        # Split S_list
        S1_in = S_list[1:k1]
        S2_in = S_list[k1+1:end]

        # Track cycles added by the first call
        num_cycles_before_call1 = length(C_basis)

        # --- Recursive call 1 ---
        S1_out = extend_cb!(C_basis, S1_in, k1, graph_info)
        
        num_cycles_after_call1 = length(C_basis)
        cycles_from_call1 = C_basis[num_cycles_before_call1 + 1 : num_cycles_after_call1]

        # --- Update step ---
        T2_list = update(S1_out, S2_in, cycles_from_call1, graph_info)

        # --- Recursive call 2 ---
        S2_out = extend_cb!(C_basis, T2_list, k2, graph_info)

        # Combine the final S lists from both branches
        # This assumes the returned S vectors are the ones corresponding to the found basis elements.
        return vcat(S1_out, S2_out)
    end
end

# --- Main Function ---

"""
    fast_mcb(num_vertices, edges, centroids)

Computes a Minimum Weight Cycle Basis using the FAST-MCB algorithm structure.
Edge weights are Euclidean distances between centroids.

Args:
    num_vertices (Int): Number of vertices (centroids).
    edges (Vector{Tuple{Int, Int}}): List of directed edges (u, v). Edges are treated
                                     as undirected for basis computation but directionality
                                     might be relevant for specific paper details omitted here.
    centroids (Vector{Vector{Float64}}): List of centroid coordinates.

Returns:
    Vector{Vector{Int}}: A list of cycles, where each cycle is represented as a
                         sequence of vertex indices. Returns empty list on error or if no cycles found.
"""
function fast_mcb(num_vertices::Int, edges::Vector{Tuple{Int, Int}}, centroids::Vector{<:AbstractVector{<:Real}})
    num_edges = length(edges)
    if num_edges == 0
        return Vector{Vector{Int}}[] # Return empty basis if no edges
    end
     num_centroids = length(centroids)
     # Simple check: Ensure centroids cover at least the max vertex index used in edges.
     max_vertex_idx_in_edges = 0
     if !isempty(edges)
         max_vertex_idx_in_edges = maximum(max(u, v) for (u, v) in edges)
     end
 
     # Validate inputs
     if num_centroids < max_vertex_idx_in_edges
         error("Number of centroids ($num_centroids) is less than the maximum vertex index used in edges ($max_vertex_idx_in_edges).")
     end
     # Use the maximum index found in edges or num_vertices if it's larger (covers isolated vertices)
     effective_num_vertices = max(num_vertices, max_vertex_idx_in_edges)
     if effective_num_vertices > num_centroids
         @warn "Effective number of vertices ($effective_num_vertices) is greater than provided centroids ($num_centroids). Check vertex numbering."
         # Decide how to handle this: error or proceed with caution? Error is safer.
         error("Mismatch between vertex count and centroid data.")
     end
     if num_vertices < max_vertex_idx_in_edges
        @warn "Provided num_vertices ($num_vertices) is less than max index in edges ($max_vertex_idx_in_edges). Using $max_vertex_idx_in_edges."
        num_vertices = max_vertex_idx_in_edges
     end


    println("Starting FAST-MCB Implementation...")

    # 1. Create edge map and calculate weights
    edge_list = Vector{Tuple{Int, Int}}(undef, num_edges)
    edge_weights = Vector{Float64}(undef, num_edges)
    edge_to_index = Dict{Tuple{Int, Int}, Int}() 
    undirected_edge_set = Set{Tuple{Int,Int}}() # To build undirected graph correctly

    for (j, edge) in enumerate(edges)
        u, v = edge
        # Validate indices against centroid list length (already partially done above)
         if u < 1 || u > num_centroids || v < 1 || v > num_centroids
             error("Edge vertex index ($u or $v) out of bounds [1, $num_centroids] based on centroids list.")
         end
         
        edge_list[j] = edge
        # Store original directed edge for mapping, but check undirected for graph building
        edge_to_index[edge] = j 
        
        # Use canonical representation (min_idx, max_idx) for undirected check
        u_std, v_std = minmax(u, v)
        push!(undirected_edge_set, (u_std, v_std))

        # Calculate weight (Euclidean distance)
        weight = euclidean_distance(centroids[u], centroids[v])
        edge_weights[j] = weight
    end
    
    # 2. Build weighted undirected graph for shortest paths
    g = SimpleWeightedGraph(num_vertices) 
    temp_edge_weights = Dict{Tuple{Int,Int}, Float64}() # Store weights for undirected graph
    
    for j in 1:num_edges
        u, v = edge_list[j]
        weight = edge_weights[j]
        u_std, v_std = minmax(u, v)
        edge_key = (u_std, v_std)
        
        # Use the minimum weight if the edge appears in both directions (unlikely with the 'edges' function logic)
        current_weight = get(temp_edge_weights, edge_key, Inf)
        temp_edge_weights[edge_key] = min(current_weight, weight)
    end
    
    for (edge, weight) in temp_edge_weights
        add_edge!(g, edge[1], edge[2], weight)
    end

    # Container for graph information
    graph_info = (
        graph=g,
        weights=edge_weights, # Weights of original directed edges
        edge_list=edge_list, # Original directed edges
        edge_to_index=edge_to_index, # Map directed edge tuple to index
        num_vertices=nv(g), 
        num_edges=num_edges
    )

    # 3. Initialize S_j vectors (GF(2) sparse vectors)
    S_initial = Vector{SparseVector{Int, Int}}(undef, num_edges)
    for j in 1:num_edges
        S_initial[j] = sparsevec([j], [1], num_edges) # Vector j has 1 at index j
    end

    # 4. Initialize cycle basis C
    C_basis_vectors = Vector{SparseVector{Int, Int}}() 

    # 5. Call the recursive function
    println("Calling extend_cb!...")
    try
        # The second argument S_initial will be modified by the process if using in-place updates in `update`
        # Let's pass a copy if we want to preserve the original S_initial
        final_S_list = extend_cb!(C_basis_vectors, deepcopy(S_initial), num_edges, graph_info)
        println("extend_cb! finished.")
    catch e
        println("Error during extend_cb!: $e")
        showerror(stdout, e, catch_backtrace())
        println()
        return Vector{Vector{Int}}[] # Return empty on error
    end


    # 6. Convert basis vectors (GF(2)) to lists of vertex sequences
    println("Converting basis vectors to vertex sequence lists...")
    C_basis_vertex_sequences = Vector{Vector{Int}}()
    failed_conversions = 0
    for cycle_vector in C_basis_vectors
        # Find edge indices for this cycle
        indices, _ = findnz(cycle_vector)
        if isempty(indices) continue end # Skip empty cycles

        # Convert edge indices to vertex sequence
        vertex_sequence = edge_indices_to_vertex_sequence(indices, graph_info)

        if vertex_sequence !== nothing && !isempty(vertex_sequence)
            push!(C_basis_vertex_sequences, vertex_sequence)
        else
            failed_conversions += 1
            @warn "Failed to convert one cycle vector to a valid vertex sequence."
            # Option: Store edge indices instead? Or just skip? Skipping for now.
            # push!(C_basis_vertex_sequences, indices) # Fallback to edge indices?
        end
    end
     if failed_conversions > 0
         println("Warning: $failed_conversions cycle(s) could not be converted to simple vertex sequences.")
     end

    println("FAST-MCB finished. Found $(length(C_basis_vertex_sequences)) basis cycles as vertex sequences.")
    
    # Optional: Calculate total weight of the basis found (using edge indices before conversion)
     total_weight = 0.0
     for cycle_vector in C_basis_vectors
          indices, _ = findnz(cycle_vector)
          if !isempty(indices)
            cycle_weight = sum(graph_info.weights[idx] for idx in indices)
            total_weight += cycle_weight
          end
     end
     println("Total weight of the computed basis (sum of cycle edge weights): $total_weight")


    return C_basis_vertex_sequences # Return list of vertex sequences
end

end # module FastMCB
