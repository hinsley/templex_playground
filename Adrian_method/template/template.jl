module Template

using Clustering
using LinearAlgebra
using Statistics
using Graphs
using SparseArrays

"""
  cluster(trajectory, k)

Cluster the trajectory (except for the first and last point) into k clusters in
1-jet space.

Returns a tuple:
- A vector of tuples (centroid, point_idxs).
- The raw assignments vector from k-means.
"""
function cluster(trajectory, k)
  dim = size(trajectory, 2)
  points = trajectory[2:end-1, :]
  # Use next - prev for velocity.
  velocities = trajectory[3:end, :] - trajectory[1:end-2, :]
  # Create jets by concatenating points and velocities.
  jets = hcat(points, velocities)'

  # Cluster the jets.
  result = kmeans(jets, k)
  
  # Extract centroids and assignments.
  centroids = result.centers
  assignments = result.assignments
  
  # Group point indices by cluster. Indices are relative to the 'points' array.
  point_idxs = [findall(==(i), assignments) for i in 1:k]
  
  # Format cluster info.
  cluster_info = [(centroids[1:dim, i], point_idxs[i]) for i in 1:k]

  # Return cluster info and raw assignments.
  return cluster_info, assignments
end

"""
  cluster_connectivity(assignments, k)

Construct a sparse matrix representing the transitions between clusters based on
the time-ordered assignments.

Args:
- assignments: A vector of integers representing the cluster assignment for each 
               point (typically trajectory[2:end-1]).
- k: The total number of clusters.

Returns:
- A SparseMatrixCSC where `M[i, j]` contains the number of transitions observed
  from cluster `i` to cluster `j`.
"""
function cluster_connectivity(assignments, k)
  n_assignments = length(assignments)
  if n_assignments < 2
      # Need at least two points to have a transition.
      return spzeros(Int, k, k) 
  end

  # I, J store the row and column indices for transitions.
  # V stores the count for each transition (initially 1).
  rows = Int[]
  cols = Int[]
  vals = Int[] 

  for i in 1:(n_assignments - 1)
    from_cluster = assignments[i]
    to_cluster = assignments[i+1]
    
    # Only record transitions between different clusters.
    if from_cluster != to_cluster
      # Record the transition.
      push!(rows, from_cluster)
      push!(cols, to_cluster)
      push!(vals, 1) # Each transition counts as 1.
    end
  end

  # Create the sparse matrix, summing counts for duplicate transitions.
  # The size is k x k.
  connectivity_matrix = sparse(rows, cols, vals, k, k)

  return connectivity_matrix
end

"""
  cycle_is_2cell(cycle_vertices, edges)

Check if a cycle, given as a list of vertices, cannot be traversed consistently
using the provided directed edges.

A cycle `[v1, v2, ..., vn]` is considered a 2-cell depending on whether either
all directed edges `(v1, v2), (v2, v3), ..., (vn, v1)` exist, or all reversed
directed edges `(v2, v1), (v3, v2), ..., (v1, vn)` exist in the `edges` list.
This function returns true if NEITHER of these conditions holds.

Args:
- cycle_vertices: A vector of vertex indices representing the cycle (output from
  `FastMCB`).
- edges: A vector of tuples `(u, v)` representing the directed edges in the
  graph.

Returns:
- A boolean indicating whether the cycle vertices fail to form a consistently
  oriented cycle.
"""
function cycle_is_2cell(
  cycle_vertices::Vector{Int},
  edges::Vector{Tuple{Int, Int}}
)
  n = length(cycle_vertices)
  # A cycle needs at least 3 vertices.
  if n < 3
    # Should this be true or false? If it cannot be a 2-cell, maybe true?
    # Let's assume for now it's not considered a failure in the same way, return false.
    return false 
  end

  # Create a set for efficient edge lookup.
  directed_edges = Set(edges)

  # Check forward consistency: v_i -> v_{i+1}
  forward_consistent = true
  for i in 1:n
    u = cycle_vertices[i]
    # Use mod1 for 1-based indexing wrap-around.
    v = cycle_vertices[mod1(i + 1, n)] 
    if (u, v) ∉ directed_edges
      forward_consistent = false
      break
    end
  end

  if forward_consistent
    return false # It IS consistently oriented (forward), so return false.
  end

  # Check backward consistency: v_{i+1} -> v_i
  backward_consistent = true
  for i in 1:n
    u = cycle_vertices[i]
    v = cycle_vertices[mod1(i + 1, n)]
    if (v, u) ∉ directed_edges
      backward_consistent = false
      break
    end
  end

  # Return true ONLY if NEITHER forward nor backward traversal was consistent.
  return !backward_consistent 
end

"""
  edges(connectivity_matrix, threshold=30)

Extract directed edges from a connectivity matrix based on transition counts.

An edge `(i, j)` is created if the number of transitions from cluster `i` to 
cluster `j` (`connectivity_matrix[i, j]`) meets two conditions:
1. It is greater than or equal to the `threshold`.
2. It is greater than or equal to the number of transitions in the reverse 
   direction (`connectivity_matrix[j, i]`).

Args:
- connectivity_matrix: A SparseMatrixCSC representing directed cluster
  transitions.
- threshold: The minimum number of transitions required for an edge.

Returns:
- A vector of 2-tuples, where each tuple `(from_cluster, to_cluster)`
  represents a directed edge satisfying the conditions.
"""
function edges(connectivity_matrix, threshold=30)
  # Ensure the input is a SparseMatrixCSC for efficient indexing.
  if !(connectivity_matrix isa SparseMatrixCSC)
    connectivity_matrix = sparse(connectivity_matrix)
  end

  # findnz returns row indices, column indices, and values of non-zero elements.
  rows, cols, vals = findnz(connectivity_matrix)
  
  filtered_edges = Tuple{Int, Int}[]

  for k in eachindex(rows)
    i = rows[k]
    j = cols[k]
    val_ij = vals[k] # This is connectivity_matrix[i, j]

    # Check if the transition count meets the threshold.
    if val_ij >= threshold
      # Get the transition count in the reverse direction.
      # Accessing a sparse matrix with [j, i] will return 0 if the entry doesn't
      # exist.
      val_ji = connectivity_matrix[j, i]
      
      # Check if the forward count is greater than or equal to the reverse
      # count.
      if val_ij >= val_ji
        push!(filtered_edges, (i, j))
      end
    end
  end

  return filtered_edges
end

"""
  faces(cycle_basis, edges)

Given a minimal cycle basis, return the list of 2-cells (faces) in the graph.

Args:
- cycle_basis: A vector of vectors, where each inner vector represents a cycle
  as a sequence of vertices.
- edges: A vector of 2-tuples representing directed edges.

Returns:
- A vector of vectors, where each inner vector represents a 2-cell boundary as a
  sequence of vertices.
"""
function faces(cycle_basis, edges)
  # Filter the cycle basis, keeping only the cycles that are 2-cells.
  two_cells = [cycle for cycle in cycle_basis if cycle_is_2cell(cycle, edges)]
  return two_cells
end

"""
  joining_locus(edges, faces)

Given a vector of edges, return the joining locus. The joining locus is the set
of edges to which 3 or more 2-cells are attached.

Args:
- edges: A vector of 2-tuples representing directed edges.
- faces: A vector of vectors, where each inner vector represents a 2-cell boundary
  as a sequence of vertices.

Returns:
- A vector of 2-tuples representing the directed edges in the joining locus.
"""
function joining_locus(edges::Vector{Tuple{Int, Int}}, faces::Vector{Vector{Int}})
  # Dictionary to store the count of faces attached to each undirected edge.
  # Key: Tuple{Int, Int} (canonical representation: min_vertex, max_vertex)
  # Value: Int (count)
  edge_face_count = Dict{Tuple{Int, Int}, Int}()

  # Count face attachments for each edge.
  for face_vertices in faces
    n = length(face_vertices)
    if n < 3 continue end # Skip invalid faces

    for i in 1:n
      u = face_vertices[i]
      v = face_vertices[mod1(i + 1, n)] # Next vertex, wrapping around

      # Use canonical representation for the undirected edge.
      edge_key = minmax(u, v) 

      # Increment count for this edge.
      edge_face_count[edge_key] = get(edge_face_count, edge_key, 0) + 1
    end
  end

  # Identify undirected edges belonging to the joining locus (count >= 3).
  joining_locus_keys = Set{Tuple{Int, Int}}()
  for (edge_key, count) in edge_face_count
    if count >= 3
      push!(joining_locus_keys, edge_key)
    end
  end

  # Filter the original directed edges based on the joining locus keys.
  result_edges = Vector{Tuple{Int, Int}}()
  for edge in edges
    u, v = edge
    edge_key = minmax(u, v)
    if edge_key in joining_locus_keys
      push!(result_edges, edge)
    end
  end

  return result_edges
end

end # module