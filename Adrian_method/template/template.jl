module Template

using LinearAlgebra
using Statistics
using Graphs
using ParallelKMeans
using SparseArrays
using StatsBase
include("networkx.jl")

"""
  best_complex(
    trajectory,
    k;
    rounds=10,
    maxiter=nothing,
    connectivity_threshold=30,
    max_force_boundary=0,
    print_frequencies=false,
    allow_empty_joining_locus=false
  )

Run multiple rounds of clustering and complex construction, returning the result
with the smallest joining locus (i.e., the best fitness).

Each round, the trajectory is clustered into `k` clusters using ParallelKMeans,
and the resulting complex is evaluated by the length of its joining locus (edges
shared by 3 or more faces). The fitness function is the negative of the joining
locus length (so maximizing fitness minimizes the joining locus). The best
result across all rounds is returned.

Args:
- trajectory: Matrix where rows are time points and columns are dimensions.
- k: Number of clusters.

Keyword args:
- rounds::Int=10: Number of clustering rounds to try (each with a different
  k-means initialization).
- maxiter::Union{Int, Nothing}=nothing: Maximum iterations for k-means (passed
  to ParallelKMeans).
- connectivity_threshold::Int=30: Minimum number of transitions required for an
  edge between clusters.
- max_force_boundary::Int=0: Minimum cycle length to force as a 2-cell boundary
  (passed to faces()).
- print_frequencies::Bool=false: If true, print the frequencies of joining locus
  lengths across rounds.
- allow_empty_joining_locus::Bool=false: If false (default), results with joining locus length 0 are skipped and not considered as best. If true, such results are allowed.

Returns:
- A named tuple containing the best result, with fields:
    - clusters: Vector of (centroid, point_idxs) tuples.
    - assignments: Raw assignments vector from k-means.
    - centroids: Vector of centroid coordinates.
    - connectivity_matrix: Sparse matrix of cluster transitions.
    - edges: Vector of directed edges.
    - mcb: Minimal cycle basis (vector of cycles).
    - faces: Vector of 2-cell faces (unique cycles).
    - joining_locus: Vector of edges in the joining locus.
    - joining_locus_length: Number of edges in the joining locus.
    - round: The round in which the best result was found.
"""
function best_complex(
  trajectory,
  k;
  rounds::Int=10,
  maxiter::Union{Int, Nothing}=nothing,
  connectivity_threshold::Int=30,
  max_force_boundary::Int=0,
  print_frequencies::Bool=false,
  allow_empty_joining_locus::Bool=false
)::NamedTuple
  best_fitness = -Inf
  best_result = nothing
  joining_locus_lengths = Int[]
  best_round = 0

  for round in 1:rounds
    println("[best_complex] Round $round of $rounds...")
    # Cluster the trajectory.
    clusters, assignments = cluster(trajectory, k; maxiter=maxiter)
    centroids = [c[1] for c in clusters]

    # Compute connectivity and edges.
    connectivity_matrix = cluster_connectivity(assignments, k)
    edges = Template.edges(connectivity_matrix, connectivity_threshold)

    # Skip if not enough edges for a complex.
    if length(edges) < k
      push!(joining_locus_lengths, 0)
      continue
    end

    # Build NetworkX graph and compute minimal cycle basis.
    nx_graph = NetworkX.edges_to_networkx(edges, centroids)
    mcb = NetworkX.mcb(nx_graph)

    # Compute faces (2-cells).
    two_cells = faces(mcb, edges; max_force_boundary=max_force_boundary)
    # Get unique faces.
    distinct_faces = unique(two_cells)
    # Compute joining locus.
    joining_locus_edges = joining_locus(edges, distinct_faces)
    joining_locus_length = length(joining_locus_edges)
    push!(joining_locus_lengths, joining_locus_length)

    # Optionally skip empty joining locus results.
    if !allow_empty_joining_locus && joining_locus_length == 0
      println("[best_complex] WARNING: Joining locus is empty (bad templex) in round $round. Skipping.")
      continue
    end

    # Fitness is negative joining locus length (maximize = minimize locus).
    fitness = -joining_locus_length
    if fitness > best_fitness
      best_fitness = fitness
      best_result = (
        clusters=clusters,
        assignments=assignments,
        centroids=centroids,
        connectivity_matrix=connectivity_matrix,
        edges=edges,
        mcb=mcb,
        faces=distinct_faces,
        joining_locus=joining_locus_edges,
        joining_locus_length=joining_locus_length,
        round=round
      )
      best_round = round
    end
  end

  # Optionally print frequencies of joining locus lengths.
  if print_frequencies
    freq = StatsBase.countmap(joining_locus_lengths)
    println("Joining locus length frequencies:")
    for (len, count) in sort(collect(freq), by=x->x[1])
      println("  $len: $count")
    end
  end

  return best_result
end

"""
  cluster(trajectory, k; [maxiter])

Cluster the trajectory (except for the first and last point) into k clusters in
1-jet space using the ParallelKMeans algorithm.

Args:
- trajectory: Matrix where rows are time points and columns are dimensions.
- k: Number of clusters.
- maxiter::Union{Int, Nothing} (keyword, default nothing): Maximum iterations for k-means.

Returns a tuple:
- A vector of tuples (centroid, point_idxs). Centroids are 3D spatial coordinates.
- The raw assignments vector from k-means.
"""
function cluster(trajectory, k; maxiter::Union{Int, Nothing}=nothing)
  dim = size(trajectory, 2)
  if dim < 3
      error("Trajectory must have at least 3 spatial dimensions.")
  end
  points = trajectory[2:end-1, :]
  # Use next - prev for velocity.
  velocities = trajectory[3:end, :] - trajectory[1:end-2, :]
  # Create jets by concatenating points and velocities.
  # Data format for kmeans/fastkmeans/ParallelKMeans: features (dims) as rows, observations as columns
  jets = hcat(points, velocities)' # jets is (2*dim) x num_points

  # Cluster the jets using ParallelKMeans.
  local result # Ensure scope
  if maxiter === nothing
    result = ParallelKMeans.kmeans(jets, k) # Use default max_iters
  else
    result = ParallelKMeans.kmeans(jets, k, max_iters=maxiter) # Pass maxiter as max_iters
  end

  # Extract assignments and centroids from the result struct
  assignments = result.assignments
  centroids = result.centers # centroids is (2*dim) x k

  # Group point indices by cluster. Indices are relative to the 'points' array (1 to num_points-2).
  point_idxs = [findall(==(i), assignments) for i in 1:k]

  # Format cluster info, extracting only the spatial dimensions (first `dim`) for the centroid.
  # Ensure centroids is accessed correctly (it's (2*dim) x k)
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
- max_force_boundary::Int (keyword, default 0): If the number of vertices `n`
  in `cycle_vertices` is less than this value, the function will
  immediately return `true` (treating it as not consistently orientable/not a 2-cell).

Returns:
- A boolean indicating whether the cycle vertices fail to form a consistently
  oriented cycle.
"""
function cycle_is_2cell(
  cycle_vertices::Vector{Int},
  edges::Vector{Tuple{Int, Int}};
  max_force_boundary::Int=0
)
  n = length(cycle_vertices)

  # Force small cycles to be boundaries of 2-cells.
  if n <= max_force_boundary
    return true
  end

  # A cycle needs at least 3 vertices for orientation checks.
  if n < 3
    return false # Treat triangles and smaller as potentially orientable by default
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
- max_force_boundary::Int (keyword, default 0): If the number of vertices `n`
  in `cycle_basis` is less than this value, the function will immediately return
  the cycle (treating it as not consistently orientable, i.e., a boundary of a
  2-cell).

Returns:
- A vector of vectors, where each inner vector represents a 2-cell boundary as a
  sequence of vertices.
"""
function faces(cycle_basis, edges; max_force_boundary::Int=0)
  # Filter the cycle basis, keeping only the cycles that are 2-cells.
  two_cells = [
    cycle for cycle in cycle_basis if cycle_is_2cell(
      cycle,
      edges;
      max_force_boundary=max_force_boundary
    )
  ]
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

"""
  mcb(edges, k)

Constructs a minimal cycle basis for a graph represented by edges.

Args:
- edges: A vector of 2-tuples representing directed edges.
The function treats these as undirected for cycle basis calculation.
- k: The number of vertices in the graph.

Returns:
- A vector of vectors, where each inner vector represents a cycle as a sequence 
  of vertices.
"""
function mcb(edges, k)
  # Create an undirected graph with k vertices.
  # SimpleGraph should be available via `using Graphs`.
  g = SimpleGraph(k)

  # Add edges one by one to the undirected graph.
  # add_edge! should be available via `using Graphs`.
  for edge in edges
    u, v = edge
    add_edge!(g, u, v) # Graphs.jl handles potential duplicate edges automatically.
  end

  # Compute the cycle basis of the undirected graph.
  # cycle_basis should be available via `using Graphs`.
  basis = cycle_basis(g)

  # Return the basis.
  return basis
end

end # module