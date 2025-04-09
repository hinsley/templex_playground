using Pkg
Pkg.activate(".")
Pkg.instantiate()

using GLMakie, OrdinaryDiffEq, Statistics, Clustering, DataStructures
using LinearAlgebra, KernelDensity, NearestNeighbors

# Define the Lorenz system
function lorenz!(du, u, p, t)
    σ, ρ, β = p
    x, y, z = u
    
    du[1] = σ * (y - x)
    du[2] = x * (ρ - z) - y
    du[3] = x * y - β * z
end

# Parameters for the Lorenz system
p = (10.0, 28.0, 8/3)
u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 1e3)

# Calculate the trajectory of the Lorenz system
prob = ODEProblem(lorenz!, u0, tspan, p)
trajectory = solve(prob, Tsit5(), saveat=0.01)

# Extract points from the trajectory
all_points = hcat(trajectory.u...)'
total_points = size(all_points, 1)

# Skip initial transient (first 10% of trajectory)
transient_skip = Int(ceil(0.1 * total_points))
points = all_points[(transient_skip+1):end, :]
num_points = size(points, 1)

println("Total trajectory points: $total_points")
println("Skipping initial transient: $transient_skip points")
println("Using $num_points points for analysis")

# Perform k-means clustering to create cellular complex
k = 128  # Number of clusters
kmeans_result = kmeans(points', k)
cluster_centers = kmeans_result.centers'
cluster_assignments = assignments(kmeans_result)

# Function to compute density estimate at each centroid
function compute_density(centers, points, assignments)
    density = zeros(size(centers, 1))
    
    # Count points in each cluster
    for a in assignments
        density[a] += 1
    end
    
    return density
end

# Compute point density at each centroid
centroid_density = compute_density(cluster_centers, points, cluster_assignments)

# Create a graph representation of the cellular complex
# Connect centers if there are consecutive points in different clusters
edges = Set{Tuple{Int, Int}}()

for i in 1:(num_points-1)
    c1, c2 = cluster_assignments[i], cluster_assignments[i+1]
    if c1 != c2
        edge = c1 < c2 ? (c1, c2) : (c2, c1)
        push!(edges, edge)
    end
end

# Convert edges to array
edges_array = collect(edges)

# Build the adjacency list
adjacency_list = Dict{Int, Set{Int}}()
for (i, j) in edges_array
    if !haskey(adjacency_list, i)
        adjacency_list[i] = Set{Int}()
    end
    if !haskey(adjacency_list, j)
        adjacency_list[j] = Set{Int}()
    end
    push!(adjacency_list[i], j)
    push!(adjacency_list[j], i)
end

# Create a KD-tree for efficient nearest neighbor search
kdtree = KDTree(cluster_centers')

# Find all 2-cells using a nearest neighbors approach
function find_all_2cells(centers, kdtree, adjacency_list, n_neighbors=6)
    n_centers = size(centers, 1)
    all_2cells = Vector{Vector{Int}}()
    
    for i in 1:n_centers
        # Find nearest neighbors of this center
        idxs, dists = knn(kdtree, centers[i, :], n_neighbors+1, true)
        neighbors = idxs[2:end]  # Skip self
        
        # Find triangles involving this center and its neighbors
        for j in 1:length(neighbors)
            for l in (j+1):length(neighbors)
                n1, n2 = neighbors[j], neighbors[l]
                
                # Check if these three points form a valid triangle
                # Both neighbors must be connected to the center in the graph
                if n1 in get(adjacency_list, i, Set{Int}()) && 
                   n2 in get(adjacency_list, i, Set{Int}())
                    # And the neighbors should be connected to each other
                    if n2 in get(adjacency_list, n1, Set{Int}())
                        push!(all_2cells, sort([i, n1, n2]))
                    end
                end
            end
        end
    end
    
    # Remove duplicates
    unique_2cells = unique(all_2cells)
    return unique_2cells
end

# Find 2-cells where each centroid participates in at least one
all_2cells = find_all_2cells(cluster_centers, kdtree, adjacency_list, 8)
println("Found $(length(all_2cells)) 2-cells from nearest neighbor approach")

# For centroids not in any 2-cell, create artificial 2-cells with nearest neighbors
centroids_in_2cells = Set(vcat(all_2cells...))
isolated_centroids = setdiff(Set(1:k), centroids_in_2cells)

if !isempty(isolated_centroids)
    println("Found $(length(isolated_centroids)) isolated centroids, creating artificial 2-cells")
    
    for i in isolated_centroids
        idxs, dists = knn(kdtree, cluster_centers[i, :], 3, true)
        # Create a triangle with the two closest neighbors
        push!(all_2cells, sort([i, idxs[2], idxs[3]]))
    end
    
    println("After adding artificial 2-cells: $(length(all_2cells)) 2-cells total")
end

# Count how many faces each edge is part of
edge_face_count = Dict{Tuple{Int, Int}, Int}()
for cell in all_2cells
    n = length(cell)
    for i in 1:n
        j = (i % n) + 1  # Next vertex in the cell
        edge = cell[i] < cell[j] ? (cell[i], cell[j]) : (cell[j], cell[i])
        edge_face_count[edge] = get(edge_face_count, edge, 0) + 1
    end
end

# Identify joining locus with a more strict threshold (3+ attached faces)
joining_locus = Set{Tuple{Int, Int}}()
for edge in edges_array
    if get(edge_face_count, edge, 0) >= 3
        push!(joining_locus, edge)
    end
end

println("Identified $(length(joining_locus)) edges in the joining locus ($(round(100*length(joining_locus)/length(edges_array), digits=1))% of all edges)")

# Print distribution of face counts
face_counts = collect(values(edge_face_count))
if !isempty(face_counts)
    println("Face count distribution: min=$(minimum(face_counts)), max=$(maximum(face_counts)), mean=$(mean(face_counts))")
    countmap = Dict{Int, Int}()
    for count in face_counts
        countmap[count] = get(countmap, count, 0) + 1
    end
    for count in sort(collect(keys(countmap)))
        percentage = round(100 * countmap[count] / length(edges_array), digits=1)
        println("  $count faces: $(countmap[count]) edges ($(percentage)%)")
    end
end

begin
  # Create the 3D visualization
  fig = Figure(resolution=(1200, 900))
  ax = Axis3(fig[1, 1], aspect=:data)

  # Plot the trajectory as a transparent line (excluding transient)
  # traj_plot = lines!(ax, points[:, 1], points[:, 2], points[:, 3], 
  #       linewidth=1.5, color=(:blue, 0.3), label="Lorenz Trajectory")

  # Plot the cluster centers
  centers_plot = scatter!(ax, cluster_centers[:, 1], cluster_centers[:, 2], cluster_centers[:, 3], 
          color=:red, markersize=15, label="Cluster Centers")

  # Plot the regular edges connecting clusters
  regular_edges = setdiff(edges_array, joining_locus)
  
  # Initialize plot elements with default values
  reg_edges_plot = nothing
  join_locus_plot = nothing
  
  # Plot the regular edges
  if !isempty(regular_edges)
      first_edge = first(regular_edges)
      i, j = first_edge
      
      reg_edges_plot = lines!(ax, [cluster_centers[i, 1], cluster_centers[j, 1]],
                 [cluster_centers[i, 2], cluster_centers[j, 2]],
                 [cluster_centers[i, 3], cluster_centers[j, 3]],
                 color=(:green, 0.7), linewidth=2.5, label="Regular Edges")
      
      # Plot the rest of the regular edges
      for edge in Iterators.drop(regular_edges, 1)
          i, j = edge
          lines!(ax, [cluster_centers[i, 1], cluster_centers[j, 1]],
                  [cluster_centers[i, 2], cluster_centers[j, 2]],
                  [cluster_centers[i, 3], cluster_centers[j, 3]],
                  color=(:green, 0.7), linewidth=2.5)
      end
  else
      # Create dummy plot for legend if no regular edges
      reg_edges_plot = lines!(ax, [0, 0], [0, 0], [0, 0], 
                            color=(:green, 0.7), linewidth=2.5, 
                            label="Regular Edges", visible=false)
  end
  
  # Plot the joining locus (edges with 3+ attached faces) with different color and thickness
  if !isempty(joining_locus)
      first_edge = first(joining_locus)
      i, j = first_edge
      
      join_locus_plot = lines!(ax, [cluster_centers[i, 1], cluster_centers[j, 1]],
                 [cluster_centers[i, 2], cluster_centers[j, 2]],
                 [cluster_centers[i, 3], cluster_centers[j, 3]],
                 color=(:magenta, 0.9), linewidth=4.0, label="Joining Locus (3+ faces)")
      
      # Plot the rest of the joining locus edges
      for edge in Iterators.drop(joining_locus, 1)
          i, j = edge
          lines!(ax, [cluster_centers[i, 1], cluster_centers[j, 1]],
                  [cluster_centers[i, 2], cluster_centers[j, 2]],
                  [cluster_centers[i, 3], cluster_centers[j, 3]],
                  color=(:magenta, 0.9), linewidth=4.0)
      end
  else
      # If no joining locus is found, create a dummy plot for the legend
      join_locus_plot = lines!(ax, [0, 0], [0, 0], [0, 0], 
                             color=(:magenta, 0.9), linewidth=4.0, 
                             label="Joining Locus (3+ faces)", visible=false)
  end

  # Create mesh representations of the 2-cells (faces)
  face_meshes = []

  # Flag to track visibility of 2-cells
  faces_visible = Observable(false)

  # Build mesh for each 2-cell
  if !isempty(all_2cells)
      for cell in all_2cells
          # Create polygon vertices
          polygon_points = [Point3f(cluster_centers[i, 1], cluster_centers[i, 2], cluster_centers[i, 3]) for i in cell]
          
          # For polygons with more than 3 vertices, we need to triangulate
          if length(polygon_points) == 3
              # Triangle - render directly
              mesh_face = mesh!(ax, polygon_points, color=(:cyan, 0.6), 
                            visible=faces_visible, shading=:NoShading)
              push!(face_meshes, mesh_face)
          else
              # Polygon with more than 3 vertices - triangulate using the fan method
              center = sum(polygon_points) / length(polygon_points)
              for i in 1:length(polygon_points)
                  triangle_points = [
                      center,
                      polygon_points[i],
                      polygon_points[i % length(polygon_points) + 1]
                  ]
                  mesh_face = mesh!(ax, triangle_points, color=(:cyan, 0.6), 
                                visible=faces_visible, shading=:NoShading)
                  push!(face_meshes, mesh_face)
              end
          end
      end
      
      # Add a mesh to the legend
      mesh_legend = mesh!(ax, [Point3f(0, 0, 0), Point3f(0, 0, 0), Point3f(0, 0, 0)], 
                         color=(:cyan, 0.6), visible=false, label="2-Cells (Faces)")
  end

  # Add a text label showing the current state
  toggle_label = Observable("2-Cells: OFF")
  
  # Add a label for the 2-cell visibility state
  Label(fig[1, 1, Bottom()], toggle_label, 
         fontsize=20, padding=(0, 0, 10, 0))
  
  # Create a button to toggle face visibility
  toggle_btn = Button(fig[1, 1, TopRight()], label="Toggle 2-Cells")
  
  # Set up the callback for the toggle button
  on(toggle_btn.clicks) do n
      faces_visible[] = !faces_visible[]
      toggle_label[] = faces_visible[] ? "2-Cells: ON" : "2-Cells: OFF"
  end
  
  # Also allow keyboard shortcut (press 't' to toggle)
  on(events(fig).keyboardbutton) do event
      if event.action == Keyboard.press && event.key == Keyboard.t
          faces_visible[] = !faces_visible[]
          toggle_label[] = faces_visible[] ? "2-Cells: ON" : "2-Cells: OFF"
      end
      return Consume(false)
  end

  # Add a legend - use the plot objects directly
  axislegend(ax)

  # Set axis labels
  ax.xlabel = "X"
  ax.ylabel = "Y"
  ax.zlabel = "Z"
  ax.title = "Lorenz System with Cellular Complex (k=$k, skipped $transient_skip transient points)"

  # Display the figure
  display(fig)
end