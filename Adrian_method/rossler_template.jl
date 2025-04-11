using Pkg
Pkg.activate(".")
Pkg.instantiate()

using GLMakie, OrdinaryDiffEq, Statistics, DataStructures
using LinearAlgebra, KernelDensity, NearestNeighbors, GeometryBasics

include("template/template.jl")
using .Template

include("fast_mcb.jl")
using .FastMCB

# Define the Rössler system
function rossler!(du, u, p, t)
    a, b, c = p
    x, y, z = u
    
    du[1] = -y - z
    du[2] = x + a * y
    du[3] = b + z * (x - c)
end

# Parameters for the Rössler system (chaotic)
a = 0.2
b = 0.2
c = 5.7
p = (a, b, c)
u0 = [0.0, 1.0, 0.0] # Typical initial condition
tspan = (0.0, 1e5) # Longer time span might be needed for good attractor coverage

# Calculate the trajectory of the Rössler system
prob = ODEProblem(rossler!, u0, tspan, p)
trajectory = solve(prob, Tsit5(), saveat=0.03, abstol=1e-8, reltol=1e-8) # Use a slightly larger saveat and tolerance

# Extract points from the trajectory
all_points = hcat(trajectory.u...)'
total_points = size(all_points, 1)

# Skip initial transient (adjust percentage if needed, e.g., 20%)
transient_skip = Int(ceil(0.2 * total_points))
points = all_points[(transient_skip+1):end, :]
num_points = size(points, 1)

println("Total trajectory points: $total_points")
println("Skipping initial transient: $transient_skip points")
println("Using $num_points points for analysis")

# Perform clustering using the template module.
k = 64  # Number of clusters (may need adjustment for Rössler)
@time clusters, assignments = Template.cluster(points, k)

# Extract centroids from clustering results.
centroids = [c[1] for c in clusters] # This is now a Vector{Vector{Float64}} of 3D points.
cluster_indices = [c[2] for c in clusters]

# Get the connectivity matrix.
connectivity_matrix = Template.cluster_connectivity(assignments, k)
undirected_connectivity_matrix = (connectivity_matrix + connectivity_matrix') # Note: This is unused later

# Set threshold for connectivity between clusters & get edges.
connectivity_threshold = 50 # May need adjustment for Rössler
edges = Template.edges(connectivity_matrix, connectivity_threshold)

# For debugging. If you get bidirectional connections, increase the
# connectivity threshold.
begin
  # Check for bidirectional connections in the connectivity matrix.
  bidirectional_connections = []
  for i in 1:size(connectivity_matrix, 1)
      for j in (i+1):size(connectivity_matrix, 2)
          # Check using the threshold defined for edges
          if connectivity_matrix[i,j] >= connectivity_threshold && connectivity_matrix[j,i] >= connectivity_threshold
              # Ensure the edge (i, j) or (j, i) was actually selected by Template.edges
              # This check is complex as Template.edges favors one direction
              # A simpler check: are both (i,j) and (j,i) above threshold?
               push!(bidirectional_connections, (i,j))
          end
      end
  end

  # Print information about bidirectional connections.
  if !isempty(bidirectional_connections)
      println("Found potential bidirectional connections (both directions >= threshold) between clusters:")
      for (i,j) in bidirectional_connections
          println("  Cluster $i ↔ Cluster $j (Check if both (i,j) and (j,i) are in 'edges')")
      end
  else
      println("No potential bidirectional connections found between clusters.")
  end
end

include("template/networkx.jl")
using .NetworkX

nx_graph = NetworkX.edges_to_networkx(edges, centroids)
mcb = NetworkX.mcb(nx_graph)

# Get the minimal cycle basis (vertex sequences).
mcb = FastMCB.fast_mcb(k, edges, centroids)

# Filter for 2-cells (faces).
two_cells = Template.faces(mcb, edges)
println("Found $(length(two_cells)) 2-cells (raw) out of $(length(mcb)) basis cycles.")

# Get unique faces before counting for joining locus
distinct_face_sequences = unique(two_cells) # Renamed variable
println("Found $(length(distinct_face_sequences)) unique 2-cell sequences.") # Use new name

# Calculate the joining locus using unique faces.
joining_locus_edges = Template.joining_locus(edges, distinct_face_sequences) # Pass renamed variable
joining_locus_set = Set(joining_locus_edges)
println("Found $(length(joining_locus_edges)) directed edges in the joining locus.")

# Plot.
begin
  fig = Figure()

  # Layout: Axis on top, controls below
  ax = Axis3(fig[1, 1], 
    xlabel = "x", 
    ylabel = "y", 
    zlabel = "z",
    title = "Rössler Template") # Updated title
  
  controls_grid = GridLayout(fig[2, 1], tellwidth = false)

  # Slider for selecting MCB cycle (to highlight edges)
  sl_label = Label(controls_grid[1, 1], "MCB Cycle Index:")
  sl = Slider(controls_grid[1, 2], range = 0:length(mcb), startvalue = 0)
  # Display slider value
  sl_val_label = Label(controls_grid[1, 3], string(sl.value[]), justification = :left)
  on(sl.value) do val
      sl_val_label.text = string(val)
  end

  # Toggle for face visibility
  toggle_label = Label(controls_grid[2, 1], "Show 2-Cell Faces:")
  face_toggle = Toggle(controls_grid[2, 2], active = true, halign = :left) # Start visible, align left
  faces_visible = face_toggle.active

  # Adjust column widths: Labels auto-size, slider/toggle column expands
  colsize!(controls_grid, 1, Auto())
  colsize!(controls_grid, 2, Relative(0.5)) # Make the middle column expand
  colsize!(controls_grid, 3, Auto())


  # --- Plot static elements --- 
  scatter!(ax, [c[1] for c in centroids], [c[2] for c in centroids], [c[3] for c in centroids], 
      color = :green, markersize = 16, label = "Centroids")

  # Prepare and plot edge segments and directional markers.
  black_segment_points = Vector{Point3f}() # For regular edges
  green_segment_points = Vector{Point3f}() # For joining locus edges
  red_marker_points = Vector{Point3f}()
  blue_marker_points = Vector{Point3f}()

  for edge_tuple in edges
      u, v = edge_tuple
      centroid_u = centroids[u]
      centroid_v = centroids[v]
      point_start = Point3f(centroid_u)
      point_end = Point3f(centroid_v)
      
      # Add to appropriate segment list based on joining locus membership
      if edge_tuple in joining_locus_set
          push!(green_segment_points, point_start)
          push!(green_segment_points, point_end)
      else
          push!(black_segment_points, point_start)
          push!(black_segment_points, point_end)
      end
      
      # Directional markers are plotted for all edges
      point_1_3 = centroid_u + (1/3) * (centroid_v - centroid_u)
      point_2_3 = centroid_u + (2/3) * (centroid_v - centroid_u)
      push!(red_marker_points, Point3f(point_1_3))
      push!(blue_marker_points, Point3f(point_2_3))
  end

  # Plot black edges (not in joining locus)
  if !isempty(black_segment_points)
      linesegments!(ax, black_segment_points, color = :black, linewidth = 1, label = "Edges")
  end
  # Plot green edges (in joining locus)
  if !isempty(green_segment_points)
      linesegments!(ax, green_segment_points, color = :red, linewidth = 5, label = "Joining Locus") # Thicker line
  end
  
  # Plot directional markers (remain the same)
  if !isempty(red_marker_points)
      scatter!(ax, red_marker_points, color = :red, markersize = 8, label = "Edge Start (+)")
  end
  if !isempty(blue_marker_points)
      scatter!(ax, blue_marker_points, color = :blue, markersize = 8, label = "Edge End (-)")
  end

  # --- Plot all 2-cell faces (meshes) --- 
  face_plots = [] # Store mesh plot objects
  for cycle_vertices in distinct_face_sequences # Use unique faces
      if length(cycle_vertices) >= 3
          pts = [Point3f(centroids[v]) for v in cycle_vertices]
          num_pts = length(pts) # Use local pts variable
          triangles = TriangleFace{Int}[] # Use TriangleFace
          p1_idx = 1
          for i in 2:(num_pts - 1)
              push!(triangles, TriangleFace(p1_idx, i, i + 1)) # Use TriangleFace
          end
          # Plot the mesh and store the plot object
          p = mesh!(ax, pts, triangles, color = (:cyan, 0.6), label=nothing) # Add label=nothing
          push!(face_plots, p)
      end
  end

  # Connect face visibility toggle to the mesh plots
  on(faces_visible) do is_visible
      for p in face_plots
          p.visible = is_visible
      end
  end
  # Set initial visibility
  notify(faces_visible)

  # --- Plot dynamic highlighted MCB cycle (edges) --- 
  mcb_lines_points_obs = Observable(Point3f[]) # Observable for line points
  mcb_lines_color_obs = Observable(:purple) # Observable for line color

  mcb_lines_plot = lines!(ax, mcb_lines_points_obs, 
                          color = mcb_lines_color_obs, 
                          linewidth = 5, 
                          label = "Selected MCB Cycle")

  # Update lines plot based on slider
  on(sl.value) do idx
      line_points = Point3f[]
      is_valid = false
      color = :purple

      if idx >= 1 && idx <= length(mcb)
          cycle_vertices = mcb[idx]
          if length(cycle_vertices) >= 3
              pts = [Point3f(centroids[v]) for v in cycle_vertices] # Calculate local pts
              push!(pts, pts[1]) # Close the loop using local pts
              line_points = pts    # Assign local pts to observable input
              is_valid = true
              # Check if it's NOT a 2-cell (using your original logic for cycle_is_2cell)
              is_2cell = Template.cycle_is_2cell(cycle_vertices, edges)
              color = is_2cell ? :purple : :orange
          end
      end

      mcb_lines_points_obs[] = line_points
      mcb_lines_color_obs[] = color
      mcb_lines_plot.visible = is_valid
  end

  # Trigger initial update for lines plot
  notify(sl.value)

  axislegend(ax)
  display(fig)
end