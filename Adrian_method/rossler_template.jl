using Pkg
Pkg.activate(".")
Pkg.instantiate()

using GLMakie, OrdinaryDiffEq, Statistics, DataStructures
using LinearAlgebra, KernelDensity, NearestNeighbors, GeometryBasics

include("template/template.jl")
include("template/networkx.jl")
using .Template
using .NetworkX

# Define the Rössler system
function rossler!(du, u, p, t)
    a, b, c = p
    x, y, z = u
    
    du[1] = -y - z
    du[2] = x + a * y
    du[3] = b + z * (x - c)
end

# Parameters for the Rössler system (chaotic)
a = 0.38
b = 0.2
c = 4.822
p = (a, b, c)
u0 = [1.0, 1.0, 1.0] # Typical initial condition
tspan = (0.0, 1.5e4) # Longer time span might be needed for good attractor coverage

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

# --- Prepare subset for plotting last 10% of *original* trajectory --- 
# Calculate number of points for the last 10% of the full trajectory
num_plot_points = min(10000, total_points)
# Calculate the starting index within the 'points' array (post-transient)
# to show the equivalent duration at the end.
plot_start_idx_in_points = max(1, num_points - num_plot_points + 1)

# Create the subset of points and assignments to plot
plot_points_subset = points[plot_start_idx_in_points:end, :]

println("Total trajectory points: $total_points")
println("Skipping initial transient: $transient_skip points")
println("Using $num_points points for analysis")
println("Plotting last segment equivalent to $num_plot_points points (from index $plot_start_idx_in_points in analysis points)")

# Perform clustering using the template module.
k = 80  # Number of clusters (may need adjustment for Rössler)
connectivity_threshold = 20 # May need adjustment for Rössler
max_force_boundary = 4
rounds = 10

# Use best_complex to optimize the template construction.
best = Template.best_complex(
  points,
  k;
  rounds=rounds,
  connectivity_threshold=connectivity_threshold,
  max_force_boundary=max_force_boundary,
  print_frequencies=true
)

println("Best joining locus length: $(best.joining_locus_length) (found in round $(best.round))")

# Get the assignments corresponding to the plotted subset
plot_assignments_subset = best.assignments[plot_start_idx_in_points:end]

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
  sl = Slider(controls_grid[1, 2], range = 0:length(best.mcb), startvalue = 0)
  # Display slider value
  sl_val_label = Label(controls_grid[1, 3], string(sl.value[]), justification = :left)
  on(sl.value) do val
      sl_val_label.text = string(val)
  end

  # Toggle for face visibility
  toggle_label = Label(controls_grid[2, 1], "Show 2-Cell Faces:")
  face_toggle = Toggle(controls_grid[2, 2], active = true, halign = :left) # Start visible, align left
  faces_visible = face_toggle.active

  # Toggle for trajectory points visibility
  traj_toggle_label = Label(controls_grid[3, 1], "Show Trajectory Points:")
  traj_toggle = Toggle(controls_grid[3, 2], active = false, halign = :left) # Start INvisible, align left
  traj_visible = traj_toggle.active

  # Adjust column widths: Labels auto-size, slider/toggle column expands
  colsize!(controls_grid, 1, Auto())
  colsize!(controls_grid, 2, Relative(0.5)) # Make the middle column expand
  colsize!(controls_grid, 3, Auto())

#   sl = controls_grid.sliders[1]

  # Define placeholder geometry for when no cycle is shown
  placeholder_point = isempty(best.centroids) ? Point3f(0,0,0) : Point3f(best.centroids[1])
  placeholder_points = fill(placeholder_point, 3) 
  placeholder_triangles = [GeometryBasics.TriangleFace(1, 2, 3)]

  # -- Plot dynamic MCB cycle mesh --
  # Plot an initial, invisible placeholder mesh
  mcb_mesh = mesh!(ax, placeholder_points, placeholder_triangles,
      color = (:purple, 0.5), 
      label = "Selected MCB Cycle",
      visible = false) # Start invisible

  # Update mesh plot attributes when slider changes
  on(sl.value) do idx
      is_valid_cycle = false
      
      if idx >= 1 && idx <= length(best.mcb)
          cycle_vertices_indices = best.mcb[idx]
          if length(cycle_vertices_indices) >= 3 # Need at least 3 vertices for a face
              # Get centroid points for the cycle vertices.
              cycle_vertices_geom = [Point3f(best.centroids[v]) for v in cycle_vertices_indices]
              
              # Simple triangulation: Fan from the first vertex.
              num_pts = length(cycle_vertices_geom)
              triangles = GeometryBasics.TriangleFace{Int}[] 
              p1_idx = 1 
              for i in 2:(num_pts - 1)
                  p2_idx = i
                  p3_idx = i + 1
                  push!(triangles, GeometryBasics.TriangleFace(p1_idx, p2_idx, p3_idx))
              end
              
              # Check if it fails the 2-cell criteria (your original logic)
              is_not_2cell = Template.cycle_is_2cell(cycle_vertices_indices, best.edges)
              color_tuple = is_not_2cell ? (:orange, 0.5) : (:purple, 0.5)

              # Update plot data and attributes
              mcb_mesh.vertices = cycle_vertices_geom
              mcb_mesh.faces[] = triangles # Update faces observable/attribute
              mcb_mesh.color = color_tuple
              is_valid_cycle = true
          end
      end

      # Update visibility
      mcb_mesh.visible = is_valid_cycle
  end

  # Trigger initial update based on startvalue
  notify(sl.value)

  # --- Plot static elements --- 
  # Plot Trajectory Points (last segment, controlled by toggle)
  scatter!(ax, plot_points_subset[:, 1], plot_points_subset[:, 2], plot_points_subset[:, 3], 
      color = plot_assignments_subset, # Use subset assignments
      colormap = :turbo, 
      colorrange = (1, k), # Explicitly set the color range
      markersize = 2, 
      label = "Trajectory Points (Last $num_plot_points)", # Updated label
      visible = traj_visible) # Link visibility to toggle observable
      
  # Plot Trajectory Path (last segment, controlled by same toggle)
  lines!(ax, plot_points_subset[:, 1], plot_points_subset[:, 2], plot_points_subset[:, 3], # Use subset points
      color = (:black, 0.3), # Faint black line
      linewidth = 1,
      label = nothing, # Don't add path to legend
      visible = traj_visible) # Link visibility to toggle observable
      
  # Plot Centroids (always visible)
  scatter!(ax, [c[1] for c in best.centroids], [c[2] for c in best.centroids], [c[3] for c in best.centroids], 
      color = :red, # Changed centroids to red for better contrast
      markersize = 10, # Slightly smaller centroids
      label = "Centroids")

  # Prepare and plot edge segments and directional markers.
  black_segment_points = Vector{Point3f}() # For regular edges
  green_segment_points = Vector{Point3f}() # For joining locus edges
  red_marker_points = Vector{Point3f}()
  blue_marker_points = Vector{Point3f}()

  for edge_tuple in best.edges
      u, v = edge_tuple
      centroid_u = best.centroids[u]
      centroid_v = best.centroids[v]
      point_start = Point3f(centroid_u)
      point_end = Point3f(centroid_v)
      
      # Add to appropriate segment list based on joining locus membership
      if edge_tuple in Set(best.joining_locus)
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
  for cycle_vertices in best.faces # Use unique faces
      if length(cycle_vertices) >= 3
          pts = [Point3f(best.centroids[v]) for v in cycle_vertices]
          num_pts = length(pts) 
          triangles = GeometryBasics.TriangleFace{Int}[] # Explicitly use GeometryBasics.TriangleFace
          p1_idx = 1
          for i in 2:(num_pts - 1)
              push!(triangles, GeometryBasics.TriangleFace(p1_idx, i, i + 1)) # Explicitly use GeometryBasics.TriangleFace
          end
          # Plot the mesh and store the plot object
          p = mesh!(ax, pts, triangles, color = (:cyan, 0.6), label=nothing)
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

  # Add a legend (optional, but helpful)
  axislegend(ax)
  display(fig)
end