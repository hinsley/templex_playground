using Pkg
Pkg.activate(".")
Pkg.instantiate()

using GLMakie, OrdinaryDiffEq, Statistics, DataStructures
using LinearAlgebra, KernelDensity, NearestNeighbors, GeometryBasics

include("template/template.jl")
include("template/networkx.jl")
using .Template
using .NetworkX

# Define the Shimizu-Morioka system
function shimizu_morioka!(du, u, p, t)
    alpha, lambda = p
    x, y, z = u
    du[1] = y
    du[2] = x - lambda * y - x * z
    du[3] = -alpha * z + x^2
end

# Parameters for the Shimizu-Morioka system (chaotic with lacuna)
alpha = 0.5
lambda = 0.85
p = (alpha, lambda)
u0 = [1.0, 1.0, 1.0]
tspan = (0.0, 1.5e4)

# Calculate the trajectory of the Shimizu-Morioka system
prob = ODEProblem(shimizu_morioka!, u0, tspan, p)
trajectory = solve(prob, Tsit5(), saveat=0.003, abstol=1e-8, reltol=1e-8)

# Extract points from the trajectory
all_points = hcat(trajectory.u...)'
total_points = size(all_points, 1)

# Skip initial transient (adjust percentage if needed, e.g., 20%)
transient_skip = Int(ceil(0.2 * total_points))
points = all_points[(transient_skip+1):end, :]
num_points = size(points, 1)

sparsification_len = 0.03
points = Template.sparsify(points, sparsification_len)
num_points = size(points, 1)

percent = 100 * num_points / (total_points - transient_skip)
println("Sparsified trajectory has $num_points points (" * string(round(percent, digits=1)) * "% of post-transient points)")

# --- Prepare subset for plotting last N points of the sparsified trajectory ---
num_plot_points = min(10000, num_points)
plot_start_idx_in_points = max(1, num_points - num_plot_points + 1)
plot_points_subset = points[plot_start_idx_in_points:end, :]

println("Total trajectory points: $total_points")
println("Skipping initial transient: $transient_skip points")
println("Using $num_points points for analysis")
println("Plotting last segment equivalent to $num_plot_points points (from index $plot_start_idx_in_points in analysis points)")

k = 80
connectivity_threshold = 10
max_force_boundary = 4
rounds = 1

best = Template.best_complex(
  points,
  k;
  rounds=rounds,
  connectivity_threshold=connectivity_threshold,
  max_force_boundary=max_force_boundary,
  print_frequencies=true
)

println("Best joining locus length: $(best.joining_locus_length) (found in round $(best.round))")

plot_assignments_subset = best.assignments[plot_start_idx_in_points:end]

begin
  fig = Figure()

  ax = Axis3(
    fig[1, 1], 
    xlabel = "x", 
    ylabel = "y", 
    zlabel = "z",
    title = "Shimizu-Morioka Template"
  )
  
  controls_grid = GridLayout(fig[2, 1], tellwidth = false)

  sl_label = Label(controls_grid[1, 1], "MCB Cycle Index:")
  sl = Slider(controls_grid[1, 2], range = 0:length(best.mcb), startvalue = 0)
  sl_val_label = Label(controls_grid[1, 3], string(sl.value[]), justification = :left)
  on(sl.value) do val
    sl_val_label.text = string(val)
  end

  toggle_label = Label(controls_grid[2, 1], "Show 2-Cell Faces:")
  face_toggle = Toggle(controls_grid[2, 2], active = true, halign = :left)
  faces_visible = face_toggle.active

  traj_toggle_label = Label(controls_grid[3, 1], "Show Trajectory Points:")
  traj_toggle = Toggle(controls_grid[3, 2], active = false, halign = :left)
  traj_visible = traj_toggle.active

  colsize!(controls_grid, 1, Auto())
  colsize!(controls_grid, 2, Relative(0.5))
  colsize!(controls_grid, 3, Auto())

  placeholder_point = isempty(best.centroids) ? Point3f(0,0,0) : Point3f(best.centroids[1])
  placeholder_points = fill(placeholder_point, 3) 
  placeholder_triangles = [GeometryBasics.TriangleFace(1, 2, 3)]

  mcb_mesh = mesh!(
    ax,
    placeholder_points,
    placeholder_triangles,
    color = (:purple, 0.5), 
    label = "Selected MCB Cycle",
    visible = false
  )

  on(sl.value) do idx
    is_valid_cycle = false
    
    if idx >= 1 && idx <= length(best.mcb)
      cycle_vertices_indices = best.mcb[idx]
      if length(cycle_vertices_indices) >= 3
        cycle_vertices_geom = [Point3f(best.centroids[v]) for v in cycle_vertices_indices]
        num_pts = length(cycle_vertices_geom)
        triangles = GeometryBasics.TriangleFace{Int}[] 
        p1_idx = 1 
        for i in 2:(num_pts - 1)
          p2_idx = i
          p3_idx = i + 1
          push!(triangles, GeometryBasics.TriangleFace(p1_idx, p2_idx, p3_idx))
        end
        is_not_2cell = Template.cycle_is_2cell(cycle_vertices_indices, best.edges)
        color_tuple = is_not_2cell ? (:orange, 0.5) : (:purple, 0.5)
        mcb_mesh.vertices = cycle_vertices_geom
        mcb_mesh.faces[] = triangles
        mcb_mesh.color = color_tuple
        is_valid_cycle = true
      end
    end
    mcb_mesh.visible = is_valid_cycle
  end

  notify(sl.value)

  scatter!(
    ax,
    plot_points_subset[:, 1],
    plot_points_subset[:, 2],
    plot_points_subset[:, 3], 
    color = plot_assignments_subset, 
    colormap = :turbo, 
    colorrange = (1, k), 
    markersize = 2, 
    label = "Trajectory Points (Last $num_plot_points)", 
    visible = traj_visible
  )
  lines!(
    ax,
    plot_points_subset[:, 1],
    plot_points_subset[:, 2],
    plot_points_subset[:, 3],
    color = (:black, 0.3), 
    linewidth = 1,
    label = nothing, 
    visible = traj_visible
  )
  scatter!(
    ax,
    [c[1] for c in best.centroids],
    [c[2] for c in best.centroids],
    [c[3] for c in best.centroids], 
    color = :green,
    markersize = 10,
    label = "Centroids"
  )

  black_segment_points = Vector{Point3f}()
  red_segment_points = Vector{Point3f}()
  red_marker_points = Vector{Point3f}()
  blue_marker_points = Vector{Point3f}()

  for edge_tuple in best.edges
    u, v = edge_tuple
    centroid_u = best.centroids[u]
    centroid_v = best.centroids[v]
    point_start = Point3f(centroid_u)
    point_end = Point3f(centroid_v)
    if edge_tuple in Set(best.joining_locus)
        push!(red_segment_points, point_start)
        push!(red_segment_points, point_end)
    else
        push!(black_segment_points, point_start)
        push!(black_segment_points, point_end)
    end
    point_1_3 = centroid_u + (1/3) * (centroid_v - centroid_u)
    point_2_3 = centroid_u + (2/3) * (centroid_v - centroid_u)
    push!(red_marker_points, Point3f(point_1_3))
    push!(blue_marker_points, Point3f(point_2_3))
  end

  if !isempty(black_segment_points)
    linesegments!(ax, black_segment_points, color = :black, linewidth = 1, label = "Edges")
  end
  if !isempty(red_segment_points)
    linesegments!(ax, red_segment_points, color = :red, linewidth = 5, label = "Joining Locus")
  end
  if !isempty(red_marker_points)
    scatter!(ax, red_marker_points, color = :red, markersize = 8, label = "Edge Start (+)")
  end
  if !isempty(blue_marker_points)
    scatter!(ax, blue_marker_points, color = :blue, markersize = 8, label = "Edge End (-)")
  end

  face_plots = []
  for cycle_vertices in best.faces
    if length(cycle_vertices) >= 3
      pts = [Point3f(best.centroids[v]) for v in cycle_vertices]
      num_pts = length(pts)
      triangles = GeometryBasics.TriangleFace{Int}[]
      p1_idx = 1
      for i in 2:(num_pts - 1)
        push!(triangles, GeometryBasics.TriangleFace(p1_idx, i, i + 1))
      end
      p = mesh!(ax, pts, triangles, color = (:cyan, 0.6), label=nothing)
      push!(face_plots, p)
    end
  end

  on(faces_visible) do is_visible
    for p in face_plots
      p.visible = is_visible
    end
  end
  notify(faces_visible)

  axislegend(ax)
  display(fig)
end