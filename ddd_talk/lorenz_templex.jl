using Pkg
Pkg.activate(".")
using Colors
using DifferentialEquations
using DynamicalSystems
using GLMakie
using Graphs
using GraphPlot
using LinearAlgebra
using Random
using Statistics

# Set random seed for reproducibility.
Random.seed!(0)

##### Lorenz system solution.

# Define the Lorenz system parameters.
function lorenz!(du, u, p, t)
    σ, ρ, β = p
    x, y, z = u
    
    du[1] = σ * (y - x)
    du[2] = x * (ρ - z) - y
    du[3] = x * y - β * z
end

# Set up parameters and initial conditions.
p = [10.0, 28.0, 8/3]  # σ, ρ, β
u0 = [1.0, 0.0, 0.0] # Initial state.
tspan = (0.0, 4e2) # Time span.
transient_length = 1e2 # Length of transient to remove from solution.

# Create and solve the problem.
prob = ODEProblem(lorenz!, u0, tspan, p)
sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)

# Remove transient from solution points.
transient_i = findfirst(t -> t > transient_length, sol.t)
us = sol.u[transient_i:end]

# Create a 3D plot of the solution.
fig_lorenz_solve = Figure()
ax = Axis3(fig_lorenz_solve[1, 1], xlabel="x", ylabel="y", zlabel="z")

# Extract solution components.
xs = [u[1] for u in us]
ys = [u[2] for u in us]
zs = [u[3] for u in us]

# Plot the trajectory.
# lines!(ax, xs, ys, zs, color=:blue, linewidth=1)

# Display the figure.
display(fig_lorenz_solve)

##### BraMAH complex determination.

d = 2 # Expected dimension of the attractor. Note: We can actually get this from the Lyapunov dimension if we don't know it in advance. Remove this?
m_min = Int(round(length(us)/12)) # Minimum number of points in a patch.
m_max = Int(round(length(us)/3)) # Maximum number of points in a patch.

N = length(us[1]) # Dimension of the phase space.

remaining_points = Set(1:length(us)) # Track indices of points not yet in a patch
patches = [] # Store all patches
patch_centers = [] # Store centers of patches

# Compute optimal patch size.
function compute_m_0(X_maximal, m_min, m_max, d)
    m0_values = collect(m_min:m_max)
    num_m0 = length(m0_values)
    sigma_values = zeros(num_m0, d)
    
    # Collect singular values for each m_0
    for (idx, m_0) in enumerate(m0_values)
        # Compute candidate X matrix.
        X = X_maximal[1:m_0, :] / sqrt(m_0)
        # Compute SVD of X matrix.
        F = svd(X)
        # Store the first d singular values.
        sigma_values[idx, :] = F.S[1:d]
    end
    
    best_m_0 = m_min
    best_R2_min = -Inf
    
    # Perform linear regression for each possible m_0
    for idx in 1:num_m0
        m0_subset = m0_values[1:idx]
        sigma_subset = sigma_values[1:idx, :]
        R2_values = zeros(d)
        for j = 1:d
            y = sigma_subset[:, j]
            x = m0_subset
            n = length(x)
            x_mean = mean(x)
            y_mean = mean(y)
            s_xy = sum((x .- x_mean) .* (y .- y_mean))
            s_xx = sum((x .- x_mean).^2)
            # Compute slope and intercept
            slope = s_xy / s_xx
            intercept = y_mean - slope * x_mean
            # Predicted values
            y_pred = slope .* x .+ intercept
            # Compute R²
            ss_tot = sum((y .- y_mean).^2)
            ss_res = sum((y .- y_pred).^2)
            R2 = 1 - ss_res / ss_tot
            R2_values[j] = R2
        end
        # Find the minimum R² among the d singular values
        R2_min = minimum(R2_values)
        # Update the best m_0 if R² is improved
        if R2_min > best_R2_min
            best_R2_min = R2_min
            best_m_0 = m0_values[idx]
        end
    end
    return best_m_0, best_R2_min
end

# Compute patches.
while !isempty(remaining_points)
  println("Computing patch #", length(patches)+1, "...")
  # Select random point from remaining points as patch center.
  center_idx = rand(remaining_points)
  x_p = us[center_idx]
  push!(patch_centers, x_p)
  
  # Compute distances from x_p to all remaining points.
  distances = [norm(x_p - u) for u in us]
  
  # Sort remaining points by distance from x_p.
  distance_order = sortperm(distances)
  sorted_points = us[distance_order]
  
  # Take up to m_max closest points.
  patch_size = min(m_max, length(sorted_points))
  
  # Form X matrix from patch points.
  X_maximal = reduce(vcat, transpose.(sorted_points[1:patch_size]))
  
  # Find optimal patch size.
  best_m_0, best_SSE = compute_m_0(X_maximal, m_min, patch_size, d)
  
  # Store final patch points.
  push!(patches, distance_order[1:best_m_0])
  
  # Remove patch points from remaining points.
  setdiff!(remaining_points, distance_order[1:best_m_0])

  # Print percentage of points processed
  total_points = length(us)
  points_processed = total_points - length(remaining_points)
  percentage = round(100 * points_processed / total_points, digits=1)
  println("Processed ", percentage, "% of points.")
end

# Find boundaries of d-cells from overlaps between patches.
cell_boundaries = [[] for _ in 1:length(patches)]
for (i, patch_indices) in enumerate(patches)
  patch = us[patch_indices]
  for j in i+1:length(patches)
    intersection = intersect(patch_indices, patches[j])
    # Check if patch meets with another with greater index.
    if length(intersection) >= d
      # Select d random representative points from the intersection to form a face.
      representative_points = rand(intersection, d)
      # Introduce the face to the faces vector for each d-cell.
      cell_boundaries[i] = [cell_boundaries[i]..., representative_points...]
      cell_boundaries[j] = [cell_boundaries[j]..., representative_points...]
    end
  end
end

cell_perimeters = [[] for _ in 1:length(patches)]
for (i, cell_boundary) in enumerate(cell_boundaries)
  # Skip if cell boundary is empty.
  if isempty(cell_boundary)
    continue
  end

  # Perform SVD on centered patch points.
  patch_center = patch_centers[i]
  patch_points = us[patches[i]]
  X = reduce(vcat, transpose.(patch_points .- Ref(patch_center))) / sqrt(length(patch_points))
  F = svd(X)

  # Take first d singular vectors as tangent plane basis.
  tangent_basis = F.V[:, 1:d]
  
  # Project boundary points onto tangent plane.
  boundary_points = us[cell_boundary]
  centered_boundary = boundary_points .- Ref(patch_center)
  projected_points = [tangent_basis' * point for point in centered_boundary]

  # Re-order projected points in counter-clockwise order.
  projected_points_ordering = sortperm(projected_points, by=x -> atan(x[2], x[1]))
  
  # Store projected points.
  cell_perimeters[i] = boundary_points[projected_points_ordering]
end

# TODO: Finish constructing cellular complex.

### Alternative begin.
cells = []
for (i, patch_idxs) in enumerate(patches)
  cell = []
  for point_idx in patch_idxs
    # Check if point_idx is not in any cell.
    if !any(cell_points -> point_idx in cell_points, cells)
      push!(cell, point_idx)
    end
  end
  push!(cells, cell)
end
### Alternative end.

# Calculate patch-flow adjacency matrix.
flow_matrix = zeros(length(cells), length(cells))
for (i, cell_indices) in enumerate(cells)
  # Iterate over each point in the cell.
  for j in cell_indices
    # Is the forward flow of the point still within the cell?
    if !(j+1 in cell_indices)
      # Determine which cells the point flows into.
      for (k, kth_cell_indices) in enumerate(cells)
        if j+1 in kth_cell_indices
          # Add edge to adjacency matrix.
          flow_matrix[i, k] += 1
        end
      end
    end
  end
end

# Compute the joining locus from the flow, ignoring any topological information.
# incoming_traj_min = 210 # Minimum number of incoming cells needed to be considered a joining cell.
# joining_locus = findall(col -> sum(col) >= incoming_traj_min, eachcol(flow_matrix))

# Plot all cells with random colors.
for (i, cell_idxs) in enumerate(cells)
  cell = us[cell_idxs]
  color = RGB(rand(), rand(), rand())
  # scatter!(
  #   ax,
  #   [u[1] for u in cell],
  #   [u[2] for u in cell],
  #   [u[3] for u in cell],
  #   color=color,
  #   markersize=10,
  #   alpha=0.4
  # )
  
  # Plot patch center with black outline.
  center = patch_centers[i]
  scatter!(
    ax,
    [center[1]],
    [center[2]], 
    [center[3]],
    color=:black,
    markersize=15
  )

  # Plot patch index on patch center.
  text!(
    ax,
    [center[1]],
    [center[2]],
    [center[3]],
    text=string(i),
    color=color,
    fontsize=48
  )

  # # Plot cell.
  # # Extract x,y,z coordinates of perimeter points.
  cell_xs = [point[1] for point in cell_perimeters[i]]
  cell_ys = [point[2] for point in cell_perimeters[i]]
  cell_zs = [point[3] for point in cell_perimeters[i]]
  # Plot the cell boundary as a closed polygon.
  # Add first point again to close the loop.
  push!(cell_xs, cell_xs[1])
  push!(cell_ys, cell_ys[1]) 
  push!(cell_zs, cell_zs[1])
  lines!(ax, cell_xs, cell_ys, cell_zs, color=color, linewidth=2)
  # Create triangulation for mesh
  n = length(cell_xs)
  center_x = mean(cell_xs)
  center_y = mean(cell_ys) 
  center_z = mean(cell_zs)
  # Draw lines between consecutive perimeter points
  for j in 1:n-1
    lines!(ax, [cell_xs[j], cell_xs[j+1]], [cell_ys[j], cell_ys[j+1]], [cell_zs[j], cell_zs[j+1]], color=color, linewidth=2)
  end
end
display(fig_lorenz_solve)
