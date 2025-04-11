# networkx.jl

module NetworkX

using LinearAlgebra
using PythonCall

# Import the networkx Python library
# This assumes networkx is installed in the Python environment PythonCall uses.
const nx = pyimport("networkx")

"""
    edges_to_networkx(edges)

Convert a list of edges to a NetworkX graph.

Args:
- edges: A vector of 2-tuples representing directed edges.

Returns:
- A NetworkX digraph object containing the edges.
"""
function edges_to_networkx(edges, centroids)
  # Convert the flow_graph to a networkx graph.
  nx_graph = nx.DiGraph()

  for edge in edges
    nx_graph.add_edge(
      edge[1],
      edge[2],
      length=norm(centroids[edge[2]] - centroids[edge[1]])
    )
  end

  return nx_graph
end

"""
    mcb(nx_graph)

Compute the minimum cycle basis of an undirected NetworkX graph.

Args:
- nx_graph: A NetworkX digraph object.

Returns:
- A list of cycles, where each cycle is represented as a list of vertices.
"""
function mcb(nx_graph)
  undirected_nx_graph = nx_graph.to_undirected()
  python_mcb = nx.minimum_cycle_basis(undirected_nx_graph, weight="length")
  # Convert Python list of lists to Julia Vector{Vector{Int}} using PythonCall.
  julia_mcb = pyconvert(Vector{Vector{Int}}, python_mcb)

  return julia_mcb
end

# Or just let users access functions via NetworkX.nx
# E.g., NetworkX.nx.Graph(), NetworkX.nx.cycle_basis(py_graph)

end # module NetworkX