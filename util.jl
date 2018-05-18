using CSV
using Graphs
#############################################################################################################

#This function converts an adjacency matrix into a directed graph.
function adjacencyMatrixToDigraph(A::SparseMatrixCSC{Float64, Int64})

  n = size(A,1)

  g = inclist(collect(1:n),is_directed=true)
  
  edges = findnz(A)
  
  m = length(edges[1])

  capacities = zeros(m)
  
  for i = 1:m
    add_edge!(g,edges[1][i],edges[2][i])
    capacities[i] = edges[3][i]
  end
 
  return g,capacities
 
end


#############################################################################################################

#This function outputs the adjacency matrix of the largest connected subgraph from
#the graph defined in the file in "path".
#It also deletes self-loops.
function extractGraphFromFile(path; separator='\t')
  # read in data and let indices start from 1
  #iris = readtable(string(path, "data.txt"), header = true, separator = separator, eltypes = [Int64, Int64])
  iris = CSV.read(path, delim = separator)
  iris = convert(Array{Float64,2}, iris)
  if size(iris,2)==2
    w = ones(length(iris[:,1]))
  else
    w = iris[:,3]
  end
  
  iris = convert(Array{Int64,2}, iris[:,1:2]+1)
  
  # transfer these edges into sparse adjacency matrix, weight of each edge is 1
  n = maximum(iris)
  mat = sparse(vcat(iris[:,1], iris[:,2]), vcat(iris[:, :2], iris[:, :1]), vcat(w, w), n, n, *)
  
  # get the largest connected subgraph in original graph.
  matCopy = adjacencyMatrixToDigraph(mat)[1]
  indices = strongly_connected_components(matCopy)
  sort!(indices, by = x -> length(x), rev = true)
  firstComponentIndices = sort!(indices[1])
  A = mat[firstComponentIndices, firstComponentIndices]
  
  # delete self-loops
  A = A - sparse(diagm(diag(A)))
  
  A = convert(SparseMatrixCSC{Float64,Int64}, A)
  return A
end

################################################################################
################################################################################
#Function to find a maximal independent set, given an adjacency matrix
function maximalIS(A::SparseMatrixCSC{Float64, Int64}; variant="mindeg")
  n = size(A,1)
  if n == 0
    return []
  else
    if variant=="random"
	  v = rand(1:n)
	elseif variant=="mindeg"
	  v = indmin(sum(A,1))
	end
	nbr_v = collect(findnz(A[v,:])[1])
	rem_v = setdiff( collect(1:n) , union(v , nbr_v))
	IS_rem = maximalIS(A[rem_v, rem_v], variant=variant)
    return union(v, rem_v[IS_rem])
  end
end
