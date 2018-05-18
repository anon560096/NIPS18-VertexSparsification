using StatsBase
include("util.jl")
include("EdmondsKarp.jl")

global INF = 1000000000

#Wrapper function that adds a super-source and a super-sink and then runs min-st-cut
function minSTCut(a::SparseMatrixCSC{Float64, Int64}, terminal_set_S::Array{Int64,1}, terminal_set_T::Array{Int64,1})

  A = a[:,:]
  n = size(A,1)
  
  super_S_row = [INF*Int(i in terminal_set_S) for i in 1:n+2]
  super_T_row = [INF*Int(i in terminal_set_T) for i in 1:n+2]

  A = vcat(A,super_S_row[1:n]',super_T_row[1:n]')
  A = hcat(A,super_S_row,super_T_row)
  
  G,capacities = adjacencyMatrixToDigraph(A)
  
  
  parity,value = min_st_cut(G,n+1,n+2,capacities)
  
  if parity[n+1] == true
  	parity = .~parity
  end
  parity = convert(Array{Int8},parity+1)
  
  return parity[1:n], value
end

#Greedy 2(1-1/k) approx k-way cut
function minSTkWayCut(a::SparseMatrixCSC{Float64, Int64}, terminal_sets::Array{Array{Int64,1},1}, k::Int64)
  
  A = a[:,:]
  n = size(A,1)
  
  if k==2
  	return minSTCut(a,terminal_sets[1],terminal_sets[2])
  end
  
  parity_array = []
  value_array = []

  for i = 1:k
  	parity,value = minSTCut(A,terminal_sets[i],reduce(vcat, terminal_sets[setdiff(1:k,i)]))
  	append!(parity_array,[parity])
  	append!(value_array,[value])
  end
  


  maxvalue = maximum(value_array)
  value = sum(value_array) - maxvalue

  maxindex = findfirst(value_array.==maxvalue)

  parity = ones(n).*maxindex

  for i in setdiff(1:k,maxindex)
  	parity[parity_array[i].==1] = i
  end

  parity = convert(Array{Int8},parity)
  return parity, value

end

#Function to contract all nonterminals. Returns contracted graph,
#and for each vertex in the contracted graph, the list of vertices absorbed into it.
function contractVertices(a::SparseMatrixCSC{Float64, Int64}, terminals::Array{Int64,1}, sparsifier_size::Int64)
  
  A=a[:,:]
  
  n = size(A,1)
  
  absorbed_vertices = sparse(eye(Bool,n)) #intitially each vertex has absorbed only itself
  
  nonterminals = setdiff((1:n), terminals)

  removable_vertices = sample(nonterminals, n - sparsifier_size, replace=false, ordered=true)
  remaining_vertices = setdiff((1:n), removable_vertices)
  
  for v in removable_vertices
    #pick a neighbor (absorber_v) with probability proportional to edge weights
    nbr_v, wt_v = findnz(A[v,:])
	absorber_v = sample(nbr_v, Weights(wt_v))
	
	#Do the absorption process
	A[absorber_v,:] += A[v,:]
	A[:,absorber_v] += A[:,v]
	A[:,v] = 0
	A[v,:] = 0
	A[absorber_v,absorber_v] = 0
	
	#Transfer v's absorbed vertices to absorber_v
	absorbed_vertices[absorber_v, findnz(absorbed_vertices[v,:])[1]] = true
	absorbed_vertices[v,:] = false
  end
  
  A=dropzeros(A)
  absorbed_vertices=dropzeros(absorbed_vertices)
	
  return A, absorbed_vertices
end

function contractMinSTkWayCut(A::SparseMatrixCSC{Float64,Int64}, terminal_sets::Array{Array{Int64,1},1}, k::Int64, sparsifier_size::Int64)
  
  n = size(A,1)

  terminals = reduce(vcat,terminal_sets)
  
  #Contract enough nonterminals
  B, absorbed_vertices = contractVertices(A, terminals, sparsifier_size)
  
  contracted_parity, contracted_value = minSTkWayCut(B,terminal_sets,k)

  contracted_sets = []
  Aparity = zeros(Int8, n)
 
  for i = 1:k
  	append!(contracted_sets, [findnz(absorbed_vertices[contracted_parity.==i,:])[2]])
  end
  
  for i = 1:k
  	Aparity[contracted_sets[i]]=i
  end
  
  return Aparity, contracted_value
end