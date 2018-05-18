using DataFrames
include("ContractGraph.jl")

function testApproxRatio(path, k_factor_array, sparsification_factor_array, no_of_terminals)
  file = string(path, "data.txt")
  if path == "random"
    println("Using complete graph with random integer weights")
    A = convert(SparseMatrixCSC{Float64,Int64}, sparse(rand(1:100,100,100)))
    A = A + A'
    for i in 1:100
      A[i,i]=0
    end
	dropzeros(A)
  else
    println("Input file ", file)
    A = extractGraphFromFile(file)
  end
  
  n = size(A,1)
  m = nnz(A)
  
  println("Number of vertices: ", n)
  println("Number of edges: ", m)
  
  df = DataFrame(Input = String[], Vertices = Int64[], Edges = Int64[], kFactor = Int64[], SparsificationFactor = Int64[], NumberOfTerminals = Int64[], ActualCutValue = Float64[], ApproxCutValue_MEAN = Float64[], CutQuality_MEAN = Float64[])

  for k in k_factor_array
	println("k factor: ", k)
    
	if no_of_terminals == -1
		println("Reading terminals from file")
		
	else
		terminals=sample(1:n, no_of_terminals, replace = false)
			
		terminal_sets = Array{Int64,1}[]
		x = sample(1:no_of_terminals, k-1,replace = false, ordered = true)

		x=append!([0],x)
		x=append!(x,[no_of_terminals])

		for i = 1:k
			append!(terminal_sets,[terminals[x[i]+1:x[i+1]]])
		end
	end
	
	println("Computing actual minST-k-waycut")
	actual_value = minSTkWayCut(A, terminal_sets, k)[2]

	for sparsification_factor in sparsification_factor_array[sparsification_factor_array.<= div(n,no_of_terminals)]
		println("Sparsification factor: ", sparsification_factor)

		sparsifier_size = div(n,sparsification_factor)
	
	
		println("Computing minST-",k,"-way-cut using contractions")
		
		approx_value_array = zeros(100)
		for i in 1:100
		  approx_value_array[i] = contractMinSTkWayCut(A,terminal_sets, k, sparsifier_size)[2]
		end
		
		approx_value_mean = mean(approx_value_array)
		
		println("Average-case quality of cut: ", approx_value_mean/actual_value)
		
		dfrow = (path, n,m, k, sparsification_factor, no_of_terminals, actual_value, approx_value_mean, approx_value_mean/actual_value)
		
		push!(df,dfrow)        		
	end
  end

  return df
end

"""
df = DataFrame(Input = String[], Vertices = Int64[], Edges = Int64[], kFactor = Int64[], SparsificationFactor = Int64[], NumberOfTerminals = Int64[], ActualCutValue = Float64[], ApproxCutValue_MEAN = Float64[], CutQuality_MEAN = Float64[])
  
for no_of_terminals in (2,4,8,16,32)
  for sparsification_factor in (2:2:10)
    for i in (1:1)
	  tup2 = testApproxRatio("random",sparsification_factor,no_of_terminals)
	  dir = "CompleteGraph"
	  dfrow = (dir,tup2...)
	  push!(df,dfrow)
	end
  end
end

CSV.write("results_randomgraph.csv",df)
"""

#df2 = DataFrame(Input = String[], Vertices = Int64[], Edges = Int64[], kFactor = Int64[], SparsificationFactor = Int64[], NumberOfTerminals = Int64[], ActualCutValue = Float64[], ApproxCutValue_MEAN = Float64[], CutQuality_MEAN = Float64[])

function run_all_tests(dataset::String)
	sparsification_factor_array = [2,4,8,16,32,64]
	k_factor_array = [2,4,8]
	no_of_terminals = 16
	"""
	if terminal_input == "given"
		println("Terminals provided in input")
		no_of_terminals_ = -1
	else
		println("Using 16 random terminals")
		no_of_terminals = 16
	end
	"""
	df = testApproxRatio(string("data/",dataset,"/"), k_factor_array, sparsification_factor_array, no_of_terminals)
	CSV.write(string("JuliaResults/results_",dataset,".csv"),df)  
end

