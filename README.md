# NIPS18experiments
Code for experiments were written in C++ and Julia.

Data files are in the folder data/.

Running C++ code:

The code is in the file ExperimentsAll.cpp.

Compile with flags:

```
g++ ExperimentaAll.cpp -O3 -std=c++11
```

(this only works with C++ 11)

It runs by piping in the file containing a list of edges,
with end points space separated, followed possibly by weights
(which it sets to 1 if no weights are given).

It skips first line by default, which usually contains either
formatting information, or global meta info.
The code currently has unspecified behavior when the number
of terminals exceeds the number of vertices in the graph,
this will be fixed in a more stable version.

For example, running

```
./a.out <data/MNs-celegansneural/data.txt 
```

Outputs:

```
READING INPUT
graph with 297 vertices and 2359 edges
20 terminals: mean 2-way cut quality:1.93887,mean 4-way cut quality:2.54169,mean 8-way cut quality:2.73129,
40 terminals: mean 2-way cut quality:1.86747,mean 4-way cut quality:2.1401,mean 8-way cut quality:2.49944,
60 terminals: mean 2-way cut quality:2.68877,mean 4-way cut quality:2.01032,mean 8-way cut quality:1.85305,
80 terminals: mean 2-way cut quality:1.6418,mean 4-way cut quality:1.66332,mean 8-way cut quality:1.80395,
100 terminals: mean 2-way cut quality:1.64831,mean 4-way cut quality:1.56423,mean 8-way cut quality:1.39303,
120 terminals: mean 2-way cut quality:1.1575,mean 4-way cut quality:1.55303,mean 8-way cut quality:1.49347,
140 terminals: mean 2-way cut quality:1.17413,mean 4-way cut quality:1.40192,mean 8-way cut quality:1.33164,
160 terminals: mean 2-way cut quality:1.43052,mean 4-way cut quality:1.38893,mean 8-way cut quality:1.36536,
180 terminals: mean 2-way cut quality:1.1837,mean 4-way cut quality:1.29539,mean 8-way cut quality:1.34592,
200 terminals: mean 2-way cut quality:1.25319,mean 4-way cut quality:1.08145,mean 8-way cut quality:1.28484,

148 terminals: mean 2-way cut quality:1.15743,mean 4-way cut quality:1.36748,mean 8-way cut quality:1.64212,
74 terminals: mean 2-way cut quality:1.6568,mean 4-way cut quality:1.75227,mean 8-way cut quality:1.6558,
37 terminals: mean 2-way cut quality:1.81707,mean 4-way cut quality:2.25735,mean 8-way cut quality:2.31218,
18 terminals: mean 2-way cut quality:2.69665,mean 4-way cut quality:2.63411,mean 8-way cut quality:2.56592,
9 terminals: mean 2-way cut quality:2.33248,mean 4-way cut quality:2.86641,mean 8-way cut quality:2.85794,
```

Running Julia code:

Commands to run the code:
```
include("test.jl")
run_all_tests("<dir_name>")
```
Here, <dir_name> refers to the name of the directory containing the dataset (in a file named data.txt).

It reads all edges as a list, and sets weights to 1 if unspecified. The code also only extracts the largest connected component of the input graph to work with.

It runs the algorithm with parameters set as follows:
Number of terminals: 16
k Values: 2,4,8
sparsification factors: 2,4,8,16,32,64

On running, the output is generated in a csv file with name "julia_results_<dir_name>.csv".
