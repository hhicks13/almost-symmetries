# almost-symmetry runtime comparison

We are interested in the behavior of this algorithm particularly during the matching phase. This github repo contains a faithful implementation of
the KAS symmetry algorithm described in Kneuven and Ostrowski 2017, however it is designed to be distinctly anti-parallel to emulate a single thread of execution for the purpose of data collection and analysis. The data presented shows that indeed most of the
computational effort is carried out within the script. These results confirm our suspicions that the cost-matrix construction and the matching/refinement are expensive.
In addition to lending itself to analysis, it is built with the networkx library which can be easily interfaced with tensorflow for any related work on neural networks. Furthermore, the program automatically generates its own data at random for the purpose of simulation
and data gathering, and can easily be run from the command line. It takes in 3 arguents: number of vertices, % of edges present (sparsity of the matrix), and the budget parameter. These profiles were taken over a wide range of possible configurations of budget,vertex
amount, and edge density. 




