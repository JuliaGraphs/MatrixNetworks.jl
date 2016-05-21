Diffusions design document
==========================

Should work for 

* PageRank
* Heat Kernel
* Katz
* Single-source PageRank (ACL)
* Single-source Heat kernel (hkrelax)
* Single-source general (diffrelax)
* And be as efficient and general as possible.

Idea:

High level-routines

aclpagerank(A, alpha, Int, tol) 
aclpagerank(A, alpha, Set, tol)
aclpagerank(A, alpha, Dict, tol) # weighted set
aclpagerank(A, alpha, DegreeWeightedSet(Set), tol)
# these all check the symmetry/undirected nature
aclpagerank(SymmetryChecked(A), alpha, DegreeWeightedSet(Set), tol)
# this skips that...

pagerank(A, alpha) # tol is machine precision
pagerank(A, alpha, tol) # is a variable

personalized_pagerank(A,alpha,Set) # personalized on a set (uniformly)
personalized_pagerank(A,alpha,Int) # personalized on a node
personalized_pagerank(A,alpha,sparsevec)
personalized_pagerank(A,alpha,Dict)
personalized_pagerank(A,alpha,DegreeWeightedSet(Set))  

Low level-routines

pagerank_power!
aclpagerank_push!

Eventually, 

we'd like some type of intermediate interface:

P = ACLPageRankProblem(A) # does symmetry checking 
P.solve(alpha, v, ...) 

PageRankDiffusion
- This is a parametric type.

    


