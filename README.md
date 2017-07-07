# sranges_validation
Independent implementation of stratigraphic range tree density of Stadler et al.

`src/srvalidate/SRTreeSimulator.java` is a class that simulates trees from the stratigraphic range fossilized birth-death process. 
The model implemented here is the budding speciation model of theorem 7, producing oriented sampled trees.

`src/srvalidate/SRTreeDensity.java` is a class that implements the computation of the probability density of a given oriented sampled tree,
as described in theorem 7.
 
Mathematical model details can be found:
 
"The fossilized birth-death model for the analysis of stratigraphic range data under different speciation concepts" (2017)
by Tanja Stadler, Alexandra Gavryushkina, Rachel C.M. Warnock, Alexei J. Drummond, Tracy A. Heath
https://arxiv.org/abs/1706.10106
