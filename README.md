# Source code from the paper "Mapping robust multiscale communities in chromosome contact networks"

_Anton Holmgren, Dolores Bernenko, Ludvig Lizana_

[https://arxiv.org/abs/2212.08456](https://arxiv.org/abs/2212.08456)

In this paper, we investigate how community detection methods can find groups of DNA regions in the human chromosome where internal interactions are higher than expected.
To do this, we use the community detection method [GenLouvain for Matlab](https://github.com/GenLouvain/GenLouvain), adapted with the fractal globule null model from polymer physics.

## Dependencies

- [GenLouvain](https://github.com/GenLouvain/GenLouvain)
- [Solution landscape clustering](https://github.com/mapequation/partition-validation)
- [Significance clustering](https://github.com/mapequation/multilevel-significance-clustering)
- C++ compiler
- Rust
- Python

# How to use

The idea for this paper is to detect scales for which the community solution landscape is robust.
As community detection is often implemented using a stochastic search algorithm sensitive to starting conditions, one often performs multiple searches to find the best (highest quality) solution.
Then, one studies the obtained partitions similarity and quality.
If all partitions are similar, the solution landscape is peaked.
If partitions are dissimilar but the quality is similar, then the solution landscape is degenerate.
Often, for HiC data, there is a middle ground where there are multiple peaks of varying quality (local maxima).

To study this phenomenon, we use a solution landscape clustering method to group similar solutions together (solution landscape clustering code).
Next, we study the variability in each "peak" to see if some DNA regions are more easy to robustly assign to communities than others (significance clustering code).

