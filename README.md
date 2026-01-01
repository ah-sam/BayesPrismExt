# BayesPrismExt

A modified version of BayesPrism package (https://github.com/Danko-Lab/BayesPrism.git)

## New run options

Two optional arguments were added to the main run functions (`run.prism` and `run.prism.st`) to control chain saving and Z coefficient-of-variation computation:

- `save.chain` (character): controls what (if anything) is written to an HDF5 chain file. Options:
	- `"none"` (default): no chains are saved (fastest).
	- `"theta"`: save posterior `theta` chains (cell-type fractions) only.
	- `"all"`: save both `theta` and `Z` chains (cell-type / cell-state expression matrices).

- `h5.file` (character): path to the HDF5 file used when `save.chain` is not `"none"`. If `NULL` and `save.chain != "none"` the default `gibbs_chain.h5` in the working directory is used.

- `compute.z.cv` (logical): whether to compute coefficient-of-variation for `Z` (`FALSE` by default). `theta.cv` is still always computed; computing `Z.cv` adds extra work and memory, so set `compute.z.cv=TRUE` only when you need per-gene Z CVs.

When `save.chain` is set to `"theta"` or `"all"`, chains are stored in HDF5 datasets:

- `theta`: dimensions `(n_samples, n_iterations, n_celltypes)`
- `Z` (only when `save.chain == "all"`): dimensions `(n_samples, n_iterations, n_genes, n_celltypes)`

The HDF5 implementation uses chunking and compression; adjust `h5.file` or the `init.h5.gibbs` helper in `R/run_gibbs.R` if you need different compression levels or chunk sizes.