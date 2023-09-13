# BigSurR (Basic Informatics and Gene Statistics from Unnormalized Reads R)

BigSur is a tool for single cell RNA sequencing analysis which interfaces with Seurat objects to both
- Select statistically significant features to use in clustering.
- Identify statistically significant correlations between genes.

The principles surrounding these analyses are described in ["Leveraging gene correlations in single cell transcriptomic data"][1].
[1]:"https://www.biorxiv.org/content/10.1101/2023.03.14.532643v1"

To change:
- Add documentation
- Add example usage
- Add more assertions to prevent improper use
- Change import structure to remove warnings from unused function conflicts
- Change minimum Fano factor from a set number to a quantile-based one.
