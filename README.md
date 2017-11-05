# phenoClust
This package implements the PhenoGraph algorithm, which is based on the Louvain algorithm.

## Installation
To install directly from github requires `devtools`.
```R
devtools::install_github('yishaishimoni/phenoClust')
```

This will install a package called `phenoGraph` on your system.

## Usage
```R
require(phenoGraph)
require(ggplot2)
iris_unique <- unique(iris) # Remove duplicates
data <- as.matrix(iris_unique[,1:4])
Rphenograph_out <- phenoClust(t(data))
modularity(G=Rphenograph_out$G, C = Rphenograph_out$C)
iris_unique$phenograph_cluster <- factor(Rphenograph_out$C)
ggplot(iris_unique, aes(x=Sepal.Length, y=Sepal.Width, col=Species, shape=phenograph_cluster)) + geom_point(size = 3)+theme_bw()
```

## References
The example is adapted from another independent implementation of phenoGraph at `https://github.com/JinmiaoChenLab/Rphenograph`

The algorithm is taken from 
Levine JH, Simonds EF, Bendall SC, Davis KL, Amir ED, Tadmor MD, et al. Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis. Cell. Elsevier Inc.; 2015; 1–14. doi:10.1016/j.cell.2015.05.047
