# bioCausal
A pipeline for causal analysis on biology data.

bioCausal provides an easy-to-use pipeline for feature selection, latent variable analysis and causal analysis from biological data

## Install Dependencies
```{r}
dependencies = c("ggnetwork", "igraph", "intergraph", "ggpubr", "ggplot2", "glmnet", "limma", "qvalue", "RColorBrewer", "pheatmap"")
BiocManager::install(dependencies)

devtools::install_github("wgmao/PLIER")
```

## Installation
```{r}
devtools::install_github("zhu-yh1/bioCausal/bioCausal")
```

## Example
### Load package
```{r}
library(bioCausal)
```
### Load data (and) metadata
```{r}
data = read.table($PATH_TO_DATA$)
# scale data
data = tscale(data)

metadata = read.table($PATH_TO_METADATA$)
```
### Select focal variable
```{r}
# FV could be a molecular feature
FV = data[1,]

# FV could also be a metadata feature
FV = metadata[1,]
```
### Fearure selection
1. Feature selection using FV
```{r}
# provide contrast if FV is binomial
selectionRes = featureSelection(data = data, focalVariable = FV, class = "binomial", limmaCutoff = 30, glmnetCutoff = 20, contr="A-B")
selected = selectionRes$selected
```
2. Plot selection result
```{r}
# yvalue can be "pval" or "adj.pval".
plotSelection(rownames(data), selectionRes, sign = T, yvalue="adj.pval")
```
### Latent variable analysis (optional)
```{r}
decompRes = simpleDecomp(tscale(data), k = 20, adaptive.frac = 0.01)
```
### Causal analysis
1. Get input data
```{r}
# combine molecular feature, LV, and metadata as input
input = rbind(data[selected,], tscale(decompRes$B[lv.use,], metadata)
```

2. run causal analysis using notears
```
# notears
ntres = notearsInterceptMultiLoss(t(input), lambda1=0.01)

# notears with contraints
# no_edge is a feature*feature 0,1 matrix with 1 specifying no edge constraint
ntres.constraint = notearsInterceptMultiLoss(t(input), lambda1=0.01, no_edge = no_edge)
```
3. visualization
```{r}
# set gvarType and gvarShape
gvarType = feature_grp1
gvarShape = feature_grp2
net = network_visualize(ntres.constraint$graph, gvarType = gvarShape, gvarShape = gvarType)
net$p
```
