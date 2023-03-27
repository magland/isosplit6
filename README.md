# isosplit6

Isosplit: A Non-Parametric Method for Unimodal Clustering

Isosplit is a non-parametric clustering method that does not require adjustable parameters nor parametric assumptions about the underlying cluster distributions. The only assumption is that clusters are unimodal and separated from one another by hyperplanes of relatively low density. The technique uses a variant of Hartigan's dip statistic and isotonic regression in its kernel operation.

Motivation: Many clustering algorithms require the tuning of parameters for each application or dataset, making them unsuitable for automated procedures that involve clustering. Some techniques require an initial estimate of the number of clusters, while density-based techniques typically require a scale parameter. Other parametric methods, such as mixture modeling, make assumptions about the underlying cluster distributions.

# Installation and usage

```bash
pip install isosplit6
```

```python
from isosplit6 import isosplit6

# Prepare a N x M Numpy Array
# N = number of features
# M = number of dimensions
features = ...

cluster_labels = isosplit6(features)
```
