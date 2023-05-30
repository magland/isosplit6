# isosplit6

Isosplit is a non-parametric clustering method that does not require adjustable parameters nor parametric assumptions about the underlying cluster distributions. The only assumption is that clusters are unimodal and separated from one another by hyperplanes of relatively low density. The technique uses a variant of Hartigan's dip statistic and isotonic regression in its kernel operation.

Motivation: Many clustering algorithms require the tuning of parameters for each application or dataset, making them unsuitable for automated procedures that involve clustering. Some techniques require an initial estimate of the number of clusters, while density-based techniques typically require a scale parameter. Other parametric methods, such as mixture modeling, make assumptions about the underlying cluster distributions.

Isosplit is used by the [MountainSort](https://github.com/magland/mountainsort5) spike sorting algorithm.

[preprint](https://arxiv.org/abs/1508.04841)

# Installation and usage

```bash
pip install isosplit6
```

```python
from isosplit6 import isosplit6

# Prepare a N x M Numpy Array
# N = number of observations
# M = number of features
features = ...

cluster_labels = isosplit6(features)
```
