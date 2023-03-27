import numpy as np
import isosplit6_cpp


def isosplit6(X):
	X = X.astype(np.float64, copy=False, order='C') #copies only if type changes, but note that we require C order, which is essential for the interface
	N = X.shape[0]
	M = X.shape[1]
	labels_out=np.zeros([N]).astype(np.int32)
	isosplit6_cpp.isosplit6_fn(labels_out, X)
	return labels_out

def isocut6(X):
	X = X.astype(np.float64, copy=False, order='C') #copies only if type changes, but note that we require C order, which is essential for the interface
	N = X.shape[0]
	output=np.zeros((2,), dtype=np.float64)
	isosplit6_cpp.isocut6_fn(output, X)
	dipscore = output[0]
	cutpoint = output[1]
	return dipscore, cutpoint