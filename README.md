# Code for Relative Gravity Adjustment

Supports loading of CG5, CG6 data and using WLS Inversion to complete gravity adjustment. Similar to the method described by Hwang et al., (2002), but applied to double closed-loop circuits so the system is easily overdetermined.

* Useful to find solutions for the free-air gradient
* Useful to discover gravity differences between benchmarks in a gravity network

Open `invert.py` to view examples. The provided example data is from a free-air gradient measurement using both CG5 and CG6 instruments.

## Class DataLoader

Loads data from disk given a filetype ("CG5", "CG6") and a filepath:

    data = DataLoader.load(type, filepath)

Returns a DataWrapper class.

## Class DataWrapper

Wraps CG5, CG6 data and allows for the inversion. Inversion is started by passing a degree (drift) between 1 - 3 and the anchor  name. If no anchor is specified, the first measurement is taken.

    result = data.invert(degree, anchor=anchor)

Returns an InversionResult class.

## Class InversionResult

Wrapper for the inverted results. Two methods exist for plotting for showing the solution and the residuals respectively.

    result.plot()
    result.plotResiduals()

