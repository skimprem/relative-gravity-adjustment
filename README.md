# Code for Review of Relative Hawaiian Gravity Campaigns (2009 - 2017)

This repository contains the code for review of relative gravity campaigns between 2009 and 2017. There is code for estimation of vertical deformation (`insar/`) and relative gravity adjustment (`hawaii.py`). The raw InSAR & gravity data files are not included in this repository.

# Code for Relative Gravity Adjustment

Supports loading of CG5, CG6 data and using WLS Inversion to complete gravity adjustment. Similar to the method described by Hwang et al., (2002), but applied to double closed-loop circuits so the system is easily overdetermined.

* Useful to find solutions for the free-air gradient
* Useful to discover gravity differences between benchmarks in a gravity network

Open `example.py` to view examples. The provided example data is from a free-air gradient measurement using both CG5 and CG6 instruments. This program can do tidal (ETERNA 3.4 or Longman, 1959) and ocean loading corrections but may require some additional dependencies (e.g., [Pygtide](https://github.com/hydrogeoscience/pygtide)) to be installed. The `ocl/` directory contains information on ocean loading corrections and relies on HARDISP inside the `hardisp/` directory for the parameter evaluation.

## Class DataLoader

Loads data from disk given a filetype ("CG5", "CG6") and a filepath:

    data = DataLoader.load(type, filepath)

Returns a DataWrapper class.

## Class DataWrapper

Wraps CG5, CG6 data and allows for the inversion. Inversion is started by passing a polynomial degree (drift) and the anchor name. If no anchor is specified, the first measurement of the circuit is taken. A tide correction ("default", "ETERNA", "Longman") can be requested and the effect of ocean loading can be removed. Setting locations for benchmarks is required when you want to do tidal corrections and can be done by:

    data.setLocations("locations/stations.csv")

The inversion result can be obtained:

    result = data.invert(degree, anchor=anchor, tide="default", loading=False)

Where the .csv is a tab delimited list of stations with latitudes and longitudes. See the example file.

Returns an `InversionResult` class.

## Class InversionResult

Wrapper for the inverted results. Two methods exist for plotting for showing the solution and the residuals respectively. To get the actual results call the `.differences` property. The method `.save` can be used to write inversion results to a file.

    result.plot(filename)
    result.plotResiduals(filename)
    result.differences
    result.save(filename)
