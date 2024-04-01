# Ocean Loading Parameters

Ocean loading parameters can be requested from the [free ocean loading provider](http://holt.oso.chalmers.se/loading/). The file `ocl-input.dat` illustrates an example of the free form input required for the web interface for gravity stations on Hawaii. The resulting text file should be saved as `harmonics.txt`, and contains a compilation of the parameters for all stations.

The parameter file can be split using `createcreate-harmonics.py`, and will create one file per station in the `harmonics/` directory. These files are used by the program HARDISP and calling when an ocean loading correction is requested in the Python code.
