import numpy as np
import matplotlib.patheffects as path_effects
import pandas as pd
from dateutil.parser import parse
from matplotlib import cm
import xarray as xr
from matplotlib.legend_handler import HandlerTuple
from pyproj import Proj
from matplotlib.colors import DivergingNorm
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os
import sys

class InsarWrapper():

  """
  class InsarWrapper
  Wraps results from InSAR vert. deformation analysis
  """

  # Projection parameters
  myProj = Proj("+proj=utm +zone=5n, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  # GPS Stations and locations
  GPSMap = {
    "UWEV": (-155.291, 19.421),
    "CNPK": (-155.305, 19.392),
    "BYRL": (-155.260, 19.412),
    "DEVL": (-155.237, 19.369),
    "CRIM": (-155.274, 19.395),
    "OUTL": (-155.280, 19.387),
    "PUHI": (-155.251, 19.385),
    "MANE": (-155.271, 19.336),
    "KOSM": (-155.313, 19.360),
    "HOVL": (-155.279, 19.404),
    "NPIT": (-155.281, 19.412),
    "AHUP": (-155.262, 19.376)
  }


  def __init__(self, trace):

    """
    def InsarWrapper.__init__
    Initializes the class and saves traces
    """

    self.trace = trace
    self.deformationField = self.calculateDeformation(trace)


  def plotDEM(self, levels):
  
    """
    def InsarWrapper.plotDEM
    Loads and plots the digital elevation model
    """
  
    # Limit of map
    XLIM = (257000, 266000)
    YLIM = (2139000, 2152000)
  
    # Load the DEM
    dataset = xr.open_dataset("../dem/kilauea.nc")
    grid = np.meshgrid(dataset.variables["lon"][:], dataset.variables["lat"][:])
    x, y = self.myProj(*grid)
    z = np.array(dataset.variables["elev"][:])
    # Set invalid numbers to zero (sea)
    z[np.isnan(z)] = 0
  
    # DEM contouring
    plt.contourf(
      x, y, z,
      levels=levels,
      cmap=cm.Greys,
      zorder=-2
    )
   
    CS = plt.contour(
      x, y, z,
      levels=levels,
      colors="black",
      linewidths=1,
      zorder=0
    )
   
    # Plot styling 
    plt.xlim(*XLIM)
    plt.ylim(*YLIM)
    plt.xlabel("UTM Easting (m)")
    plt.ylabel("UTM Northing (m)")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.gca().grid(b=True, which="major", color="white", linewidth=0.5)


  def decomposeVectors(self, los1, los2, asc, desc):
  
    """
    def decomposeVectors
    Decomposes the line of sight vectors from the ascending and descending
    traces to vertical / horizontal deformation
    """
  
    # Angles to radians
    asc = np.radians(asc)
    desc = np.radians(desc)
    clos1 = np.array(los1.__xarray_dataarray_variable__)
    clos2 = np.array(los2.__xarray_dataarray_variable__)
  
    # We should solve simple linear equations:
    # LOS_DESC = Ue * sin(λ1) + Uu * cos(λ1)
    # LOS_ASC = -Ue * sin(λ2) + Uu * cos(λ2)
  
    # In vector form:
    #
    # [LOS_DSC] = [ sin(λ1) cos(λ1)] [Ue]
    # [LOS_ASC] = [-sin(λ2) cos(λ2)] [Uu]
    #
    # Create the inverse matrix to solve for Ue and Uu
    inverse = np.linalg.inv(
      np.array([
        [np.sin(desc), np.cos(desc)],
        [-np.sin(asc), np.cos(asc)],
      ])
    )
  
    # Multiply out by hand
    Ue = inverse[0][0] * clos1 + inverse[0][1] * clos2
    Uu = inverse[1][0] * clos1 + inverse[1][1] * clos2
  
    # To xArray
    Ue = xr.DataArray(Ue, dims=("N", "E"), coords=los1.coords)
    Uu = xr.DataArray(Uu, dims=("N", "E"), coords=los1.coords)
  
    return Ue, Uu
  

  def getGPS(self, start, end):
  
    """
    def getGPS
    Returns GPS results and vertical deformation
    """
  
    # Set up UTM projection parameters
    names = list()
    cs = list()
    xs = list()
    ys = list()
  
    # The GPS values are estimated by hand because of rapid deformation
    # moving averages do not work and daily estimates are too scattered
    df = pd.read_csv("GPSManual.csv", delimiter="\t")
    st = start.isoformat()[:10]
    en = end.isoformat()[:10]
  
    measured = dict()
  
    for name in self.GPSMap:
  
      # Get the proper rows
      first = df[df["DATE"] == st][name].iloc[0]
      second = df[df["DATE"] == en][name].iloc[0]
  
      # Difference between the dates in centimeters
      diff = 100 * (second - first)
  
      # Save
      measured[name] = diff
  
      # Skip NaN
      if np.isnan(diff):
        continue
  
      # Coordinates in UTM and save everything
      xProj, yProj = self.myProj(*self.GPSMap[name])
      names.append(name)
      xs.append(xProj)
      ys.append(yProj)
      cs.append(diff)
  
    return xs, ys, names, cs, measured


  def estimateDeformation(self, deformationField, lng, lat):
  
    """
    def InsarWrapper.estimateDeformation
    Estimated deformation at a longitude, latitude for a given deformation field
    """
  
    x, y = self.myProj(lng, lat)
  
    # Use the mean of a 1x1 km area
    area = deformationField.sel(
      E=slice(x - 250, x + 250),
      N=slice(y + 250, y - 250),
    )
  
    mean = np.nanmean(area)

    # If not found try a bigger area: not ideal but better to get a height
    if np.isnan(mean):
      area = deformationField.sel(
        E=slice(x - 1000, x + 1000),
        N=slice(y + 1000, y - 1000),
      )
      mean = np.nanmean(area)

    return mean


  def predictGPS(self, deformationField):
  
    """
    def InsarWrapper.predictGPS
    Predicts deformation values at position of GPS receivers using insar field
    """

    predicted = dict()
  
    for station in self.GPSMap:
      predicted[station] = self.estimateDeformation(deformationField, *self.GPSMap[station])
  
    return predicted


  def anchorInsar(self, predicted, measured):

    """
    def InsarWrapper.anchorInsar
    Anchors the predicted GPS values to the measured GPS values and records the required offset
    """

    def insarObjective(x):

      """
      def InsarWrapper.anchorInsar.insarObjective
      Minimizes the squared residuals between InSAR and GPS measurements
      """

      residuals = list() 

      for name in self.GPSMap:

        # Must have a value for the InSAR predicted and GPS measured points
        if np.isnan(predicted[name]) or np.isnan(measured[name]):
          continue

        # Save residuals and minimize the squared sum
        residuals.append((predicted[name] + x) - measured[name])

      return np.sum(np.square(residuals))


    # Scipy minimize
    solution = minimize(
      insarObjective,
      [0],
      method="Nelder-Mead",
      options={"maxiter": 10000}
    )

    # This is the required offset of the InSAR interferogram
    return solution.x[0]


  def calculateDeformation(self, trace):

    """
    def InsarWrapper.calculateDeformation
    Calculates the vertical deformation field based on LOS data
    """

    # Get the start and end of the coverage
    start = parse(trace["start"])
    end = parse(trace["end"])
 
    # Load the ascending and descending traces from disk
    asc = xr.load_dataset("netcdf/asc/%s.nc" % trace["asc"])
    desc = xr.load_dataset("netcdf/desc/%s.nc" % trace["desc"])
 
    # Decompose LOS deformation vectors
    # UU is vertical
    UE, UU = self.decomposeVectors(desc, asc, asc=trace["ascA"], desc=trace["descA"])

    # Predict the value at the GPS receiver locations from InSAR data
    prediction = self.predictGPS(UU)

    # Get the observed GPS data
    x, y, s, c, measured = self.getGPS(start, end)

    # Calculate optimal offset between InSAR and GPS measurements
    offset = self.anchorInsar(prediction, measured)

    # Add the optimal offset to the InSAR field data. Now it is "least squares" anchored to GPS observations
    # Return coordinates too
    return asc.coords["E"], asc.coords["N"], UU + offset


  def plot(self, metadata):
  
    """
    def plotInsar
    Plots the InSAR and GPS deformation results on the map
    """
  
    # Styling
    plt.style.use("seaborn")
  
    # Plot the DEM below InSAR results
    self.plotDEM(21)
  
    # This is the deformation field / coordinates
    E, N, UU = self.deformationField

    # Same scale on y-axis
    vmax = np.round(np.nanmax(np.abs(UU)), 0)
    vmin = -vmax
    levels = np.linspace(vmin, vmax, 25)
  
    # Plot deformation field
    plt.contourf(E, N, UU, cmap=cm.seismic, norm=DivergingNorm(0), zorder=-1, vmin=vmin, levels=levels, vmax=vmax)
    bar = plt.colorbar()
    bar.set_label("Vertical Deformation (cm)")
    bar.outline.set_edgecolor("black")
    bar.outline.set_linewidth(1)
  
    # Get the start and end of the coverage
    start = parse(self.trace["start"])
    end = parse(self.trace["end"])

    x, y, s, c, measured = self.getGPS(start, end)

    # Plot the measured GPS data
    plt.scatter(x, y, c=c, zorder=100, edgecolor="black", cmap=cm.seismic, norm=DivergingNorm(0), vmin=vmin, vmax=vmax, linewidth=1, label="GPS Receivers", s=40)
    for (xi, yi, name, col) in zip(x, y, s, c):
      text = plt.annotate("%s (%d)" % (name, round(col)), (xi, yi + 100), color="white", ha="center", va="bottom", fontsize=8)
      text.set_path_effects([path_effects.Stroke(linewidth=2, foreground="black"),
                             path_effects.Normal()])
  

    xben = list()
    yben = list()
    cben = list()

    for (name, lng, lat) in zip(metadata["BM"], metadata["Longitude"], metadata["Latitude"]):
      c = self.estimateDeformation(self.deformationField[2], lng, lat)
      xx, yy = self.myProj(lng, lat)
      xben.append(xx)
      yben.append(yy)
      cben.append(c)

    plt.scatter(xben, yben, c=cben, edgecolor="black", linewidth=1, cmap=cm.seismic, norm=DivergingNorm(0), vmin=vmin, vmax=vmax, s=20)

    cmap = cm.get_cmap("seismic", 21)
    a = plt.scatter(np.nan, np.nan, color=cmap(0.2), edgecolor="black", linewidth=1, s=20)
    b = plt.scatter(np.nan, np.nan, color=cmap(0.4), edgecolor="black", linewidth=1, s=20)
    c = plt.scatter(np.nan, np.nan, color=cmap(0.5), edgecolor="black", linewidth=1, s=20)
    d = plt.scatter(np.nan, np.nan, color=cmap(0.6), edgecolor="black", linewidth=1, s=20)
    e = plt.scatter(np.nan, np.nan, color=cmap(0.8), edgecolor="black", linewidth=1, s=20)

    a2 = plt.scatter(np.nan, np.nan, color=cmap(0.2), edgecolor="black", linewidth=1, s=40)
    b2 = plt.scatter(np.nan, np.nan, color=cmap(0.4), edgecolor="black", linewidth=1, s=40)
    c2 = plt.scatter(np.nan, np.nan, color=cmap(0.5), edgecolor="black", linewidth=1, s=40)
    d2 = plt.scatter(np.nan, np.nan, color=cmap(0.6), edgecolor="black", linewidth=1, s=40)
    e2 = plt.scatter(np.nan, np.nan, color=cmap(0.8), edgecolor="black", linewidth=1, s=40)

    plt.legend([(a2, b2, c2, d2, e2), (a, b, c, d, e)], ("GPS Receivers", "Gravity Benchmarks"), scatterpoints=1, numpoints=2, handler_map={tuple: HandlerTuple(ndivide=None)}, frameon=True, prop={"size": 8})


    start = start.date()
    end = end.date()
    # Save file
    plt.tight_layout()
    plt.title("Absolute Vertical Deformation\nBetween %s and %s\nGPS and InSAR Results" % (start, end))
    plt.savefig("../figures/deformation/%s %s deformation.pdf" % (start, end), bbox_inches="tight")
    plt.close()


  def getBenchmarkPredictions(self, metadata):

    """
    def InsarWrapper.getBenchmarkPredictions
    Gets predictions for the benchmark
    """

    dhs = list()
  
    # Calculate the deformation at P1 (anchor)
    P1Def = self.estimateDeformation(self.deformationField[2], -155.301304, 19.437016)

    # Estimate the vertical deformation at each benchmark
    for (name, lng, lat) in zip(metadata["BM"], metadata["Longitude"], metadata["Latitude"]):
      dhs.append(self.estimateDeformation(self.deformationField[2], lng, lat))

    return np.round(dhs - P1Def, 1)



if __name__ == "__main__":

  """
  def __main__
  InSAR / GPS processing for height estimation
  """

  traces = [
    {"P1": -0.9, "desc": "x20091204_x20100609", "descA": 31.2, "asc": "x20091205_x20100610", "ascA": 33.2, "start": "2009-12-04", "end": "2010-06-09"},
    {"P1": 2.7, "desc": "x20100609_x20110402", "descA": 31.2, "asc": "x20100610_x20110403", "ascA": 33.2, "start": "2010-06-09", "end": "2011-04-02"},
    {"P1": 2.4, "desc": "c20110310_c20120515", "descA": 41.5, "asc": "c20110311_c20120528", "ascA": 38.8, "start": "2011-03-10", "end": "2012-05-15"},
    {"P1": -0.5, "desc": "c20120515_c20121030", "descA": 41.5, "asc": "c20120520_c20121023", "ascA": 38.8, "start": "2012-05-15", "end": "2012-10-30"},
    {"P1": 2.3, "desc": "c20121030_c20150909", "descA": 41.5, "asc": "c20121108_c20150914", "ascA": 38.8, "start": "2012-10-30", "end": "2015-09-09"},
    {"P1": 0.4, "desc": "c20150906_c20170423", "descA": 41.5, "asc": "c20150911_c20170416", "ascA": 38.8, "start": "2015-09-06", "end": "2017-04-23"}
  ]

  # The station metadata
  metadata = pd.read_csv(
     "../locations/stations.csv",
     delimiter="\t"
  )

  # Save all estimated benchmark heights
  benchmarkHeights = [list(metadata["BM"])]

  for thing in traces:
    wrapper = InsarWrapper(thing)
    wrapper.plot(metadata)
    benchmarkHeights.append(wrapper.getBenchmarkPredictions(metadata))

  pd.DataFrame(np.array(benchmarkHeights).T).to_csv("./heights.csv", header=["BM", "2009-2010", "2010-2011", "2011-2012", "2012-2012", "2012-2015", "2015-2017"], sep="\t", index=False)
