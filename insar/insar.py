import numpy as np
import matplotlib.patheffects as path_effects
import pandas as pd
from dateutil.parser import parse
from matplotlib import cm
import xarray as xr
from pyproj import Proj
from matplotlib.colors import DivergingNorm
import matplotlib.pyplot as plt
import os
import sys

def plotDEM(levels):

  """
  def plotDEM
  Loads and plots the digital elevation model
  """

  # Limit of map
  XLIM = (257000, 266000)
  YLIM = (2139000, 2152000)

  myProj = Proj("+proj=utm +zone=5n, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  dataset = xr.open_dataset("../dem/kilauea.nc")
  grid = np.meshgrid(dataset.variables["lon"][:], dataset.variables["lat"][:])
  x, y = myProj(*grid)
  z = np.array(dataset.variables["elev"][:])
  # Set invalid numbers to zero (sea)
  z[np.isnan(z)] = 0

  dem = (x, y, z)

  # DEM contouring
  plt.contourf(
    *dem,
    levels=levels,
    cmap=cm.Greys,
    zorder=-2
  )
 
  CS = plt.contour(
    *dem,
    levels=levels,
    colors="black",
    linewidths=1,
    zorder=0
  )
 
  # Set linewidth 
  plt.gca().grid(b=True, which="major", color="white", linewidth=0.5)
  plt.tight_layout()
  plt.xlim(*XLIM)
  plt.ylim(*YLIM)
  plt.xlabel("UTM Easting (m)")
  plt.ylabel("UTM Northing (m)")
  plt.gca().set_aspect("equal", adjustable="box")


def decomposeVectors(los1, los2, asc, desc):

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

  Ue = xr.DataArray(Ue, dims=("N", "E"), coords=los1.coords)
  Uu = xr.DataArray(Uu, dims=("N", "E"), coords=los1.coords)

  return Ue, Uu


def getExtent(dataset):

  """
  def getExtent
  Returns the extent of the dataset
  """

  E = dataset.coords["E"]
  N = dataset.coords["N"]

  return [E[0], E[-1], N[-1], N[0]]


def getGPS(start, end):

  """
  def getGPS
  Returns GPS results and vertical deformation
  """

  # Longitudes & Latitudes
  GPSMap = {
    "UWEV": [-155.291, 19.421],
    "CNPK": [-155.305, 19.392],
    "BYRL": [-155.260, 19.412],
    "DEVL": [-155.237, 19.369],
    "CRIM": [-155.274, 19.395],
    "OUTL": [-155.280, 19.387],
    "PUHI": [-155.251, 19.385],
    "MANE": [-155.271, 19.336],
    "KOSM": [-155.313, 19.360],
    "HOVL": [-155.279, 19.404],
    "NPIT": [-155.281, 19.412],
    "AHUP": [-155.262, 19.376]
  }

  # Set up UTM projection parameters
  myProj = Proj("+proj=utm +zone=5n, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  names = list()
  cs = list()
  xs = list()
  ys = list()

  # GPS values are estimated by hand
  df = pd.read_csv("GPSManual.csv", delimiter="\t")
  st = start.isoformat()[:10]
  en = end.isoformat()[:10]

  # Get correct estimate
  for name in GPSMap:

    first = df[df["DATE"] == st][name].iloc[0]
    second = df[df["DATE"] == en][name].iloc[0]

    # Difference to centimeters
    diff = 100 * (second - first)

    # Skip NaN
    if np.isnan(diff):
      continue

    # Coordinates in UTM and save everything
    xProj, yProj = myProj(*GPSMap[name])
    names.append(name)
    xs.append(xProj)
    ys.append(yProj)
    cs.append(diff)

  return xs, ys, names, cs


def plotInsar(insar):

  """
  def plotInsar
  Plots the InSAR and GPS deformation results on the map
  """

  plt.style.use("seaborn")

  # Plot the DEM below InSAR results
  plotDEM(21)

  # Get the start and end of the coverage
  start = parse(insar["start"])
  end = parse(insar["end"])

  # Load the ascending and descending traces from disk
  asc = xr.load_dataset("netcdf/asc/%s.nc" % insar["asc"])
  desc = xr.load_dataset("netcdf/desc/%s.nc" % insar["desc"])

  # Decompose deformation
  UE, UU = decomposeVectors(desc, asc, asc=insar["ascA"], desc=insar["descA"])

  # Same scale on y-axis
  vmax = np.round(np.nanmax(np.abs(UU)), 0)
  vmin = -vmax
  levels = np.linspace(vmin, vmax, 25)

  plt.contourf(asc.coords["E"], asc.coords["N"], UU, cmap=cm.seismic, norm=DivergingNorm(0), zorder=-1, vmin=vmin, levels=levels, vmax=vmax)
  bar = plt.colorbar()
  bar.set_label("Vertical Deformation (cm)")
  bar.outline.set_edgecolor("black")
  bar.outline.set_linewidth(1)

  ## Get GPS data
  x, y, s, c = getGPS(start, end)
  plt.scatter(x, y, c=c, zorder=100, edgecolor="black", cmap=cm.seismic, norm=DivergingNorm(0), vmin=vmin, vmax=vmax, linewidth=1)
  for (xi, yi, name, col) in zip(x, y, s, c):
    text = plt.annotate("%s (%d)" % (name, round(col)), (xi, yi + 100), color="white", ha="center", va="bottom")
    text.set_path_effects([path_effects.Stroke(linewidth=2, foreground="black"),
                           path_effects.Normal()])

  myProj = Proj("+proj=utm +zone=5n, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  xs = list()
  ys = list()
  cs = list()

  for (lng, lat, name) in zip(metadata["Longitude"], metadata["Latitude"], metadata["BM"]):

    x, y = myProj(lng, lat)

    area = UU.sel(
      E=slice(x-500, x+500),
      N=slice(y+500, y-500),
    )

    # Subtract the estimated deformation value at P1
    mn = np.nanmean(area.data) - insar["P1"]

    xs.append(x)
    ys.append(y)
    cs.append(np.round(mn, 1))

  # Save file
  plt.tight_layout()
  plt.title("Absolute Vertical Deformation\nBetween %s and %s\nGPS and InSAR Results" % (insar["start"], insar["end"]))
  plt.savefig("../figures/deformation/%s %s deformation.pdf" % (insar["start"], insar["end"]), bbox_inches="tight")
  plt.close()

  return cs


if __name__ == "__main__":

  """
  def __main__
  InSAR processing
  """

  # The station metadata
  metadata = pd.read_csv(
     "../locations/stations.csv",
     delimiter="\t"
  )

  lookUp = [
    {"P1": -0.9, "desc": "x20091204_x20100609", "descA": 31.2, "asc": "x20091205_x20100610", "ascA": 33.2, "start": "2009-12-04", "end": "2010-06-09"},
    {"P1": 2.7, "desc": "x20100609_x20110402", "descA": 31.2, "asc": "x20100610_x20110403", "ascA": 33.2, "start": "2010-06-09", "end": "2011-04-02"},
    {"P1": 2.4, "desc": "c20110310_c20120515", "descA": 41.5, "asc": "c20110311_c20120528", "ascA": 38.8, "start": "2011-03-10", "end": "2012-05-15"},
    {"P1": -0.5, "desc": "c20120515_c20121030", "descA": 41.5, "asc": "c20120520_c20121023", "ascA": 38.8, "start": "2012-05-15", "end": "2012-10-30"},
    {"P1": 2.3, "desc": "c20121030_c20150909", "descA": 41.5, "asc": "c20121108_c20150914", "ascA": 38.8, "start": "2012-10-30", "end": "2015-09-09"},
    {"P1": 0.4, "desc": "c20150906_c20170423", "descA": 41.5, "asc": "c20150911_c20170416", "ascA": 38.8, "start": "2015-09-06", "end": "2017-04-23"}
  ]

  # Save all estimated benchmark heights
  benchmarkHeights = [list(metadata["BM"])]

  for thing in lookUp:
    cs = plotInsar(thing)
    benchmarkHeights.append(cs)

  pd.DataFrame(np.array(benchmarkHeights).T).to_csv("./heights.csv", header=["BM", "2009-2010", "2010-2011", "2011-2012", "2012-2012", "2012-2015", "2015-2017"], sep="\t", index=False)
