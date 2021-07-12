import os
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.ticker as plticker
import matplotlib.pyplot as plt

from matplotlib.legend_handler import HandlerTuple
from pyproj import Proj
from matplotlib import cm


def loadResult(campaign, instrument, filename):

  """
  def loadResult
  Loads results from the results directory
  """

  filepath = os.path.join("results", campaign, instrument, filename)

  return pd.read_csv(
    filepath,
    delimiter="\t",
    comment="#"
  )


def plotDEM(myProj):

  """
  def plotDEM
  Plots a Kilaueau DEM 
  """

  levels = 21
  dataset = xr.open_dataset("dem/kilauea.nc")

  # Must meshgrid for contourf
  x, y = myProj(*np.meshgrid(dataset.variables["lon"][:], dataset.variables["lat"][:]))
  z = np.array(dataset.variables["elev"][:])
  # Set invalid numbers to zero (sea)
  z[np.isnan(z)] = 0

  # DEM contouring
  plt.contourf(
    x, y, z,
    levels=levels,
    cmap=cm.Greys,
    zorder=0
  )

  plt.contour(
    x, y, z,
    levels=levels,
    colors="black",
    linewidths=1,
    zorder=0
  )

  # Set linewidth 
  plt.gca().grid(b=True, which="major", color="white", linewidth=0.75)


def compareMap(collector, instrument, one, two, scale):

  """
  def compareMap
  Plots differences between two campaigns on a map
  """

  # Limit of map
  XLIM = (257000, 266000)
  YLIM = (2139000, 2152000)

  plt.style.use("seaborn")

  # Create UTM projection
  myProj = Proj("+proj=utm +zone=5n, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  plotDEM(myProj)

  filename = "%s %s %s" % (one, two, instrument)

  plt.title("Change in Gravity Difference with Anchor P1 \n Between %s and %s \n Instrument %s" % (one, two, instrument))

  dataOne = collector[one][instrument]
  dataTwo = collector[two][instrument]

  locations = pd.read_csv("locations/stations.csv", delimiter="\t")
  xProj, yProj = myProj(locations["Longitude"], locations["Latitude"])

  benchmarks = list()
  missing = list()
  colors = list()

  for station, longitude, latitude in zip(locations["BM"], xProj, yProj):

    if station == "P1":
      b = plt.scatter(longitude, latitude, marker="*", s=200, color="yellow", linewidth=1, edgecolor="black")
      continue
    elif not station in dataOne or not station in dataTwo:
      missing.append((longitude, latitude))
      continue

    # Calculate the difference
    dg, std = difference(dataOne, dataTwo, station)

    # Read height
    dg += getHeight(one, two, station)

    benchmarks.append((longitude, latitude))
    colors.append(dg)

  # Create map
  cmap = cm.get_cmap("seismic", 21)
  a = plt.scatter(*np.array(missing).T, color="white", s=15, cmap=cmap, edgecolor="black", linewidth=1)
  plt.scatter(*np.array(benchmarks).T, c=colors, cmap=cmap, edgecolor="black", linewidth=1)

  # Colorbar
  bar = plt.colorbar()
  bar.mappable.set_clim(-scale, scale)
  bar.outline.set_edgecolor("black")
  bar.outline.set_linewidth(1)

  # Neat layout of the legend! This is unimportant and just visually pleasing
  c = plt.scatter(np.nan, np.nan, color=cmap(0.2), edgecolor="black", linewidth=1)
  d = plt.scatter(np.nan, np.nan, color=cmap(0.4), edgecolor="black", linewidth=1)
  e = plt.scatter(np.nan, np.nan, color=cmap(0.5), edgecolor="black", linewidth=1)
  f = plt.scatter(np.nan, np.nan, color=cmap(0.6), edgecolor="black", linewidth=1)
  g = plt.scatter(np.nan, np.nan, color=cmap(0.8), edgecolor="black", linewidth=1)
  plt.legend([a, b, (c, d, e, f, g)], ("Unavailable", "Anchor (P1)", "Campaign"), scatterpoints=1, numpoints=2, handler_map={tuple: HandlerTuple(ndivide=None)}, frameon=True, prop={"size": 8})

  # Figure setup
  plt.tight_layout()
  plt.xlim(*XLIM)
  plt.ylim(*YLIM)
  plt.xlabel("UTM Easting (m)")
  plt.ylabel("UTM Northing (m)")
  plt.gca().set_aspect("equal", adjustable="box")

  plt.savefig("figures/map/%s.pdf" % filename, bbox_inches="tight")
  plt.close()


def difference(one, two, station):

  """
  def difference
  Reads heights from a .CSV file
  """

  # Difference between two campaigns
  dg = two[station][0]  - one[station][0]

  # 2 times standard deviation from differences sqrt(σ1^2 + σ2^2)
  std = 2 * np.sqrt(two[station][1] ** 2 + one[station][1] ** 2)

  return dg, std


def getHeight(one, two, station):

  """
  def getHeight
  Reads heights from a .CSV file
  """

  df = pd.read_csv(
    "heights.csv",
    delimiter="\t",
    comment="#"
  )

  # Undo the effect of height change (negative height will give positive gravity)
  return 3.086 * df[df["Benchmark"] == station].iloc[0]["%s %s" % (one, two)]


def compareDistribution(collector, instrument, one, two, scale):

  """
  def compareDistribution
  Plots distribution of data against the benchmarks
  """

  plt.style.use("seaborn")

  plt.figure(figsize=(5, 8))
  plt.title("Change in Gravity Difference with Anchor P1 \n Between %s and %s \n Instrument %s" % (one, two, instrument))

  # Output filename
  filename = "%s %s %s" % (one, two, instrument)

  # Fetch the results
  dataOne = collector[one][instrument]
  dataTwo = collector[two][instrument]

  locations = pd.read_csv("locations/stations.csv", delimiter="\t")

  for station in sorted(locations["BM"]):

    if station == "P1":
      continue

    if not station in dataOne or not station in dataTwo:
      plt.barh(station, 0, edgecolor="black", xerr=0, capsize=0, linewidth=0, color=color, error_kw=dict(capthick=0))
      continue

    dg, std = difference(dataOne, dataTwo, station)

    # Height correction
    dg += getHeight(one, two, station)

    if (np.abs(dg) - std) < 0 or dg == 0:
      color = "white"
    elif dg > 0:
      color = "red"
    elif dg < 0:
      color = "blue"

    #plt.bar(station, dg, edgecolor="black", yerr=std, capsize=2, linewidth=1, color=color, error_kw=dict(capthick=1))
    plt.barh(station, dg, edgecolor="black", xerr=std, capsize=2, linewidth=1, color=color, error_kw=dict(capthick=1))

  plt.xlim(-scale, scale)
  plt.margins(y=0)
  plt.ylabel("Benchmark")
  plt.xlabel("Change in Gravity Difference (µGal)")
  plt.yticks(ha="left", position=(-0.13, 0))
  plt.tight_layout()

  plt.savefig("figures/dist/%s.pdf" % filename, bbox_inches="tight")
  plt.close()


def relativeTo(one, onesd, two, twosd):

  return one + two, np.sqrt(onesd ** 2 + twosd ** 2)

if __name__ == "__main__":

  """
  def __main__
  Main entry point to calculate differences between two campaigns and plot the results
  to a geographical map or distribution plot
  """

  collector = {}

  for campaign in os.listdir("results"):

    collector[campaign] = {"578": {}, "579": {}}

    for instrument in os.listdir(os.path.join("results", campaign)):

      for filename in os.listdir(os.path.join("results", campaign, instrument)):

        if filename.startswith("."):
          continue

        # Read optimization results
        df = loadResult(campaign, instrument, filename)

        # Go over each row
        for _, row in df.iterrows(): 

          # Extract the benchmark, dg, and standard deviation
          benchmark = row["Benchmark"]
          gravity = row["Gravity (µGal)"]
          sd = row["SD (µGal)"]
          anchor = row["Anchor"]

          # Sometimes a circuit is measured relative to HVO1
          # convert to P1 THROUGH its difference with HVO41: this is not ideal..
          if anchor == "HVO41":

            if filename == "578_2012-06-22.dat":
              gravity, sd = relativeTo(gravity, sd, 29891.0, 5.67)
            if filename == "578_2011-03-23.dat":
              gravity, sd = relativeTo(gravity, sd, 29799.0, 2.61)
            if filename == "578_2012-11-27.dat":
              gravity, sd = relativeTo(gravity, sd, 29947.0, 2.76)

            if filename == "579_2012-06-22.dat":
              gravity, sd = relativeTo(gravity, sd, 29915.0, 9.1)
            if filename == "579_2011-03-23.dat":
              gravity, sd = relativeTo(gravity, sd, 29811.0, 5.9)
            if filename == "579_2012-11-27.dat":
              gravity, sd = relativeTo(gravity, sd, 29965.0, 2.87)

          # Stations that are measured twice in a campaign: let's pick the best ones and skip the following:
          if (filename == "578_2010-07-15.dat" and benchmark == "HVO41") or\
             (filename == "579_2010-07-15.dat" and benchmark == "HVO41") or\
             (filename == "578_2010-06-30.dat" and benchmark == "21YY") or\
             (filename == "579_2010-06-30.dat" and benchmark == "21YY") or\
             (filename == "578_2010-07-09.dat" and benchmark == "HVO33") or\
             (filename == "579_2010-07-09.dat" and benchmark == "HVO33") or\
             (filename == "578_2010-07-15.dat" and benchmark == "204YY") or\
             (filename == "579_2010-07-15.dat" and benchmark == "204YY") or\
             (filename == "578_2010-06-29.dat" and benchmark == "T3973") or\
             (filename == "579_2010-06-29.dat" and benchmark == "T3973") or\
             (filename == "578_2010-07-16.dat" and benchmark == "79-511") or\
             (filename == "578_2012-10-24.dat" and benchmark == "HVO41") or\
             (filename == "579_2012-10-24.dat" and benchmark == "HVO41") or\
             (filename == "578_2017-04-24.dat" and benchmark == "HVO35") or\
             (filename == "579_2017-04-24.dat" and benchmark == "HVO35") or\
             (filename == "578_2017-04-21.dat" and benchmark == "93YY") or\
             (filename == "579_2017-04-21.dat" and benchmark == "93YY") or\
             (filename == "578_2017-04-21.dat" and benchmark == "HVO25") or\
             (filename == "579_2017-04-21.dat" and benchmark == "HVO25") or\
             (filename == "578_2017-04-20.dat" and benchmark == "113YY") or\
             (filename == "579_2017-05-03.dat" and benchmark == "113YY") or\
             (filename == "578_2015-09-15.dat" and benchmark == "HVO49") or\
             (filename == "579_2015-09-15.dat" and benchmark == "HVO49"):
            continue

          # Update the dictionary
          collector[campaign][instrument][benchmark] = (gravity, sd)


  # Start, end, scale of graph
  combinations = [
    ("2009-Dec", "2010-Jun", 80),
    ("2010-Jun", "2011-Mar", 100),
    ("2011-Mar", "2012-Jun", 250),
    ("2012-Jun", "2012-Nov", 200),
    ("2012-Nov", "2015-Sept", 200),
    ("2015-Sept", "2017-Apr", 180)
  ]

  for instrument in ["578", "579"]:
    for combination in combinations:
      compareMap(collector, instrument, *combination)
      #compareDistribution(collector, instrument, *combination)
