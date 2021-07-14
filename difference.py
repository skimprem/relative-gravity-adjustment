import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

from matplotlib.legend_handler import HandlerTuple
from pyproj import Proj
from matplotlib import cm

def relativeTo(one, onesd, two, twosd):

  return one + two, np.sqrt(onesd ** 2 + twosd ** 2)


def loadResult(campaign, instrument, filename):

  """
  def ResultWrapper.loadResult
  Loads results from the results directory
  """

  filepath = os.path.join("results", campaign, instrument, filename)

  return pd.read_csv(
    filepath,
    delimiter="\t",
    comment="#"
  )

class ResultWrapper():

  myProj = Proj("+proj=utm +zone=5n, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  def __init__(self, data):

    self.data = data
    self.dem = self.loadDEM("dem/kilauea.nc")


  def loadDEM(self, filepath):

    dataset = xr.open_dataset(filepath)
    grid = np.meshgrid(dataset.variables["lon"][:], dataset.variables["lat"][:])
    x, y = self.myProj(*grid)
    z = np.array(dataset.variables["elev"][:])
    # Set invalid numbers to zero (sea)
    z[np.isnan(z)] = 0

    return (x, y, z)


  def plotDEM(self, levels):
  
    """
    def ResultWrapper.plotDEM
    Plots a Kilaueau DEM 
    """
  
    # DEM contouring
    plt.contourf(
      *self.dem,
      levels=levels,
      cmap=cm.Greys,
      zorder=0
    )
  
    CS = plt.contour(
      *self.dem,
      levels=levels,
      colors="black",
      linewidths=1,
      zorder=0
    )
  
    # Set linewidth 
    plt.gca().grid(b=True, which="major", color="white", linewidth=0.5)


  def get(self, which, instrument):

    return self.data[which][instrument]

  def compareMap(self, instrument, one, two, scale):
  
    """
    def ResultWrapper.compareMap
    Plots differences between two campaigns on a map
    """
  
    # Limit of map
    XLIM = (257000, 266000)
    YLIM = (2139000, 2152000)
  
    plt.style.use("seaborn")
  
    self.plotDEM(21)
  
    filename = "%s %s %s" % (one, two, instrument)
    plt.title("Change in Gravity Difference with Anchor P1 \n Between %s and %s \n Instrument %s" % (one, two, instrument))
  
    dataOne = self.get(one, instrument)
    dataTwo = self.get(two, instrument)
  
    locations = pd.read_csv("locations/stations.csv", delimiter="\t")
    xProj, yProj = self.myProj(locations["Longitude"], locations["Latitude"])
  
    benchmarks = list()
    missing = list()
    colors = list()
  
    for station, longitude, latitude in zip(locations["BM"], xProj, yProj):
  
      if station == "P1":
        b = plt.scatter(longitude, latitude, marker="*", s=200, color="yellow", linewidth=1, edgecolor="black")
        continue
      elif not station in dataOne or not station in dataTwo:
        missing.append((longitude, latitude))
        text = plt.gca().annotate(station, (longitude, latitude + 75), color="white", fontsize=2, ha="center")
        text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'),
                               path_effects.Normal()])
        continue
  
      # Calculate the difference
      dg, std = self.difference(dataOne, dataTwo, station)
  
      # Read height
      dg += self.getHeight(one, two, station)

      if np.isnan(dg):
        missing.append((longitude, latitude))
        text = plt.gca().annotate(station, (longitude, latitude + 75), color="white", fontsize=2, ha="center")
        text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'),
                               path_effects.Normal()])
        continue
  
      benchmarks.append((longitude, latitude))
      colors.append(dg)
      text = plt.gca().annotate(station, (longitude, latitude + 150), color="white", fontsize=2, ha="center")
      text.set_path_effects([path_effects.Stroke(linewidth=1, foreground='black'),
                             path_effects.Normal()])
  
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


  def difference(self, one, two, station):
  
    """
    def ResultWrapper.difference
    Reads heights from a .CSV file
    """
  
    # Difference between two campaigns
    dg = two[station][0]  - one[station][0]
  
    # 2 times standard deviation from differences sqrt(σ1^2 + σ2^2)
    std = 2 * np.sqrt(two[station][1] ** 2 + one[station][1] ** 2)
  
    return dg, std


  def save(self, filepath):

    header = ["BM, UTM (E), UTM(N), ELEV"]

    combinations = [
      ("2009-Dec", "2010-Jun"),
      ("2010-Jun", "2011-Mar"),
      ("2011-Mar", "2012-Jun"),
      ("2012-Jun", "2012-Nov"),
      ("2012-Nov", "2015-Sept"),
      ("2015-Sept", "2017-Apr")
    ]

    locations = pd.read_csv("locations/stations.csv", delimiter="\t")
    xProj, yProj = self.myProj(locations["Longitude"], locations["Latitude"])
    data = [locations["BM"], xProj, yProj, locations["Elevation"]]

    for combination in combinations:
      for instrument in ["578", "579"]:

        dataOne = self.get(combination[0], instrument)
        dataTwo = self.get(combination[1], instrument)
 
        columnGravity = list()
        columnStd = list()

        for station, longitude, latitude in zip(locations["BM"], xProj, yProj):

          if station == "P1" or (station not in dataOne or station not in dataTwo):
            columnGravity.append(np.nan)
            columnStd.append(np.nan)
            continue

          dg, std = self.difference(dataOne, dataTwo, station)
          columnGravity.append(dg)
          columnStd.append(np.round(std))

        header.append("DG %s (%s %s)" % (instrument, *combination))
        header.append("STD %s (%s %s)" % (instrument, *combination))

        data.append(columnGravity)
        data.append(columnStd)

    np.savetxt(filepath,
               np.array(data).T,
               delimiter=",",
               comments="",
               fmt="%s",
               header=",".join(header))

  def getHeight(self, one, two, station):
  
    """
    def ResultWrapper.getHeight
    Reads heights from a .CSV file
    """
  
    df = pd.read_csv(
      "heights.csv",
      delimiter="\t",
      comment="#"
    )
  
    # Undo the effect of height change (negative height will give positive gravity)
    return 3.086 * df[df["Benchmark"] == station].iloc[0]["%s %s" % (one, two)]
  
  
  def compareDistribution(self, instrument, one, two, scale):
  
    """
    def ResultWrapper.compareDistribution
    Plots distribution of data against the benchmarks
    """
  
    plt.style.use("seaborn")
  
    plt.figure(figsize=(5, 8))
    plt.title("Change in Gravity Difference with Anchor P1 \n Between %s and %s \n Instrument %s" % (one, two, instrument))
  
    # Output filename
    filename = "%s %s %s" % (one, two, instrument)
  
    # Fetch the results
    dataOne = self.get(one, instrument)
    dataTwo = self.get(two, instrument)
  
    locations = pd.read_csv("locations/stations.csv", delimiter="\t")
  
    for station in sorted(locations["BM"]):
  
      if station == "P1":
        continue
  
      if not station in dataOne or not station in dataTwo:
        plt.barh(station, 0, edgecolor="black", xerr=0, capsize=0, linewidth=0, color=color, error_kw=dict(capthick=0))
        continue
  
      dg, std = self.difference(dataOne, dataTwo, station)
  
      # Height correction
      dg += self.getHeight(one, two, station)
  
      if (np.abs(dg) - std) < 0 or dg == 0:
        color = "white"
      elif dg > 0:
        color = "red"
      elif dg < 0:
        color = "blue"
  
      plt.barh(station, dg, edgecolor="black", xerr=std, capsize=2, linewidth=1, color=color, error_kw=dict(capthick=1))
  
    plt.xlim(-scale, scale)
    plt.margins(y=0)
    plt.ylabel("Benchmark")
    plt.xlabel("Change in Gravity Difference (µGal)")
    plt.yticks(ha="left", position=(-0.13, 0))
    plt.tight_layout()
  
    plt.savefig("figures/dist/%s.pdf" % filename, bbox_inches="tight")
    plt.close()


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
        result = loadResult(campaign, instrument, filename)

        # Go over each row
        for _, row in result.iterrows(): 

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


  wrapper = ResultWrapper(collector)
  #wrapper.save("compiled.csv")
  #sys.exit(0)

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
      wrapper.compareMap(instrument, *combination)
      #wrapper.compareDistribution(instrument, *combination)
