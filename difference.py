import os
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pyproj import Proj

def load(campaign, instrument, filename):

  return pd.read_csv(os.path.join("results", campaign, instrument, filename), delimiter="\t", comment="#")


def compareMap(collector, instrument, one, two):

  plt.style.use("seaborn")
  plt.figure()

  # Create UTM projection
  myProj = Proj("+proj=utm +zone=5n, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  save = "%s %s" % (one, two)

  plt.title("Change in Gravity between %s and %s (%s)" % (one, two, instrument))

  one = collector[one][instrument]
  two = collector[two][instrument]

  locations = pd.read_csv("locations/stations.csv", delimiter="\t")
  xProj, yProj = myProj(locations["Longitude"], locations["Latitude"])

  latitudes = list()
  longitudes = list()
  colors = list()

  for station, longitude, latitude in zip(locations["BM"], xProj, yProj):

    if not station in one or not station in two:
      continue

    # Calculate the difference
    dg = two[station][0]  - one[station][0]

    latitudes.append(latitude)
    longitudes.append(longitude)
    colors.append(dg)

  # Create map
  cmap = cm.get_cmap("seismic", 21)
  plt.scatter(longitudes, latitudes, c=colors, cmap=cmap, edgecolor="black", linewidth=1)
  scale = np.max(np.abs(colors))

  bar = plt.colorbar()
  bar.mappable.set_clim(-scale, scale)
  plt.tight_layout()
  plt.xlim(257000, 266000)
  plt.ylim(2139000, 2152000)
  plt.xlabel("UTM Easting (m)")
  plt.ylabel("UTM Northing (m)")
  plt.gca().set_aspect("equal", adjustable="box")

  plt.savefig("figures/map/%s.png" % save, bbox_inches="tight")
  plt.close()


def compareDistribution(collector, instrument, one, two):

  plt.style.use("seaborn")

  plt.title("Change in Gravity between %s and %s (%s)" % (one, two, instrument))

  save = "%s %s" % (one, two)

  one = collector[one][instrument]
  two = collector[two][instrument]

  locations = pd.read_csv("locations/stations.csv", delimiter="\t")

  for station in locations["BM"]:

    if not station in one or not station in two:
      continue

    # Calculate the difference
    dg = two[station][0]  - one[station][0]
    # Two sigma
    std = 2 * np.sqrt(two[station][1] ** 2 + one[station][1] ** 2)

    if (np.abs(dg) - std) < 0 or dg == 0:
      color = "white"
    elif dg > 0:
      color = "red"
    elif dg < 0:
      color = "blue"

    plt.bar(station, dg, edgecolor="black", zorder=1, yerr=std, capsize=4, linewidth=1, color=color)

  plt.xlabel("Gravity Benchmark")
  plt.ylabel("Change in Gravity (µGal)")
  plt.tight_layout()
  plt.xticks(rotation=90)

  plt.savefig("figures/dist/%s.png" % save, bbox_inches="tight")
  plt.close()


if __name__ == "__main__":

  collector = {}

  for campaign in os.listdir("results"):

    collector[campaign] = {"578": {}, "579": {}}

    for instrument in os.listdir(os.path.join("results", campaign)):

      for filename in os.listdir(os.path.join("results", campaign, instrument)):

        if filename.startswith("."):
          continue

        df = load(campaign, instrument, filename)

        for _, row in df.iterrows(): 
          benchmark = row["Benchmark"]
          gravity = row["Gravity (µGal)"]
          sd = row["SD (µGal)"]
          collector[campaign][instrument][benchmark] = (gravity, sd)

  compareDistribution(collector, "578", "2009-Dec", "2010-Jun")
  compareDistribution(collector, "578", "2010-Jun", "2011-Mar")
  compareDistribution(collector, "578", "2011-Mar", "2012-Jun")
  compareDistribution(collector, "578", "2012-Jun", "2012-Nov")
  compareDistribution(collector, "578", "2012-Nov", "2015-Sept")
  compareDistribution(collector, "578", "2015-Sept", "2017-Apr")
