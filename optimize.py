import pandas as pd
from matplotlib.colors import DivergingNorm
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import numpy as np
import sys
from scipy.optimize import minimize
from pyproj import Proj
from lib.model import getDEM,plotDEM

def getUTMCoordinates():

  # Load the DEM
  DEM = getDEM(area="kilauea")

  # Get parameters
  longitude = np.array(DEM.variables["lon"][:])
  latitude = np.array(DEM.variables["lat"][:])

  # Get the DEM elevations
  elevation = np.array(DEM.variables["elev"][:])
  elevation[np.isnan(elevation)] = 0

  # Interpolate the lat/lng grid to find station elevations
  model = interp2d(longitude, latitude, elevation)

  # Convert lat/lng to UTM
  myProj = Proj("+proj=utm +zone=5n, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  xProj, yProj = myProj(metadata["Longitude"], metadata["Latitude"])

  # Add the elevation
  zProj = []
  for (x, y) in zip(metadata["Longitude"], metadata["Latitude"]):
    zProj.append(model(x, y)[0])
  zProj = np.array(zProj)

  return xProj, yProj, zProj


def __gravity(M, dx, dy, dz):

  """
  def __gravity
  Objective function to be minimised
  """

  G = 6.67408E-11

  rs = dx ** 2 + dy ** 2 + dz ** 2

  return 1E8 * G * M * (dz / np.power(rs, 1.5))


if __name__ == "__main__":

  which = "2009-Dec 2010-Jun"
  #which = "2010-Jun 2011-Mar"
  which = "2011-Mar 2012-Jun"
  which = "2012-Jun 2012-Nov"
  which = "2012-Nov 2015-Sept"
  which = "2015-Sept 2017-Apr"

  # The station metadata
  metadata = pd.read_csv(
    "metadata/stations.csv",
    delimiter="\t"
  )

  xProj, yProj, zProj = getUTMCoordinates()
  
  # Load the gravity differences
  measurements = pd.read_csv(
    "results-578.csv",
    delimiter="\t"
  )

  # Set up the order of coordinates
  xs = []
  ys = []
  zs = []
  vals = []
  
  # Go over each benchmark: same order as projection
  for (BM, x, y, z) in zip(metadata["BM"], xProj, yProj, zProj):
  
    # Skip the anchor
    if BM == "P1":
      P1x = x
      P1y = y
      P1z = z
      continue

    # Benchmarks too close to the lava lake
    if BM in ["68-15"]:
      continue

    # Benchmarks too close to the lava lake
    if BM in ["HOVL-G", "205YY", "HVO41"]:
      continue

    # Get the measurement
    g = measurements[measurements["Station"] == BM][which].iloc[0]

    if np.isnan(g):
      continue

    xs.append(x)
    ys.append(y)
    zs.append(z)
    vals.append(g)
  
  # Convert to numpy arrays
  xs = np.array(xs)
  ys = np.array(ys)
  # Elevation is negative z
  zs = -np.array(zs)
  measuredGravity = np.array(vals)
  
  # Initial guess
  x0 = np.array([1E8, 260000, 2147000, 3000])

  def gravityObjective(x):
  
    """
    def gravityObjective
    Objective function to minimize
    """
  
    # Calculate coordinate differences
    dx = x[1] - xs
    dy = x[2] - ys
    dz = x[3] - zs
  
    predictedGravity = __gravity(x[0], dx, dy, dz)
  
    # Correct for the effect at P1
    p1dx = x[1] - P1x
    p1dy = x[2] - P1y
    p1dz = x[3] - P1z
  
    predictedGravity -= __gravity(x[0], p1dx, p1dy, p1dz)
  
    # Objective: minimize the residual squared sum
    return np.sum(np.square(measuredGravity - predictedGravity))

  # Minimize
  solution = minimize(
    gravityObjective,
    x0,
    method="Nelder-Mead",
    options={"maxiter": 10000}
  )
  
  predictedGravity = []
  
  for (x, y, z) in zip(xs, ys, zs):
  
    dx = solution.x[1] - x
    dy = solution.x[2] - y
    dz = solution.x[3] - z
  
    predictedGravity.append(__gravity(solution.x[0], dx, dy, dz))
  
  predictedGravity = np.array(predictedGravity)
    
  plotDEM(area="kilauea")
  label = "{:.2e}kg".format(solution.x[0]) + " & depth = " + str(np.round(solution.x[3])) + "m"
  plt.scatter(xs, ys, c=measuredGravity, cmap=cm.seismic, edgecolor="black", linewidth=1, norm=DivergingNorm(0))
  plt.colorbar()

  plt.scatter(solution.x[1], solution.x[2], label=label, edgecolor="black", marker="*", color="yellow", s=500, linewidth=1)
  plt.legend()
  plt.show()
