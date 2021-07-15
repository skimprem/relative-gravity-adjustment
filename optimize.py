import pandas as pd
from matplotlib.colors import DivergingNorm
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import numpy as np
import sys
from scipy.optimize import minimize
from pyproj import Proj

def __gravity(M, dx, dy, dz):

  """
  def __gravity
  Objective function to be minimised
  """

  G = 6.67408E-11

  rs = dx ** 2 + dy ** 2 + dz ** 2
  return 1E8 * G * M * (dz / np.power(rs, 1.5))


def optimize(which):

  # Load the gravity differences
  measurements = pd.read_csv(
    "compiled.csv",
    delimiter="\t"
  )

  measurements = measurements[~np.isnan(measurements[which])]
  measurements = measurements[~measurements["BM"].isin(["P1", "HOVL-G", "205YY", "HVO41"])]

  # Convert to numpy arrays
  # Elevation is negative z
  xs = np.array(measurements["UTM (E)"])
  ys = np.array(measurements["UTM (N)"])
  zs = -np.array(measurements["ELEV"])
  measuredGravity = np.array(measurements[which])

  # The initial guess for the parameters (M, x, y, z) for source position
  x0 = np.array([1E10, 260000, 2147000, 3000])

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
    p1dx = x[1] - 258376
    p1dy = x[2] - 2150798
    p1dz = x[3] - -1202.0
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

  return solution.x

if __name__ == "__main__":

  """
  def __main__
  Optimization
  """

  combinations = (
    "2009-Dec 2010-Jun",
    "2010-Jun 2011-Mar",
    "2011-Mar 2012-Jun",
    "2012-Jun 2012-Nov",
    "2012-Nov 2015-Sept",
    "2015-Sept 2017-Apr"
  )

  save = list()
  for instrument in ("578", "579"):
    for combination in combinations:
      which = "DG %s (%s)" % (instrument, combination)
      solution = optimize(which)
      save.append([which] + list(solution))

  df = pd.DataFrame(save)
  df.columns = ["Campaign", "Mass", "UTM (E)", "UTM (N)", "ELEV"]

  df.to_csv(
    "sources.csv",
    sep="\t",
    index=False
  )
