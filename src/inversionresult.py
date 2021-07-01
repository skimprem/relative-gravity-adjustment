import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

from matplotlib import cm
from datetime import datetime

__VERSION__ = "0.0.1"

class InversionResult():

  """
  Class InversionResult
  Container for results that come from the gravity adjustment inversion
  """

  def __init__(self, parent, degree, anchor, drift, x, y, s, dg, vardg, residuals, changes, stations, chi):

    self.start = parent.start
    self.end = parent.end
    self.df = parent.df
    self.filename = parent.filename
    self.x = x
    self.s = s
    self.y = y
    self.degree = degree
    self.anchor = anchor
    self.drift = np.poly1d(drift)
    self.mdg = dg
    self.stddg = vardg
    self.residuals = residuals
    self.changes = changes
    self.stations = stations
    self.chi = chi


  def save(self, filepath):

    """
    property InversionResult.save
    Saves the inversion results to a file on disk
    """

    # Create the compiled array to be stored
    data = np.array([self.changes, np.round(self.mdg), np.round(self.stddg, 2)]).T

    header = "\n".join([
      "# CG5, CG6 Relative Gravity Adjustment Export",
      "# Version: %s" % __VERSION__,
      "# Created: %s" % datetime.utcnow(),
      "# Reduced Chi Squared: %s" % self.chi,
      "# Anchor: %s" % self.anchor,
      "# Polynomial Degree: %s" % self.degree,
      "# Linear Drift Rate: %s" % self.getDriftRate(),
      "Benchmark\tGravity (µGal)\tSD (µGal)"
    ])

    # Save to file
    np.savetxt(filepath,
               data,
               delimiter="\t",
               comments="",
               fmt="%s",
               header=header)


  @property
  def differences(self):

    """
    property InversionResult.differences
    Returns a list of tuples with the gravity differences
    """

    tuples = list()

    for (x) in zip(self.changes, self.mdg, self.stddg):
      tuples.append(x)

    return {"anchor": self.anchor, "differences": tuples}


  def plotResiduals(self):

    """
    def InversionResult.plotResiduals
    Function to plot the data residuals from the model: should be normally distributed around 0
    """

    plt.style.use("seaborn")

    plt.title("Model Residuals")
    plt.xlabel("Gravity Residual (μGal)")
    plt.ylabel("Probability Density")
    plt.hist(self.residuals, 50, density=True, edgecolor="black", linewidth=1)
    plt.show()


  def plotGroup(self, benchmark, label, color):

    """
    def InversionResult.plotGroup
    Plots a single group
    """

    # Get indices of the right group
    idx = self.stations == benchmark 
    xb = self.df["Date_Time"][idx]
    yb = self.y[idx]
    sb = self.s[idx]

    plt.scatter(xb, yb, label=label, edgecolor="black", linewidth=1, zorder=3, color=color)
    plt.errorbar(xb, yb, yerr=sb, uplims=True, lolims=True, zorder=2, color="black", linewidth=1, fmt="o", capsize=2)


  def getDriftRate(self):

    """
    def InversionResult.getDriftRate
    Returns and formats the drift rate from the recovered polynomial
    """

    # Disregard the intercept
    return "%sµGal/day" % int(round(86400 * self.drift[1]))


  def plot(self, removeDrift=False):

    """
    def InversionResult.plot
    Plots the inversion result with drift / solution
    """

    plt.style.use("seaborn")

    # Fetch real time
    rx = np.array(self.df["Date_Time"])

    # Sample the polynomial
    polyy = np.polyval(self.drift, self.x)

    # Removing drift or not
    if not removeDrift:
      plt.plot(rx, polyy, color="red", linestyle="dashed", label=self.getDriftRate())
    else:
      self.y -= polyy
      plt.plot(rx, np.zeros(len(rx)), color="red", linestyle="dashed", label=self.getDriftRate())

    # Plot the anchor
    self.plotGroup(self.anchor, "%s (Anchor)" % self.anchor, "white")

    # Plot the other groups
    # Go over all solutions for the stations
    colors = cm.rainbow(np.linspace(0, 1, len(self.changes)))

    for (station, dg, w, color) in zip(self.changes, self.mdg, self.stddg, colors):
      label = "%s (%s±%sµGal)" % (station, int(round(dg)), int(round(2 * w)))
      self.plotGroup(station, label, color)

    plt.legend(frameon=True)
    plt.title("Relative Gravity Adjustment: Inversion Results (%s)" % self.filename)
    plt.xlabel("Timestamp")
    plt.ylabel("Relative Gravity (μGal)")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))

    plt.show()
