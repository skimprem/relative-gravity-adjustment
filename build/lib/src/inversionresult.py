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

  def __init__(self, parent, degree, anchor, drift, x, y, s, dg, vardg, residuals, changes, stations, chi, tare, dtare, tm):

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
    self.tare = tare or 0
    self.dtare = dtare or 0
    self.tm = tm


  def save(self, filepath, relativeTo=None):

    """
    property InversionResult.save
    Saves the inversion results to a file on disk
    """

    # Create the compiled array to be stored
    data = np.array([
      self.changes,
      np.round(self.mdg),
      np.round(self.stddg, 2),
      np.full(len(self.changes), self.anchor)
    ]).T

    header = "\n".join([
      "# CG5, CG6 Relative Gravity Adjustment Export",
      "# Version: %s" % __VERSION__,
      "# Created: %s" % datetime.utcnow(),
      "# Reduced Chi Squared: %s" % self.chi,
      "# Polynomial Degree: %s" % self.degree,
      "# Linear Drift Rate: %s" % self.getDriftRate(),
      "# Tare Index: %s (%sµGal)" % (self.tare, np.round(self.dtare)),
      "Benchmark\tGravity (µGal)\tSD (µGal)\tAnchor"
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
    return "Drift Rate: %sµGal/day" % int(round(86400 * self.drift[1]))


  @property
  def instrument(self):
    try:
      return self.filename.split("_")[0]
    except:
      return "Unknown"

  @property
  def date(self):
    try:
      return self.filename.split("_")[1].split(".")[0]
    except:
      return "Unknown"

  
  def plot(self, filepath, removeDrift=False):

    """
    def InversionResult.plot
    Plots the inversion result with drift / solution
    """

    # Styling
    plt.style.use("seaborn")

    # Fetch real time
    rx = np.array(self.df["Date_Time"])
    ry = self.y

    # Sample the polynomial
    polyy = np.polyval(self.drift, self.x)

    plt.plot(rx[0], np.nan, label=r"$\bf{Legend}$", linewidth=0)

    # Removing drift or not
    if not removeDrift:
      plt.plot(rx, polyy, color="red", linestyle="dashed", label=self.getDriftRate())
    else:
      ry -= polyy
      plt.plot(rx, np.zeros(len(rx)), color="red", linestyle="dashed", label=self.getDriftRate())

    if self.tare:
      tx = rx[np.array(self.tm, dtype=bool)]
      ty = ry[np.array(self.tm, dtype=bool)] + self.dtare
      plt.scatter(tx, ty, edgecolor="red", color="white", linewidth=1, linestyle="dotted", label="Restored Tare: %sμGal" % int(np.round(self.dtare)), zorder=10)

    # Plot the anchor
    self.plotGroup(self.anchor, "%s (Anchor)" % self.anchor, "white")

    plt.scatter(rx[0], np.mean(ry), s=0, label=r"$\bf{\Delta g\ Stations\ (μGal)}$")

    # Plot the other groups
    # Go over all solutions for the stations
    colors = cm.rainbow(np.linspace(0, 1, len(self.changes)))

    for (station, dg, w, color) in zip(self.changes, self.mdg, self.stddg, colors):
      label = "%s (%s±%s)" % (station, int(round(dg)), int(round(2 * w)))
      self.plotGroup(station, label, color)

    plt.legend(frameon=True, loc="center left", bbox_to_anchor=(1, 0.5))
    plt.title("Relative Gravity Adjustment: Inversion Results \n Instrument %s (%s)" % (self.instrument, self.date))
    plt.xlabel("Timestamp (UTC)")
    plt.ylabel("Relative Gravity (μGal)")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

    plt.savefig(filepath, bbox_inches="tight")
    plt.close()
