import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

from longman.longman import TideModel

from scipy.interpolate import interp1d
from matplotlib import cm
from datetime import datetime, timedelta

__VERSION__ = "0.0.1"

class DataWrapper():

  """
  Class DataWrapper
  Wraps the read CG5, CG6 data in a pandas dataframe
  """

  def __init__(self, df, filepath):

    # Wrap the dataframe
    self.df = df[~np.isnan(df["StdErr"])]
    self.filename = os.path.basename(filepath)


  def getAnchor(self):

    """
    def DataWrapper.getAnchor
    Returns the station name of the first measurement: if not specified this is the default anchor
    """

    return str(self.df["Station"][0])


  def extract(self):

    """
    def DataWrapper.extract
    Extracts the necessary information from the data file
    """

    benchmarks = np.array(self.df["Station"], dtype=str)

    x = np.array(self.df["Date_Time"])

    # Scale to microGal
    y = 1E3 * np.array(self.df["CorrGrav"])
    s = 1E3 * np.array(self.df["StdErr"])

    # Zero out first measurement: OK since relative
    x = 86400 * (mdates.date2num(x) - mdates.date2num(x[0]))
    y -= y[0]

    return benchmarks, x, y, s


  @property
  def start(self):

    """
    def DataWrapper.start
    Returns the start time of the data set
    """

    return mdates.date2num(self.df["Date_Time"].iloc[0])


  @property
  def end(self):

    """
    def DataWrapper.start
    Returns the start time of the data set
    """

    return mdates.date2num(self.df["Date_Time"].iloc[-1])


  def setupGravityDesign(self, stations, changes):

    """
    def DataWrapper.setupGravityDesign
    Sets up the gravity part of the design matrix for the inversion
    """

    # Broadcast to set up the design matrix
    # We want to repeat the lists to easy mask them out and get the design matrix for free
    matrix_stations = np.broadcast_to(stations, (changes.size, stations.size))
    matrix_changes = np.broadcast_to(changes, (stations.size, changes.size))

    # The gravity design matrix (binary)
    return np.array(matrix_stations.T == matrix_changes, dtype=np.int)


  def setupPolynomialDesign(self, degree, x):

    """
    def setupPolynomialDesign.setupPolynomialDesign
    Sets up the polynomial design matrix
    """

    if degree == 0:
      raise ValueError("Polynomial degree cannot be 0.")

    matrix = [np.ones(x.size), x]

    if degree == 1:
      return np.flip(np.array(matrix), axis=0)

    # Add higher orders
    for deg in range(2, degree + 1):
      matrix.append(np.power(x, deg))

    return np.flip(np.array(matrix), axis=0)


  def __invert(self, G, W, y):

    """
    def DataWrapper.__invert
    Inverts the observations to the model vector & uncertainties using WLS inversion
    """

    N = np.linalg.inv(G.T @ W @ G)
    # Get the inversion results
    lsq = N @ G.T @ W @ y
    # Reduced chi-squared (degrees of freedom = number of observations - number of model)
    dof = y.size - np.size(G, 1) - 1
    # These are the residuals from the model
    residuals = y - (G @ lsq)
    # Calculate chi squared
    rchi = (residuals.T @ W @ residuals) / dof
    # Variance of unit weight multiplied by N
    res = rchi * N
    # Extract the standard deviations
    std = np.sqrt(np.diag(res))

    return lsq, std, residuals, rchi


  def plotTide(self):

    """
    def DataWrapper.plotTide
    Plots tidal comparison between Longman, Instrument Default, ETERNA 3.4
    """

    plt.style.use("seaborn")

    model = self.getETERNA(self.df["Date_Time"].iloc[0], self.df["Date_Time"].iloc[-1])
    x = self.df["Date_Time"]
    y = self.getLongman(x)
    plt.plot(x, y, label="Longman")
    plt.plot(x, -model(mdates.date2num(x)), label="ETERNA 3.4")
    plt.plot(x, 1E3 * self.df["TideCorr"], label="Default")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    plt.legend()
    plt.show()

  def correctTide(self, y, which):

    x = self.df["Date_Time"]

    # Eliminate the default correction
    y -= 1E3 * self.df["TideCorr"]

    latitude = 19.40840
    longitude = -155.28385
    height = 1000

    if which == "ETERNA":
      return y - self.getETERNA(latitude, longitude, height, x)(mdates.date2num(x))
    elif which == "Longman":
      return y + self.getLongman(latitude, longitude, height, x)
    else:
      raise ValueError("Unknown tidal correction.")


  def invert(self, degree, anchor=None, tide="default"):

    """
    def DataWrapper.invert
    Main routine for the inversion
    """

    # Read data from file
    stations, x, y, s = self.extract()

    # Correct for the tide using ETERNA
    if tide != "default":
      y = self.correctTide(y, tide)

    # First measurement is the anchor
    if anchor is None:
      anchor = self.getAnchor()

    # Get a list of unique stations and remove the anchor itself
    # These are the stations we find gravity solutions for
    unique_stations = np.array(sorted(set(stations)))
    changes = unique_stations[unique_stations != anchor]

    # Polynomial design matrix
    Gpoly = self.setupPolynomialDesign(degree, x)
    # The gravity design matrix
    Gdg = self.setupGravityDesign(stations, changes)
    # Combine polynomial and gravity design matrices
    G = np.hstack((Gpoly.T, Gdg))

    # Weight matrix
    W = np.diag(np.reciprocal(np.square(s)))

    # Invert to model parameters & uncertainties
    lsq, std, residuals, chi = self.__invert(G, W, y)

    # These are the drift parameters
    mbeta = lsq[:degree + 1]
    # Gravity differences & uncertainties
    mdg = lsq[degree + 1:]
    stddg = std[degree + 1:]

    # Eliminate the gravity differences for plotting
    y -= Gdg @ mdg

    return InversionResult(self, degree, anchor, mbeta, x, y, s, mdg, stddg, residuals, changes, stations, chi)


  def getLongman(self, latitude, longitude, height, times):
  
    """
    Returns Longman tidal model
    Taken from https://github.com/jrleeman/LongmanTide
    """
  
    model = TideModel()  # Make a model object
  
    solution = list()

    # Not vectorized
    for time in times:
      _, _, g = model.solve_longman(latitude, longitude, height, time)
      solution.append(g)
  
    # Scale results to microGal
    return 1E3 * np.array(solution)


  def getETERNA(self, latitude, longitude, height, times):
  
    try:
      import pygtide
    except ImportError:
      raise ValueError("Tide correction can only be ETERNA if pygtide is installed.")

    start = times.iloc[0]
    end = times.iloc[-1]

    # create a PyGTide object
    pt = pygtide.pygtide(msg=False)
   
    # Start, duration and sample rate of the model
    duration = max(48, (end - start).days * 24)
    samplerate = 60
  
    # Predict the tides using ETERNA
    pt.predict(
      latitude,
      longitude,
      height,
      start,
      duration,
      samplerate,
      tidalcompo=0
    )
   
    # retrieve the results as dataframe
    data = pt.results()
  
    x = mdates.date2num(data["UTC"])
    y = 1E-1 * data["Signal [nm/s**2]"]

    return interp1d(x, y, kind="cubic")


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
    Saves inversion results to a file
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
      "Benchmark\tGravity\tSD"
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


class DataLoader():

  """
  Class DataLoader
  Loads CG5 or CG6 data from disk
  """

  def __init__(self):
    pass


  @classmethod
  def load(self, ver, filepath):

    if ver == "CG5":
      return self.readCG5Dat(filepath)
    elif ver == "CG6":
      return self.readCG6Dat(filepath)
    elif ver == "USGS":
      return self.readUSGS(filepath)
    else:
      raise ValueError("Input version must be CG5 or CG6")


  @classmethod
  def readCG6Dat(self, filepath):
  
    """
    def DataLoader.readCG6Dat
    Reads CG6 data from disk
    """

    df = pd.read_csv(filepath, comment="/", delimiter="\t", parse_dates=[[1, 2]])

    return DataWrapper(df, filepath)
  

  @classmethod
  def readUSGS(self, filepath):

    df = pd.read_csv(filepath, delimiter="\t", parse_dates=[[15, 12]])
    header = ["Date_Time", "Station", "Latitude", "Longitude", "Altitude", "CorrGrav", "StdErr", "TiltX", "TiltY", "Temp", "TideCorr", "duration", "rej", "Dec. Time+Date", "Terrain", "Accepted"]

    df.columns = header
    df["StdErr"] = df["StdErr"] / np.sqrt((df["duration"] - df["rej"]))

    return DataWrapper(df, filepath)


  @classmethod
  def readCG5Dat(self, filepath):
  
    """
    def DataLoader.readCG5Dat
    Reads CG5 data and converts it to CG6 format to use the same DataWrapper class
    """

    # Change header to match CG6
    header = ["Date_Time", "line", "Station", "altitude", "CorrGrav", "StdErr", "tiltx", "tilty", "temperature", "TideCorr", "duration", "rej", "dec", "terrain"]
  
    df = pd.read_csv(filepath, skiprows=31, delimiter="\s+", header=None, parse_dates=[[14, 11]])
    # Make some modifications to the header
    df.columns = header

    # Calculate the stderr (CG5 has stdev).. It samples at 6Hz but CG6 calculates the standard error like this
    df["StdErr"] = df["StdErr"] / np.sqrt(df["duration"] - df["rej"])
  
    return DataWrapper(df, filepath)
