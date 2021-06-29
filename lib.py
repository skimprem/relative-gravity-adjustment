import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

from matplotlib import cm

class DataWrapper():

  """
  Class DataWrapper
  Wraps the read CG5, CG6 data in a pandas dataframe
  """

  def __init__(self, df):

    # Wrap the dataframe
    self.df = df


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


  def setupPolynomialDesign(self, degree, x, y):

    """
    def setupPolynomialDesign.setupPolynomialDesign
    Sets up the polynomial design matrix (up to 3rd order)
    """

    # Support until third poly: can trivially extend to higher orders
    if degree == 1:
      return np.array([x, np.ones(y.size)])
    elif degree == 2:
      return np.array([np.square(x), x, np.ones(y.size)])
    elif degree == 3:
      return np.array([np.power(x, 3), np.square(x), x, np.ones(y.size)])


  def __invert(self, G, W, y):

    """
    def DataWrapper.__invert
    Inverts the observations to the model vector & uncertainties using WLS inversion
    """

    N = np.linalg.inv(G.T @ W @ G)
    # Get the inversion results
    lsq = N @ G.T @ W @ y
    # Reduced chi-squared (degrees of freedom = number of observations - number of model)
    dof = y.size - np.size(G, 1)
    # These are the residuals from the model
    residuals = y - (G @ lsq)
    # Calculate chi squared
    chi = residuals.T @ W @ residuals
    # Variance of unit weight multiplied by N
    res = (chi / dof) * N
    # Extract the standard deviations
    std = np.sqrt(np.diag(res))

    return lsq, std, residuals


  def invert(self, degree, anchor=None):

    """
    def DataWrapper.invert
    Main routine for the inversion
    """

    # Read data from file
    stations, x, y, s = self.extract()

    # First measurement is the anchor
    if anchor is None:
      anchor = self.getAnchor()

    # Get a list of unique stations and remove the anchor itself
    # These are the stations we find gravity solutions for
    unique_stations = np.array(sorted(set(stations)))
    changes = unique_stations[unique_stations != anchor]

    # Polynomial design matrix
    Gpoly = self.setupPolynomialDesign(degree, x, y)
    # The gravity design matrix
    Gdg = self.setupGravityDesign(stations, changes)
    # Combine polynomial and gravity design matrices
    G = np.hstack((Gpoly.T, Gdg))

    # Weight matrix
    W = np.diag(np.reciprocal(np.square(s)))

    # Invert to model parameters & uncertainties
    lsq, std, residuals = self.__invert(G, W, y)

    # These are the drift parameters
    mbeta = lsq[:degree + 1]
    # Gravity differences & uncertainties
    mdg = lsq[degree + 1:]
    stddg = std[degree + 1:]

    # Eliminate the gravity differences for plotting
    y -= Gdg @ mdg

    return InversionResult(self.df, degree, anchor, mbeta, x, y, s, mdg, stddg, residuals, changes, stations)


class InversionResult():

  """
  Class InversionResult
  Container for results that come from the gravity adjustment inversion
  """

  def __init__(self, df, degree, anchor, drift, x, y, s, dg, vardg, residuals, changes, stations):

    self.df = df
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


  def plot(self):

    """
    def InversionResult.plot
    Plots the inversion result with drift / solution
    """

    plt.style.use("seaborn")

    # Fetch real time and plot the polynomial
    rx = np.array(self.df["Date_Time"])
    plt.plot(rx, np.polyval(self.drift, self.x), color="red", linestyle="dashed")

    # Plot the anhor
    anchor = self.anchor

    idx = self.stations == anchor
    xb = rx[idx]
    yb = self.y[idx]
    sb = self.s[idx]

    plt.scatter(xb, yb, label="%s (Anchor)" % anchor, edgecolor="black", linewidth=1, zorder=3, color="white")
    plt.errorbar(xb, yb, yerr=sb, uplims=True, lolims=True, zorder=2, color="black", linewidth=1, fmt="o", capsize=2)

    colors = cm.rainbow(np.linspace(0, 1, len(self.changes)))

    # Go over all solutions for the stations
    for (station, dg, w, c) in zip(self.changes, self.mdg, self.stddg, colors):

      # Get all data that belongs to one station
      idx = self.stations == station
      xb = rx[idx]
      yb = self.y[idx]
      sb = self.s[idx]

      label = "%s (%i±%i)" % (station, np.round(dg), np.round(2 * w))

      plt.scatter(xb, yb, label=label, edgecolor="black", linewidth=1, zorder=3, color=c)
      plt.errorbar(xb, yb, yerr=sb, uplims=True, lolims=True, zorder=2, color="black", linewidth=1, fmt="o", capsize=2)

    plt.legend(frameon=True)
    plt.title("Relative Gravity Adjustment: Inversion Results")
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
    else:
      raise ValueError("Input version must be CG5 or CG6")


  @classmethod
  def readCG6Dat(self, filepath):
  
    """
    def DataLoader.readCG6Dat
    Reads CG6 data from disk
    """

    df = pd.read_csv(filepath, comment="/", delimiter="\t", parse_dates=[[1, 2]])
  
    return DataWrapper(df)
  

  @classmethod
  def readCG5Dat(self, filepath):
  
    """
    def DataLoader.readCG5Dat
    Reads CG5 data and converts it to CG6 format to use the same DataWrapper class
    """

    # Change header to match CG6
    header = ["Date_Time", "line", "Station", "altitude", "CorrGrav", "StdErr", "tiltx", "tilty", "temperature", "tide", "duration", "rej", "dec", "terrain"]
  
    # Make some modifications to the header
    df = pd.read_csv(filepath, comment="/", delimiter="\s+", parse_dates=[[14, 11]])
    df.columns = header
    # Calculate the stderr (CG5 has stdev)
    df["StdErr"] = df["StdErr"] / np.sqrt(df["duration"] - df["rej"])
  
    return DataWrapper(df)
