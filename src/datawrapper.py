import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

from src.tide import TidalModel
from src.oceanloading import OceanLoadingModel
from src.inversionresult import InversionResult

class DataWrapper():

  """
  Class DataWrapper
  Wraps the read CG5, CG6 data in a pandas dataframe
  """

  def __init__(self, inp, df, filepath):

    """
    DataWrapper.__init__
    Initializes a datawrapper
    """

    # Wrap the dataframe
    self.inp = inp
    self.df = self.filter(df)
    self.filename = os.path.basename(filepath)

    # Metadata
    self.locations = None


  def setLocations(self, filename): 

    """
    def DataWrapper.setLocation
    Sets location required for tidal correctio
    """

    import pandas as pd

    self.locations = pd.read_csv(filename, delimiter="\t")


  def getLocation(self, benchmark):

    """
    def DataWrapper.setLocation
    Sets location required for tidal correctio
    """

    row = self.locations[self.locations["BM"] == benchmark]

    if len(row) == 0:
      raise ValueError("No location specified for the requested benchmark (%s)." % benchmark)

    return row["Latitude"].iloc[0], row["Longitude"].iloc[0], 0


  def getAnchor(self):

    """
    def DataWrapper.getAnchor
    Returns the station name of the first measurement: if not specified this is the default anchor
    """

    return str(self.df["Station"].iloc[0])


  def filter(self, df):

    """
    def DataWrapper.filter
    Applies a filter to the data to remove poor data
    """

    # Only keep the accepted values
    try:
      df = df[(df["Accepted"] == 1)]
    except:
      pass

    df = df[(np.abs(df["TiltX"]) < 20) & (np.abs(df["TiltY"]) < 20)]
    df = df[df["duration"] >= 60]

    try:
      df = df[df["rej"] < 5]
    except:
      pass

    return df


  def extract(self):

    df = self.df

    benchmarks = np.array(df["Station"], dtype=str)

    x = np.array(df["Date_Time"])

    # Scale to microGal
    y = 1E3 * np.array(df["CorrGrav"], dtype=float)
    s = 1E3 * np.array(df["StdErr"])

    # Zero out first measurement: OK since relative
    x = 86400 * (mdates.date2num(x) - mdates.date2num(x[0]))
    y -= y[0]

    # Undo the instrument tide correction
    # y -= 1E3 * df["TideCorr"]

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
    return np.array(matrix_stations.T == matrix_changes, dtype=int)


  def setupTareDesign(self, stations, tare):

    Gtare = np.zeros(stations.size, dtype=int)

    # Set the respective tare indices to one
    if tare != 0:
      Gtare[tare:] = 1

    return np.array([Gtare])


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
    dof = y.size - np.size(G, 1) 
    # These are the residuals from the model
    residuals = y - (G @ lsq)
    # Calculate reduced chi squared
    rchi = (residuals.T @ W @ residuals) / dof
    # Variance of unit weight multiplied by N
    res = rchi * N
    # Extract the standard deviations
    std = np.sqrt(np.diag(res))

    return lsq, std, residuals, rchi


  def plotTide(self, location):

    """
    def DataWrapper.plotTide
    Plots tidal comparison between Longman, Instrument Default, ETERNA 3.4
    """

    plt.style.use("seaborn")

    # Create a tidal model
    tideModel = TidalModel(*location)

    x = self.df["Date_Time"]

    # Create an ocean loading model
    loadingModel = OceanLoadingModel("ocl/harmonics/HOVL-G.txt").getOceanLoadingModel(self.df["Date_Time"])

    eternaModel = tideModel.getETERNA(x)

    y = tideModel.getLongman(x)

    plt.plot(x, y, label="Longman")
    plt.plot(x, -eternaModel(mdates.date2num(x)), label="ETERNA 3.4")
    plt.plot(x, -eternaModel(mdates.date2num(x)) + loadingModel(mdates.date2num(x)), label="ETERNA 3.4 + Ocean Loading")
    plt.plot(x, 1E3 * self.df["TideCorr"], label="Default")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    plt.legend()
    plt.show()


  def correctLoading(self, y):

    """
    def DataWrapper.correctLoading
    Corrects for the ocean loading effect
    """

    x = self.df["Date_Time"]

    # Go over all the benchmarks
    for benchmark in set(self.df["Station"]):
      idx = self.df["Station"] == benchmark
      xs = x[idx]
      loadingModel = OceanLoadingModel("ocl/harmonics/%s.txt" % benchmark).getOceanLoadingModel(xs)

      # Ocean loading needs to be ADDED because it is a CORRECTION
      y[idx] += loadingModel(mdates.date2num(xs))

    return y


  def correctTide(self, y, which):

    """
    def DataWrapper.correctTide
    Corrects for the solid earth tide using ETERNA or Longman 1959
    """

    # Must be supplied

    x = self.df["Date_Time"]

    for benchmark in set(self.df["Station"]):
      
      # Fetch the location of the benchmark itself
      location = self.getLocation(benchmark)

      # All indices of this benchmark
      idx = self.df["Station"] == benchmark

      # Get the times of this benchmark
      xs = x[idx]

      # Create a model at the particular location
      model = TidalModel(*location)

      # ETERNA gives effect: SUBTRACT
      if which == "ETERNA":
        y[idx] -=  model.getETERNA(xs)(mdates.date2num(xs))
      # Longman is correction: ADD
      elif which == "Longman":
        y[idx] += model.getLongman(xs)
      else:
        raise ValueError("Unknown tidal correction requested.")

    return y


  def invert(self, degree, anchor=None, tide="default", loading=False, tare=None):

    """
    def DataWrapper.invert
    Main routine for the inversion
    """

    # Read the data from the datawrapper
    stations, x, y, s = self.extract()

    # Correct for the tide
    if tide != "default":
      y = self.correctTide(y, tide)

    # Apply the ocean loading model
    if loading:
      y = self.correctLoading(y)

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

    # Combine polynomial and gravity design matrices: different depending on whether we introduce a tare
    if tare is None:
      Gtare = [None]
      G = np.hstack((Gpoly.T, Gdg))
    else:
      # Tare matrix
      Gtare = self.setupTareDesign(stations, tare)
      G = np.hstack((Gpoly.T, Gdg, Gtare.T))

    # Weight matrix
    W = np.diag(np.reciprocal(np.square(s)))

    # Invert to model parameters & uncertainties
    lsq, std, residuals, chi = self.__invert(G, W, y)

    # These are the drift parameters
    mbeta = lsq[:degree + 1]

    if tare is None:
      # Gravity differences & uncertainties
      mdg = lsq[degree + 1:]
      stddg = std[degree + 1:]
    else:
      # Gravity differences & uncertainties
      mdg = lsq[degree + 1:-1]
      stddg = std[degree + 1:-1]

    # Eliminate the gravity differences for plotting
    y -= Gdg @ mdg

    # Eliminate the tare
    dtare = None
    if tare is not None:
      y -= Gtare[0] * lsq[-1]
      dtare = lsq[-1]

    return InversionResult(self, degree, anchor, mbeta, x, y, s, mdg, stddg, residuals, changes, stations, chi, tare, dtare, Gtare[0])
