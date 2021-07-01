import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates

from src.tide import TidalModel
from src.oceanloading import OceanLoadingModel
from src.inversionresult import InversionResult
from matplotlib import cm
from datetime import datetime, timedelta

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

    latitude = 19.40840
    longitude = -155.28385
    height = 1000

    model = TidalModel(latitude, longitude, height)

    x = self.df["Date_Time"]

    loadingModel = OceanLoadingModel("harmonics/hawaii.txt").getOceanLoadingModel(self.df["Date_Time"])

    eternaModel = model.getETERNA(x)

    y = model.getLongman(x)
    plt.plot(x, y, label="Longman")
    plt.plot(x, -eternaModel(mdates.date2num(x)), label="ETERNA 3.4")
    plt.plot(x, -eternaModel(mdates.date2num(x)) + loadingModel(mdates.date2num(x)), label="ETERNA 3.4 + Ocean Loading")
    plt.plot(x, 1E3 * self.df["TideCorr"], label="Default")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    plt.legend()
    plt.show()

  def correctLoading(self, y):

    x = self.df["Date_Time"]
    loadingModel = OceanLoadingModel("harmonics/hawaii.txt").getOceanLoadingModel(self.df["Date_Time"])

    # Ocean loading needs to be ADDED because it is a CORRECTION
    return y + loadingModel(mdates.date2num(x))


  def correctTide(self, y, which):

    x = self.df["Date_Time"]

    # Eliminate the default correction. What is given is the CORRECTION not the EFFECT so undo the correction by subtraction
    y -= 1E3 * self.df["TideCorr"]

    latitude = 19.40840
    longitude = -155.28385
    height = 1000

    model = TidalModel(latitude, longitude, height)

    if which == "ETERNA":
      return y - model.getETERNA(x)(mdates.date2num(x))
    elif which == "Longman":
      return y + model.getLongman(x)
    else:
      raise ValueError("Unknown tidal correction requested.")


  def invert(self, degree, anchor=None, tide="default", loading=False):

    """
    def DataWrapper.invert
    Main routine for the inversion
    """

    # Read data from file
    stations, x, y, s = self.extract()

    # Correct for the tide using ETERNA
    if tide != "default":
      y = self.correctTide(y, tide)

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
