import pandas as pd
import numpy as np

from src.datawrapper import DataWrapper

class DataLoader():

  """
  Class DataLoader
  Loads CG5 or CG6 data from disk
  """

  def __init__(self):
    pass


  @classmethod
  def load(self, ver, filepath):

    """
    def DataLoader.load
    Reads data from a specific format
    """

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

    from datetime import timedelta

    df = pd.read_csv(filepath, comment="/", delimiter="\t", parse_dates=[[1, 2]])

    return DataWrapper("CG6", df, filepath)


  @classmethod
  def readUSGS(self, filepath):

    """
    def DataLoader.readUSGS
    Reads USGS gravity data
    """

    df = pd.read_csv(filepath, delimiter="\t", parse_dates=[[15, 12]])
    header = ["Date_Time", "Station", "Latitude", "Longitude", "Altitude", "CorrGrav", "StdErr", "TiltX", "TiltY", "Temp", "TideCorr", "duration", "rej", "Dec. Time+Date", "Terrain", "Accepted"]

    df.columns = header
    df["StdErr"] = df["StdErr"] / np.sqrt((df["duration"] - df["rej"]))

    return DataWrapper("USGS", df, filepath)


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

    return DataWrapper("CG5", df, filepath)
