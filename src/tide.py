from scipy.interpolate import interp1d
import numpy as np
import matplotlib.dates as mdates
from longman.longman import TideModel

class TidalModel():

  """
  Class TidalModel
  """

  def __init__(self, latitude, longitude, height):

    self.latitude = latitude
    self.longitude = longitude
    self.height = height


  def getLongman(self, times):
 
    """
    Returns Longman tidal model
    Taken from https://github.com/jrleeman/LongmanTide
    """
 
    model = TideModel()  # Make a model object
 
    solution = list()

    # Not vectorized
    for time in times:
      _, _, g = model.solve_longman(self.latitude, self.longitude, self.height, time)
      solution.append(g)

    # Scale results to microGal
    return 1E3 * np.array(solution)


  def getETERNA(self, times):
 
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
      self.latitude,
      self.longitude,
      self.height,
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
