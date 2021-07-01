from subprocess import Popen, PIPE, DEVNULL
from scipy.interpolate import interp1d
from io import StringIO

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

class OceanLoadingModel():

  """
  Class OceanLoadingModel
  Container for ocean loading models. Harmonics should be supplied from the ocean tide loading provider
  that can be accessed at: http://holt.oso.chalmers.se/loading/
  """

  def __init__(self, harmonics):

    """
    def OceanLoadingModel.__init__
    Initializes the class: harmonics should point to a file that contains the site-specific harmonics from the ocean tide provider.
    The contents of the file should look similar to:

       33.97  12.84   8.19   4.25  36.51  20.31  10.84   3.19   1.66   0.61   3.05
       15.63   8.04   3.64   2.06   8.77   5.36   2.64   1.00   0.08   0.10   0.56
       33.12  10.33   6.47   3.24  12.88   6.68   3.90   1.12   1.35   0.58   0.76
      -116.8 -135.1 -140.8 -136.5   52.5   43.1   51.7   38.3 -139.5  175.9   99.2
        -7.1    5.0   -9.3   16.9 -133.2 -152.2 -134.4 -166.7   34.5   47.4 -144.5
      -101.9  -83.7 -126.5  -94.3   15.0   -2.5   11.5  -28.7 -168.5 -171.0  173.6

    For gravity make sure the input is given in nm/s^2.
    """

    self.harmonics = harmonics


  def callHARDISP(self, cmd):
  
    """
    def callHARDISP
    Subprocess call to run HARDISP on input using harmonics from filepath
    """
  
    # Pop the standard input of the harmonics file to HARDISP
    stdin = open(self.harmonics, "r")
    p = Popen(cmd, stdin=stdin, stdout=PIPE, stderr=DEVNULL)
    stdout = p.communicate()
  
    # First column of output (nm/s^2) to microGal (0.1)
    return 1E-1 * np.loadtxt(StringIO(stdout[0].decode("ascii"))).T[0]
  
  
  def getOceanLoadingModel(self, times):
  
    """
    get getOceanLoadingModel
    Returns the Ocean Loading Model
    """
  
    HARDISP_DIR = "./hardisp/HARDISP"
  
    # Get the start and end of the segment
    start = times.iloc[0]
    end = times.iloc[-1]
  
    # Determine the number of minutes between the start and ending with some padding
    minutes = int((((end - start).total_seconds()) / 60) + 10)
  
    # HARDISP input is YYYY MM DD HH MM SS NUM SAM
    # Where NUM is the number of samples and SAM is the sampling interval
    startParameters = start.strftime("%Y %m %d %H %M %S").split(" ")
  
    # This is the CMD we pass to HARDISP
    cmd = [HARDISP_DIR] + startParameters + [str(minutes), "60"]
  
    # Add the numbers (one unit in mdates is a full day: 1440 minutes in a day)
    x = mdates.date2num(start) + (np.arange(minutes) / 1440)
  
    # Call HARDISP
    y = self.callHARDISP(cmd)
  
    # Return the an interpolated model
    return interp1d(x, y, kind="cubic")
