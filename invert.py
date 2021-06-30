from lib import DataLoader
import sys
import os

if __name__ == "__main__":

  """
  Examples
  """

  for file in os.listdir("."):
    if not file.endswith(".csv"): continue
    data = DataLoader.load("USGS", file)
    result = data.invert(1, tide="ETERNA")
    result.plot(removeDrift=True)
    result = data.invert(1, tide="Longman")
    result.plot(removeDrift=True)
    result = data.invert(1)
    result.plot(removeDrift=True)

  sys.exit(0)
  data = DataLoader.load("CG6", "CG6.dat")
  result = data.invert(1, tide="ETERNA")
  result.plot(removeDrift=False)
  result = data.invert(1, tide="Longman")
  result.plot(removeDrift=False)
  result = data.invert(1)
  result.plot(removeDrift=False)
  result.plotResiduals()
  result.save("CG6-results.dat")

  data = DataLoader.load("CG5", "CG5.dat")
  result = data.invert(1, tide="ETERNA")
  result.plot(removeDrift=False)
  result = data.invert(1, tide="Longman")
  result.plot(removeDrift=False)
  result = data.invert(1)
  result.plot(removeDrift=False)
  result.plotResiduals()
  result.save("CG5-results.dat")
