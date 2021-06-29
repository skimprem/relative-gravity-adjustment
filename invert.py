from lib import DataLoader
import sys
if __name__ == "__main__":

  """
  Examples
  """

  data = DataLoader.load("CG6", "CG6.dat")

  result = data.invert(1)
  result.plot(removeDrift=False)
  result.plotResiduals()
  result.save("CG6-results.dat")

  data = DataLoader.load("CG5", "CG5.dat")
  result = data.invert(2, anchor="10.0")
  result.plot(removeDrift=False)
  result.plotResiduals()
  result.save("CG5-results.dat")
