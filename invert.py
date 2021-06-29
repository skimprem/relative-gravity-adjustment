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

  data = DataLoader.load("CG5", "CG5.dat")
  result = data.invert(2)
  result.plot(removeDrift=True)
  result.plotResiduals()
  result.save("test.dat")
