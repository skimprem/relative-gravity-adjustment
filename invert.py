from lib import DataLoader

if __name__ == "__main__":

  """
  Examples
  """

  data = DataLoader.load("CG6", "CG6.dat")

  result = data.invert(1)
  result.plot()
  result.plotResiduals()

  data = DataLoader.load("CG5", "CG5.dat")
  result = data.invert(2)
  result.plot()
  result.plotResiduals()
