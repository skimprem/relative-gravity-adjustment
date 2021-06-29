from lib import DataLoader

if __name__ == "__main__":

  """
  Examples
  """

  df = DataLoader.load("CG6", "CG6.dat")

  inv = df.invert(1)
  inv.plot()
  inv.plotResiduals()

  df = DataLoader.load("CG5", "CG5.dat")
  inv = df.invert(2)
  inv.plot()
  inv.plotResiduals()
