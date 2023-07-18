from src.dataloader import DataLoader

if __name__ == "__main__":

  """
  Example showing how to use CG5, CG6 data
  More complex analysis requires pygtide / ocl parameters
  """

  # data = DataLoader.load("CG5", "example/CG5.dat")
  # result = data.invert(2)
  # result.plot("CG5.pdf")

  # data = DataLoader.load("CG6", "example/LINE N 23 FROM CG-6_0368_9347-VG.dat .txt")
  data = DataLoader.load("CG6", "/home/roman/gitrepo/simple_grav_proc/example/vertical_gradient/CG-6_0461_VG_1253.dat")
  data.setLocations("example/locations/stations.csv")
  result = data.invert(degree=1)
  result.save("result.txt")
  result.plot('result.pdf')
