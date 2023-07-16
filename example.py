from src.dataloader import DataLoader

if __name__ == "__main__":

  """
  Example showing how to use CG5, CG6 data
  More complex analysis requires pygtide / ocl parameters
  """

  data = DataLoader.load("CG6", "/home/roman/gitrepo/simple_grav_proc/example/vertical_gradient/CG-6_0461_VG_1253.dat")
  data.setLocations("example/locations/stations.csv")
  result = data.invert(degree=1, anchor="1", tide="Longman")
  result.plot("CG-6_0461_VG_1253.pdf")
  ##result.plotResiduals()
  #result.differences
  result.save("CG-6_0461_VG_1253.txt")

  #data = DataLoader.load("CG6", "example/CG-6_0461_VG_1253.dat")
  #data.setLocations("example/locations/stations.csv")
  #result = data.invert(degree=6, anchor="1", tide="Longman")
  #result.plot("CG-6_0461_VG_1253.pdf")
  ##result.plotResiduals()
  #result.differences
  #result.save("CG-6_0461_VG_1253.txt")

  #data = DataLoader.load("CG5", "example/CG5.dat")
  #result = data.invert(1)
  #result.plot("CG5.pdf")
  #result.plotResiduals()
  #result.differences
  #result.save("CG5.txt")

  #data = DataLoader.load("CG6", "example/CG6.dat")
  #result = data.invert(1, anchor="40")
  #result.plot("CG6.pdf")
  #result.plotResiduals()
  #result.differences
  #result.save("CG6.txt")
