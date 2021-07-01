from src.dataloader import DataLoader
import sys
import os

if __name__ == "__main__":

  """
  Examples
  """

  latitude = 52.0
  longitude = 4.38
  height = 0

  data = DataLoader.load("USGS", "578_2009-12-02.csv")
  #data = DataLoader.load("CG6", "CG6.dat")
  data.setLocations("locations/stations.csv")

  #result = data.invert(1)
  #result.plot(removeDrift=False)
  #result = data.invert(1, tide="Longman")
  #result.plot(removeDrift=False)
  #result = data.invert(1, tide="ETERNA")
  #result.plot(removeDrift=False)
  result = data.invert(1, tide="ETERNA", loading=True)
  result.plot(removeDrift=False)
  result.plotResiduals()
  result.save("CG6-results.dat")

  #data = DataLoader.load("CG5", "CG5.dat")
  #data.setLocation(latitude, longitude, height)
  #result = data.invert(1)
  #result.plot(removeDrift=False)
  #result = data.invert(1, tide="Longman")
  #result.plot(removeDrift=False)
  #result.save("CG5-results.dat")
