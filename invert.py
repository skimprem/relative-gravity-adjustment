from src.dataloader import DataLoader
import sys
import os

def solve(campaign, instrument, filename):

  print("Solving %s" % filename)

  filepath = os.path.join("data", campaign, instrument, filename)

  data = DataLoader.load("USGS", filepath)
  data.setLocations("locations/stations.csv")
  result = data.invert(1, tide="Longman", loading=False)
  result.plot(os.path.join("figures", campaign, "%s.pdf" % filename))

  # Sometimes a circuit is measured relative to HVO1: convert to P1 THROUGH its difference with HVO41.
  if result.anchor == "HVO41":

    if filename == "578_2012-06-22.csv":
      result.relativeTo("P1", 29891.0, 5.67)
    if filename == "578_2011-03-23.csv":
      result.relativeTo("P1", 29799.0, 2.61)
    if filename == "578_2012-11-27.csv":
      result.relativeTo("P1", 29947.0, 2.76)

    if filename == "579_2012-06-22.csv":
      result.relativeTo("P1", 29915.0, 9.1)
    if filename == "579_2011-03-23.csv":
      result.relativeTo("P1", 29811.0, 5.9)
    if filename == "579_2012-11-27.csv":
      result.relativeTo("P1", 29965.0, 2.87)

  result.save("results/%s/%s/%s.dat" % (campaign, instrument, filename))


if __name__ == "__main__":

  """
  Examples
  """

  for campaign in os.listdir("data"):
    for instrument in os.listdir(os.path.join("data", campaign)):
      for filename in os.listdir(os.path.join("data", campaign, instrument)):
        if filename != "579_2012-10-24.csv":
          continue
        solve(campaign, instrument, filename)
