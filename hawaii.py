from src.dataloader import DataLoader
import sys
import os

def getTare(filename):

  """
  def getTare
  This helps the inversion to fix the tare.
  """

  # The returned value is the offset of the sample
  # These are definitely tares
  if filename == "579_2009-12-02.csv":
    return 20
  if filename == "578_2012-11-27.csv":
    return 72
  if filename == "579_2009-12-15.csv":
    return 34
  if filename == "578_2017-04-21.csv":
    return 66

  # Could be but.. it is really difficult to say for sure
  if filename == "578_2017-04-24.csv":
    return 51
  if filename == "579_2017-04-24.csv":
    return 60
  if filename == "578_2017-04-20.csv":
    return 40
  if filename == "579_2015-09-15.csv":
    return 56

  return None


def solve(campaign, instrument, filename):

  """
  def solve
  Calls relative gravity code and writes results to file
  """

  print("Solving %s" % filename)

  filepath = os.path.join("data", campaign, instrument, filename)

  # Load the data
  data = DataLoader.load("USGS", filepath)
  data.setLocations("locations/stations.csv")

  # Get index of when tare needs to be restored
  tare = getTare(filename)

  # Complete the inversion
  #result = data.invert(1, tide="ETERNA", loading=True, tare=tare)
  result = data.invert(1, tide="Longman", loading=False, tare=tare)

  # Plot the inversion results
  result.plot(os.path.join("figures", campaign, "%s.pdf" % filename), removeDrift=True)

  # Save the results
  result.save("results/%s/%s/%s.dat" % (campaign, instrument, filename.split(".")[0]))


if __name__ == "__main__":

  """
  Main entrypoint for relative gravity solutions
  """

  # Go over the gravity data
  for campaign in os.listdir("data"):
    for instrument in os.listdir(os.path.join("data", campaign)):
      for filename in os.listdir(os.path.join("data", campaign, instrument)):
        solve(campaign, instrument, filename)
