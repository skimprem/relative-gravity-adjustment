from src.dataloader import DataLoader

if __name__ == "__main__":

  """
  Example showing how to use CG5, CG6 data
  More complex analysis requires pygtide / ocl parameters
  """

  data = DataLoader.load("CG5", "example/CG5.dat")
  result = data.invert(2)
  result.plot("CG5.pdf")

  data = DataLoader.load("CG6", "example/CG6.dat")
  result = data.invert(2, anchor="40")
  result.plot("CG6.pdf")
