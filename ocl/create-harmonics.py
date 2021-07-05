import os

if __name__ == "__main__":

  """
  Script to rewrite ocean loading service response to harmonics per station
  """

  if not os.path.exists("harmonics.txt"):
    raise ValueError("Input file (harmonics.txt) from the ocean loading provider is missing.")

  if not os.path.exists("harmonics"):
    os.makedirs("harmonics")

  with open("results.txt", "r") as infile:

    counter = 0

    for line in infile.read().split("\n"):
  
      if line.startswith("$$"):
        continue
  
      if counter == 0:
        collection = ""
        name = line.strip()
      else:
        collection += line + "\n"
  
      counter += 1

      if counter == 7:
        counter = 0

        with open(os.path.join("harmonics", "%s.txt" % name), "w") as outfile:
          outfile.write(collection) 
