from BioAid import MMBSearchTK as MMTK

import sys

#PARAMS:

filter = sys.argv[3]
if filter == "False":
  f_name = "final_bir_locs.txt"
elif filter == "True":
  f_name = "non_commons_fromAll.txt"
else:
  print("argument not recognized. Exitting")
  exit()


path=sys.argv[1]
#End PARAMS


MMTK.createAnnotatedOutput(f_name, path, sys.argv[2])
