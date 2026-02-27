import sys

# load filenames
if (len(sys.argv) != 2):
  print("The script needs a PDB filename with solutes")
  exit(0)
else:
  filename = sys.argv[1]

# list of solute residue names included in this script
solutes = "sol1 sol2"

# parameters for extraBonds walls
# (weak repulsion for x<xlower)
k=0.1
lower = 6.0
upper=1000.0

fo = open(filename, "r")
lines = fo.readlines()
fo.close()

list = []
listrsn = []
for line in lines:
  if "ATOM" in line:
    #sline = line.split() 
    atomn=line[12:17].strip()
    resn=line[17:21].strip()
    if(resn in solutes) and (atomn[0] != "H"):
      #list.append(int(sline[1]))
      #listrsn.append(int(sline[4]))
      list.append(int(line[6:11]))
      listrsn.append(int(line[22:26]))

for i in range(len(list)-1):
  for j in range(i+1, len(list)):
    if(listrsn[i] != listrsn[j]):
      print("wall ", list[i]-1, list[j]-1, k, lower, upper)
