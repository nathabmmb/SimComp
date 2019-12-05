import re
source = open("BH_output.txt", "r")
txt = source.read()
nots = re.findall("not", txt)
output = open("nots.txt", "a")
output.write("{}\n".format(len(nots)))
output.close()
exit()