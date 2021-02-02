#!/usr/bin/env python


import sys

par=sys.argv[1]
tim=sys.argv[2]

parlines=[]
glitches={}
with open(par) as f:
    for line in f:
        if line.startswith("GL"):
            e=line.split()
            pp = e[0].split("_")
            i = int(pp[1])
            param=pp[0]
            if not i in glitches:
                glitches[i] = {}
            if param=="GLEP":
                glitches[i]['epoch'] = float(e[1])
            if param=="GLPH":
                glitches[i]['turns'] = round(float(e[1]))
                e[1] = "{}".format(float(e[1]) - glitches[i]['turns'])

            glitches[i][param] = " ".join(e[1:])

        else:
            parlines.append(line)


for ig in glitches:
    print("glitch[{}] epoch {} turns {}".format(ig,glitches[ig]['epoch'],glitches[ig]['turns']))


gg = sorted(glitches,key=lambda x: glitches[x]['epoch'])
print("Time-order is",gg)

i=1
for ig in gg:
    glitches[ig]['newid'] = i
    i+=1


timlines=[]
with open(tim) as f:
    for line in f:
        e=line.split()
        if "-pn" in e:
            epoch = float(e[2])
            ii = e.index("-pn")
            pn = int(e[ii+1])
            for ig in gg:
                if epoch > glitches[ig]['epoch']:
                    pn -= glitches[ig]['turns']
            newline = " ".join(e[:ii])+" -pn {} ".format(pn)+(" ".join(e[ii+2:]))
            timlines.append(" "+newline+"\n")
        else:
            timlines.append(line)


with open(tim+".fix","w") as f:
    f.writelines(timlines)

with open(par+".fix","w") as f:
    f.writelines(parlines)
    for ig in gg:
        for param in glitches[ig]:
            if param in ["epoch","turns","newid"]:
                continue
            f.write("{}_{} {}\n".format(param,glitches[ig]['newid'],glitches[ig][param]))



