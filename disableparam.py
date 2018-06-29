#!/usr/bin/env python

from __future__ import print_function
from sys import argv



par=argv[1]
params=argv[2:]


with open(par) as f:
    for line in f:
        e=line.split()
        if len(e) > 0:
            if e[0] in params:
                print("{:<14} {}".format(e[0],e[1]))
            else:
                print(line.strip())


