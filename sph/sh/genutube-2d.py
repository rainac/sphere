#! /usr/bin/env python

import sys
import getopt
from math import *

optstr = "vhc:l:r:L:R:d:x:"
(options, args) = getopt.gnu_getopt(sys.argv[1:], optstr)

# values must be given as floats
circ = 0.1 # width of pipe
ltop = 0.6 # height of left side
rtop = 0.6 # height of right side
down = 0.1 # bottom
left = 0.6 # left side border
right= 0.6 # right side border
dens = 2e-3

for (opt, optarg) in options:
    if opt == '-c':
        circ = float(optarg)
    if opt == '-L':
        ltop = float(optarg)
    if opt == '-R':
        rtop = float(optarg)
    if opt == '-l':
        left = float(optarg)
    if opt == '-r':
        right = float(optarg)
    if opt == '-d':
        down = float(optarg)
    if opt == '-x':
        dens = float(optarg)

l1 = (ltop - down)
l2 = (ltop - down - circ)
nl1 = int(ceil(l1 / dens))
nl2 = int(ceil(l2 / dens))

r1 = (rtop - down)
r2 = (rtop - down - circ)
nr1 = int(ceil(r1 / dens))
nr2 = int(ceil(r2 / dens))

no = int(ceil((right - left - 2*circ) / dens))
nd = int(ceil((right - left) / dens))

dl1 = l1 / nl1
dl2 = l2 / nl2

dr1 = r1 / nr1
dr2 = r2 / nr2

for i in range(nl1):
    print left, down + i*dl1, 0
for i in range(nl2):
    print left + circ, down + circ + i*dl2, 0

for i in range(nr1):
    print right, down + i*dr1, 0
for i in range(nr2):
    print right - circ, down + circ + i*dr2, 0

for i in range(no):
    print left + circ + i*dens, down + circ, 0
for i in range(nd):
    print left + i*dens, down, 0

