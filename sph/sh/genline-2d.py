#! /usr/bin/env python

import sys
import getopt
from math import sqrt, ceil;

optstr = "vhl:r:u:d:x:H:I:"
(options, args) = getopt.gnu_getopt(sys.argv[1:], optstr)

left = 0.
right = 1.
down = 0.
up = 1.
holex=0.5
holew=0
dens = 2e-3

for (opt, optarg) in options:
    if opt == '-r':
        right = float(optarg)
    if opt == '-l':
        left = float(optarg)
    if opt == '-u':
        up = float(optarg)
    if opt == '-d':
        down = float(optarg)
    if opt == '-H':
        holex = float(optarg)
    if opt == '-I':
        holew = float(optarg)
    if opt == '-x':
        dens = float(optarg)

c_ul = [left, down, .0];
c_or = [right, up, .0];

dv = (c_or[0] - c_ul[0], c_or[1] - c_ul[1], 0,0);
length = sqrt(dv[0]**2 + dv[1]**2);

nv = (dv[0]/length, dv[1]/length, 0,0);

N = int(ceil(length / dens));
ds=length/N;
for i in range(N + 1):
    print c_ul[0] + nv[0]*i*ds, c_ul[1] + nv[1]*i*ds, 0

