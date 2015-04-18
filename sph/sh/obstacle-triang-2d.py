#! /usr/bin/env python

import sys
import getopt
from math import sin, cos, sqrt, asin, acos

optstr = "vhl:r:u:d:x:H:I:"
(options, args) = getopt.gnu_getopt(sys.argv[1:], optstr)

left = 0.
right = 1.
down = 0.
up = 1.
dens = 2e-3

for (opt, optarg) in options:
    if opt == '-r':
        right = float(optarg)
    if opt == '-l':
        left = float(optarg)
    if opt == '-d':
        down = float(optarg)
    if opt == '-u':
        up = float(optarg)
    if opt == '-x':
        dens = float(optarg)

height = up - down
width = right - left
w_2 = width / 2

wallLength = sqrt(height*height + w_2*w_2);
slope = asin(height/wallLength);

N = int(wallLength / dens);

dx = wallLength/N

for i in range(N):
    print left + i*cos(slope)*dx, down + i*sin(slope)*dx, 0
    print left + w_2 + (N - i - 1)*cos(slope)*dx, down + i*sin(slope)*dx, 0
