#! /usr/bin/env python

import sys
import getopt
from primitives2d import genstair
from math import pi

optstr = "vhl:d:n:a:L:H:x:"
(options, args) = getopt.gnu_getopt(sys.argv[1:], optstr)

left = 0.
down = 0.
angle = 90
length = 0.3
height = 0.05
dens = 2e-3
num = 3

for (opt, optarg) in options:
    if opt == '-l':
        left = float(optarg)
    if opt == '-d':
        down = float(optarg)
    if opt == '-n':
        num = int(optarg)
    if opt == '-a':
        angle = float(optarg)
    if opt == '-H':
        height = float(optarg)
    if opt == '-L':
        length = float(optarg)
    if opt == '-x':
        dens = float(optarg)

angle = angle / 360.0 * 2 * pi

genstair(left, down, num, length, height, angle, dens)

