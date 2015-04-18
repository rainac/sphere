#! /usr/bin/env python

import sys
import getopt
from primitives2d import gentriang

optstr = "vhl:r:u:d:x:c:H:I:"
(options, args) = getopt.gnu_getopt(sys.argv[1:], optstr)

left = 0.
crown = 0.5
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
    if opt == '-c':
        crown = float(optarg)
    if opt == '-x':
        dens = float(optarg)

gentriang(left, down, right, up, crown, dens);

