#! /usr/bin/env python

import sys
import getopt

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

Nx = int((c_or[0] - c_ul[0]) / dens);
Ny = int((c_or[1] - c_ul[1]) / dens);

lx = c_or[0] - c_ul[0];
ly = c_or[1] - c_ul[1];

dx = lx/Nx
if Ny != 0:
    dy = ly/Ny
else:
    dy = 1

for i in range(Ny + 1):
    print c_ul[0], c_ul[1] + i*dy, c_ul[2]
    print c_or[0], c_ul[1] + i*dy, c_ul[2]

for i in range(Nx - 1):
    x = c_ul[0] + (i+1)*dx
    if x <= holex or x >= holex + holew:
        print x, c_ul[1], c_ul[2]

#for i in range(Ny + 1):

