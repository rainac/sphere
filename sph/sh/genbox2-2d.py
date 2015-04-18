#! /usr/bin/env python

import sys


def genBox(dens, c_ul, c_or):
    N = int((c_or[0] - c_ul[0]) / dens);
    
    lx = c_or[0] - c_ul[0];
    ly = c_or[1] - c_ul[1];
    
    dx = lx/N
    dy = ly/N
    
    for i in range(N):
        print c_ul[0], c_ul[1] + i*dx, c_ul[2]
        print c_or[0], c_ul[1] + i*dx, c_ul[2]

    for i in range(N):
        print c_ul[0] + i*dx, c_ul[1], c_ul[2]

if __name__ = "__main__":
    dens = 1e-3;
    if len(sys.argv) > 1:
        dens = float(sys.argv[1])
    c_ul = [0.29,   0.19, .9];
    c_or = [0.71,   0.51, .9];
    genBox(dens, c_ul, c_or)
