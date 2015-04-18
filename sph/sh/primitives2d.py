from math import sqrt, ceil, pi, sin, cos

# draw line from point (l,u) to (r, o)
def genline(l, u, r, o, ds):
   dv = (r - l, o - u);
   length = sqrt(dv[0]**2 + dv[1]**2);
   
   nv = (dv[0]/length, dv[1]/length);

   N = int(ceil(length / ds));
   ds=length/N;
   for i in range(N + 1):
       print l + nv[0]*i*ds, u + nv[1]*i*ds, 0
       

# draw triangle inscribed in box (l,u), (r, o)
# the crown is at x-location cx
def gentriang(l, u, r, o, cx, ds):
   genline(l, u, cx, o, ds)
   genline(cx, o, r, u, ds)

# draw triangle inscribed in box (l,u), (r, o)
# the crown is at x-location cx
def genstep(l, u, r, o, cx, ds):
   genline(l, u, cx, o, ds)
   genline(cx, o, r, u, ds)


# draw line from point (l,u) to (r, o)
def genlineFrom(x, y, len, phi, ds):
   N = int(ceil(len / ds))
   ds = len/N
   for i in range(N + 1):
       print x + sin(phi)*i*ds, y + cos(phi)*i*ds, 0

# draw line from point (l,u) to (r, o)
def genStep(x, y, len, height, phi, ds):
   genlineFrom(x, y, len, phi, ds);
   xend = x + sin(phi)*len
   yend = y + cos(phi)*len
   genlineFrom(xend, yend, height, phi+pi/2, ds);

# draw line from point (l,u) to (r, o)
def genStepAlt(x, y, len, height, phi, ds):
   genlineFrom(x, y, len, phi, ds);
   xend = x + sin(phi)*len
   yend = y + cos(phi)*len
   genlineFrom(xend, yend, height, phi-pi/2, ds);
   genlineFrom(xend, yend, height, phi+pi/2, ds);

def genstair(x, y, numSteps, len, height, phi, ds):
   xoff = x
   yoff = y
   for i in range(numSteps):
      if i % 2 == 0:
         genStep(xoff, yoff, len, height, phi, ds);
      else:
         genStepAlt(xoff, yoff, len, height, phi, ds);
      xoff = xoff + sin(phi)*(len)
      yoff = yoff + cos(phi)*(len)
      xoff = xoff + sin(phi+pi/2)*(height)
      yoff = yoff + cos(phi+pi/2)*(height)
