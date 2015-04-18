% -*- octave -*-

h = 1e-23;
delta = h / 100;

X = 0:delta:2*h;

factor = pi * sqrt(pi) * h^3;
alpha = 1 / factor;

WX = 4 * pi * alpha * exp(- ((X/h) .* (X/h)) ) .* X .* X;

WX = WX * delta;

sum(WX)
