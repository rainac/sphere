function P = taitseq(rho, K0, m2, rho0)
  P = -K0 * (rho0 - rho) ./ ((rho0 - rho)*m2 + rho);
