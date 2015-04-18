
rholow = 980;
rhohigh = 1020;

rhos = rholow:0.2:rhohigh;

B = 5e5;
K0 = 3.4e6;
m2 = 0;
g = 7;
rho0 = 1000;

Ps1 = gammaeq(rhos, B, g, rho0);
Ps2 = taitseq(rhos, K0, m2, rho0);

figure;
hold on;
plot(rhos, Ps1, 'r-');
plot(rhos, Ps2, 'b.');
legend('Gamma equation', 'Tait''s equation');
xlabel('Density $\rho$', 'interpreter', 'latex');
ylabel('Pressure $P(\rho)$', 'interpreter', 'latex');
title('Equations of state');
print('-depsc', 'img/eqs.eps');

