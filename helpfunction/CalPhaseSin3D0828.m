function result = CalPhaseSin3D0828(lp, x, y, z)
%lp: (theta, cycle, phase)
theta = lp(1);
cycle = lp(2);
thetaz = lp(3);
phase = lp(4);

ret = (x.*cos(theta).*cos(thetaz) + y.*sin(theta).*cos(thetaz)+z.*sin(thetaz)) * cycle + phase;
result = [sin(ret) cos(ret)];


end