function result = CalPhase(lp, x, y)
%lp: (theta, cycle, phase)
theta = lp(1);
cycle = lp(2);
phase = lp(3);

ret = (x.*cos(theta) + y.*sin(theta)) * cycle + phase;
result = checkPhase(ret);
end