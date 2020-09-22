Kp = 100;
time = 0:0.001:.1;
num = Kp;
den = [0, 1, (0.5+Kp)];
sys = tf(num,den);
xstep = step(sys,time);
figure(1)
plot(time,xstep);
u_t = -Kp*((Kp/(0.5+Kp)).*(1-exp(-(0.5+Kp).*time))-1);
figure(2)
plot(time,u_t)


Kp_2 = 50;
time_2 = 0:0.001:.1;
num = Kp_2;
den = [0, 1, (0.5+Kp_2)];
sys = tf(num,den);
xstep_2 = step(sys,time_2);
figure(3)
plot(time_2,xstep);
u_t_2 = -Kp_2*((Kp_2/(0.5+Kp_2)).*(1-exp(-(0.5+Kp_2).*time_2))-1);
figure(4)
plot(time,u_t_2)

