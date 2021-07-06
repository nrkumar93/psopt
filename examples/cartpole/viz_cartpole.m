clear
addpath /home/gaussian/Documents/softwares/CGPOPS_Release/examples/cartpole/cmake-build-debug
addpath /home/gaussian/cmu_ri_phd/phd_misc/psopt/cmake-build-debug/examples/cartpole

% run('cgpopsIPOPTSolutionCD1.m')

t = load('t.dat');
x = load('x.dat');
u = load('u.dat');

z = x;
z(1,:) = x(1,:);
z(2,:) = x(3,:);
z(3,:) = x(2,:);
z(4,:) = x(4,:);


p.m1 = 2.0;  % (kg) Cart mass
p.m2 = 0.5;  % (kg) pole mass
p.g = 9.81;  % (m/s^2) gravity
p.l = 0.5;   % (m) pendulum (pole) length

%%%% Draw Trajectory:
[p1,p2] = cartPoleKinematics(z,p);

figure(3); clf;
nFrame = 40;  %Number of frames to draw
drawCartPoleTraj(t,p1,p2,nFrame);

rmpath /home/gaussian/Documents/softwares/CGPOPS_Release/examples/cartpole/cmake-build-debug
rmpath /home/gaussian/cmu_ri_phd/phd_misc/psopt/cmake-build-debug/examples/cartpole
