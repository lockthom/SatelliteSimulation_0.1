% This is to test orbital performance, rather than the attitude performance

close all
clear all
clc

function out_orbitDeriv = orbitDeriv(t, rv_Input)

    muEarth = 398600; % km^3/s^2 

    rVec = rv_Input(1:3);
    vVec = rv_Input(4:6);
    rVal = norm(rVec);
    vVal = norm(vVec);
    
    oParam = (muEarth./(rVal.^3));

    derivative_of_position = vVec;
    derivative_of_velocity = -oParam.*rVec;

    out_orbitDeriv = [derivative_of_position;derivative_of_velocity];


end

rInit = [ 7000; 0; 0]; % km
vInit = [ 0; -5; 0]; % km/s

tVals_orbit = linspace(0, 32000, 1000);

[t_orbit,y_orbit] = ode45(@orbitDeriv, tVals_orbit, [rInit; vInit]);

y_orbit = y_orbit';

figure(1)
hold on
grid on
plot3(y_orbit(1,:), y_orbit(2,:), y_orbit(3,:),'k')
title('Orbital Path')
% xlim([-2000,6000])
% ylim([-2000,6000])
% zlim([-1000,6000])
xlabel('X')
ylabel('Y')
zlabel('Z')
view(3)

% figure(2)
% hold on
% grid on
% plot(y_orbit(1,:), y_orbit(2,:))
% title('XY Orbital Path')
% xlabel('X')
% ylabel('Y')

% figure(3)
% hold on
% grid on
% plot(t_orbit,y_orbit(1,:))
% title('X Coordinate over Time')
% xlabel('time')
% ylabel('x coordinate')
