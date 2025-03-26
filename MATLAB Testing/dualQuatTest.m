% This file is intended to make screw-like motion along one axis using dual
% quaternions. It is a conceptual test run to make sure that I
% understand how to use dual quaternions in code before I start using them
% in more complex versions of the system.

% R   -> Right
% L   -> Left
% r   -> Real part (standard)
% d   -> Dual part
% D   -> Dual quaternion (total)
% n/N -> Dual number, not dual quaternion
% v   -> Vector part

close all
clear all
clc

qr = [1;0;0;0]; % Standard part of quaternion
qd = [0;0;0;0]; % Dual part of quaternion

function out_quMult = quMult(q0,q1)


    qO1 = q0(1).*q1(1) - dot(q0(2:end),q1(2:end));

    qOV = q1(1).*q0(2:end) + q0(1).*q1(2:end) + cross(q0(2:end),q1(2:end));

    out_quMult = [qO1;qOV];

end

function out_dualAddN = dualAddN(dn_L,dn_R)
    % Inputs:
    % 1: qDn_L -> Dual number, left side
    % 2: qDn_R -> Dual number, right side
    %
    % Outputs: 
    % 1: out_dualAddN -> Added dual number
    
    dualAdd_rn = dn_L(1) + dn_R(2);
    dualAdd_dn = dn_L(1) + dn_R(2);

    out_dualAddN = [dualAdd_rn;dualAdd_dn];

end

function out_dualMultN = dualMultN(dn_L,dn_R)
    % Inputs: 
    % 1: qDn_L -> Dual number, left side
    % 2: qDn_R -> Dual number, right side
    %
    % Outputs:
    % 1: out_dualMultN -> Multiplied dual number
    
    dualMultN_rn = dn_L(1)*dn_R(1);
    dualMultN_dn = dn_L(1)*dn_R(2) + dn_L(2)*dn_R(1);

    out_dualMultN = [dualMultN_rn;dualMultN_dn];

end

function out_dualDivN = dualDivN(dn_num, dn_den)


    % Inputs:
    % 1: dn_num -> Numerator in dual number division
    % 2: dn_den -> Denominator in dual number division
    % 
    % Output:
    % 1: out_dualDivN -> Divided dual number

    if dn_den(2) == 0

        error("Zero dual part in denominator, division by zero!")

    end

    dualDivN_rn = (dn_num(1)*dn_den(1))/(dn_den(1)^2); % Cancellation?
    dualDivN_dn = (dn_den(1)*dn_num(2) - dn_num(1)*dn_den(2)/(dn_den(1)^2));

    out_dualDivN =[dualDivN_rn;dualDivN_dn];


end

function out_qConj = quConj(q)
    % Inputs:
    % 1: q -> Quaternion
    %
    % Output:
    % 1: out_qConj -> Quaternion conjugate

    out_qConj = [q(1);-q(2);-q(3);-q(4)];


end

function out_qMag = quMag(q)
    % Input:
    % 1: q -> Quaternion
    %
    % Output:
    % 1: 2-Norm Magnitude of q

    if (length(q) == 4) || (length(q) == 8)

        out_qMag = norm(quMult(q,quConj(q)));

    else

        error("Non-quaternion/dual quaternion entry. Size is wrong")
    
    end

end

function out_q = nth2quat(n,th)
    % Inputs: 
    % 1: n  -> 3x1 vector axis of rotation
    % 2: th -> Rotation angle in radians about n
    %
    % Outputs:
    % 1: out_q -> Equivalent quaternion

    q_r = cos(th/2);
    q_v = sin(th/2).*n;

    out_q = [q_r;q_v];


end

function out_quDMult = quDmult(qD_L,qD_R)
    % Inputs:
    % 1: qD_L -> Left dual quaternion
    % 2: qD_R -> Right dual quaternion
    %
    % Output:
    % 1: out_quDMult -> Dual quaternion product

    % q1*q2 = qr1*qr2 + (qr1*qd2 + qd1*qr2)epsilon;
    quDmult_r = quMult(qD_L(1:4),qD_R(1:4));
    quDmult_d = quMult(qD_L(1:4),qD_R(5:end)) + quMult(qD_L(5:end),qD_R(1:4));

    out_quDMult = [quDmult_r;quDmult_d];

end

function out_quDconj = quDconj(qD)
    % Inputs:
    % 1: qD -> Dual quaternion
    %
    % Output:
    % 1: out_quDconj
   
    out_quDconj = [quConj(qD(1:4));quConj(qD(5:end))];

end

function out_dualUnity = dualUnity(qD)
    % Inputs:
    % 1: qD -> 8x1 dual quaternion
    % 
    % Outputs:
    % 1: out_dualUnity -> 0 if unit. May be small number due to rounding.
    
    unityCheck = quMult(quConj(qD(1:4)),qD(5:end)) + quMult(quConj(qD(5:end)),qD(1:4));

    if unityCheck <= 0.001

        qD = qD./norm(qD);
        out_dualUnity = qD;

    elseif unityCheck > 0.001

        error('Dual Quaternion out of tolerance')

    end
    
end

%% Define parameters

translation_vector = [0;0;0];
axis_of_rotation = [0;0;1];
angle_of_rotation = deg2rad(0);

starting_point = [1;0;0];
starting_quaternion = [1;0;0;0];

% Translation vector (move along z axis by 0 units)
% Rotation vector (rotate about z axis 0 degrees)

translation_quaternion = [0;translation_vector];
rotation_quaternion = nth2quat(axis_of_rotation, angle_of_rotation);
rotation_quaternion = rotation_quaternion./quMag(rotation_quaternion);

% Form initial dual quaternion
quD_r = rotation_quaternion;
quD_d = 0.5.*(quMult(translation_quaternion, rotation_quaternion));
qD = [quD_r;quD_d];

% Pull out translation component to ensure consistency.
translation_from_dual_quaternion = 2.*(quMult(quD_d,quConj(quD_r)));
error_in_translation = abs(translation_quaternion - translation_from_dual_quaternion);


%% Define dual quaternion initial derivative

% Rotation about z-axis at 0.25 rad/s
angular_velocity_init = [0;0;0.25];
linear_velocity_init = [0;0;0];

dualVel = [[0;angular_velocity_init];[0;linear_velocity_init]];

qD_dot_init = 0.5.*quDmult(qD,dualVel);

function dynamicOutput = dynamicDeriv(~,y)

    dQ = y;
    
    % Rotation about z-axis at 0.25 rad/s
    angular_velocity_init = [0;0;0.25];
    linear_velocity_init = [0;0;1];
    
    dualVel = [[0;angular_velocity_init];[0;linear_velocity_init]];

    dynamicOutput = 0.5.*quDmult(dQ,dualVel);

end

[timeVals,outVals] = ode45(@dynamicDeriv,[0,60],qD);


% Now make these values into angles and translations that are more easily
% readable to make sure that the results here are accurate as far as is
% relevant for the application here.

lengthIteration = max(size(timeVals));
eulerAngles = zeros(lengthIteration,3);
t_dQ = zeros(lengthIteration,3);

for ii = 1:lengthIteration

    outVals = outVals';

    eulerAngles(ii,:) = quat2eul(outVals(1:4,ii)');
    tempVal = 2.*(quMult(outVals(5:end,ii), outVals(1:4,ii)));
    t_dQ(ii,:) = tempVal(2:end);

end

A = 1;

% figure(1)
% hold on
% grid on
% title("Dual Quaternion Dynamics: Rotation about Z")
% plot(timeVals,outVals(:,1))
% xlabel("Time (Seconds)")
% ylabel("Scalar Pure Quaternion Component")

% figure(2)
% hold on
% grid on
% title("Dual Quaternion Dynamics: Rotation about Z")
% plot(timeVals,outVals(:,2))
% xlabel("Time (Seconds)")
% ylabel("Pure Quaternion First Component")
% 
% figure(3)
% hold on
% grid on
% title("Dual Quaternion Dynamics: Rotation about Z")
% plot(timeVals,outVals(:,3))
% xlabel("Time (Seconds)")
% ylabel("Pure Quaternion Second Component")

% figure(4)
% hold on
% grid on
% title("Dual Quaternion Dynamics: Rotation about Z")
% plot(timeVals,outVals(:,4))
% xlabel("Time (Seconds)")
% ylabel("Pure Quaternion Third Component")

figure(5)
hold on
grid on
title("Dual Quaternion Dynamics: Rotation about Z")
plot(timeVals,outVals(:,5))
xlabel("Time (Seconds)")
ylabel("Dual Quaternion Scalar Component")

figure(6)
hold on
grid on
title("Dual Quaternion Dynamics: Rotation about Z")
plot(timeVals,outVals(:,6))
xlabel("Time (Seconds)")
ylabel("Dual Quaternion First Component")

figure(7)
hold on
grid on
title("Dual Quaternion Dynamics: Rotation about Z")
plot(timeVals,outVals(:,7))
xlabel("Time (Seconds)")
ylabel("Dual Quaternion Second Component")

figure(8)
hold on
grid on
title("Dual Quaternion Dynamics: Rotation about Z")
plot(timeVals,outVals(:,8))
xlabel("Time (Seconds)")
ylabel("Dual Quaternion Third Component")

