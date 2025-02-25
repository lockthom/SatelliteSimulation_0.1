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


qr = [1;0;0;0]; % Standard part of quaternion
qd = [0;0;0;0]; % Dual part of quaternion


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

function out_quMult = quMult(qu_L, qu_R)
    % Inputs:
    % 1: qu_L -> Left quaternion
    % 2: qu_R -> RIght quaternion
    % 
    % Output:
    % 1: out_quMult -> Quaternion product

    qu_r = qu_L(1)*qu_R(1) - dot(qu_L(2:end),qu_R(2:end));
    qu_v = qu_L(1)*qu_R(2:end) + qu_R(1)*qu_L(2:end) + cross(qu_L(2:end),qu_R(2:end));

    out_quMult = [qu_r;qu_v];

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

        out_qMag = quMult(q,quConj(q));

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
    quDmult_d = quMult(qD_L(1:4),qD_R(5:end)) + quMult(qD_L(5:End),qD_R(1:4));

    out_quDMult = [quDmult_r;quDmult_d];

end

function out_quDconj = quDconj(qD)
    % Inputs:
    % 1: qD -> Dual quaternion
    %
    % Output:
    % 1: out_quDconj
   
    out_quDconj = quConj(qD(1:4)) + quConj(qD(5:end));

end

function out_dualUnity = dualUnity(qD)
    % Inputs:
    % 1: qD -> 8x1 dual quaternion
    % 
    % Outputs:
    % 1: out_dualUnity -> 0 if unit. May be small number due to rounding.

    out_dualUnity = quMult(quConj(qD(1:4)),qD(5:end)) + quMult(quConj(qD(5:end)),qD(1:4));
    
end


% Translation vector
translation_vector = [0;0;1];

% Rotation vector
rotation_vector = nth2quat([0;0;1],deg2rad(45));
rotation_vector = rotation_vector./quMag(rotation_vector);

% Form dual quaternion
quD_r = rotation_vector;
quD_d = 0.5.*(quMult(translation_vector, rotation_vector));
qD = [quD_r;quD_d];

