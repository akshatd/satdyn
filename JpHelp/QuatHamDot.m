function [qout_scal,qout_vec] = QuatHamDot(q1_Scal,q1_Vec,q2_Scal,q2_Vec)
    % implements Hamilton's product = quat dot product (NOT 4D vector inner
    % product): returns q1 (dot) q2  (Not commutative!!)
    qout_scal = q1_Scal*q2_Scal - dot(q1_Vec,q2_Vec);
    qout_vec = q1_Scal*q2_Vec + q2_Scal*q1_Vec + cross(q1_Vec,q2_Vec);

end

