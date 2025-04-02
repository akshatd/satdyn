function qi = SlerpJP(q0,q1,t)
    % implements the SLERP method.
    % t within [0,1], interpolate from q0 to q1 along 4D hypersphere
    % because this operates by treating quat as just 4D normalized vectors,
    % it does NOT matter if is scaler first/last, so long both q0 and q1
    % are consistent
    cOmega = sum(q0.*q1);
    Omega = acos(cOmega);
    qi = 1/sin(Omega)*(sin(Omega - t*Omega)*q0 + sin(t*Omega)*q1);
end

