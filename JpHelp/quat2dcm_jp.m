function DCM = quat2dcm_jp(scal_BA,vect_BA)
    % returns the DCM based on quad's scaler and vector parts
    % quaternion = (scal_BA, vect_BA) which rotates frame A to frame B
    % quaternion can be scaler first/last, does not matter
    % DCM = O_BA, the orientation matrix to change perspective X_B =
    % O_BA*X_A, where X_B is vector X resolved in frame B
    v1 = vect_BA(1); v2 = vect_BA(2); v3 = vect_BA(3);
    s = scal_BA;
    
    DCM = [(s^2 + v1^2 - v2^2 - v3^2) 2*(v1*v2 + v3*s) 2*(v1*v3 - v2*s);
           2*(v1*v2 - v3*s) (s^2 - v1^2 + v2^2 - v3^2) 2*(v1*s + v2*v3);
           2*(v1*v3 + v2*s) 2*(v2*v3 - v1*s) (s^2 - v1^2 - v2^2 + v3^2)];
end