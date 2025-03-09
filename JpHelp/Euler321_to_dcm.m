function DCM = Euler321_to_dcm(yawR,pitchR,rollR)
    % convert 321 yaw pitch roll euler sequence to DCM
    % DCM = O_BA, the orientation matrix to change perspective X_B =
    % O_BA*X_A, where X_B is vector X resolved in frame B
    cy = cos(yawR); sy = sin(yawR);
    cp = cos(pitchR); sp = sin(pitchR);
    cr = cos(rollR); sr = sin(rollR);

    DCM = [cy*cp sy*cp -sp;
           cy*sp*sr-sy*cr sy*sp*sr+cy*cr cp*sr;
           cy*sp*cr+sy*sr sy*sp*cr-cy*sr cp*cr];
end

