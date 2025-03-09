function [scal_BA,vec_BA] = dcm2quat_jp(O_BA)
    % converts O_BA to quaternion's scaler and vector parts (valid for
    % scaler first/last)
    % quaternion = (scal_BA, vect_BA) which rotates frame A to frame B
    % quaternion can be scaler first/last, does not matter
    
    tr_O_BA = trace(O_BA);
    vec_BA1 = [1 + 2*O_BA(1,1) - tr_O_BA;
                O_BA(1,2) + O_BA(2,1);
                O_BA(1,3) + O_BA(3,1);];
    scal_BA1 = O_BA(2,3) - O_BA(3,2);
    norm1 = norm([vec_BA1;scal_BA1]);

    vec_BA2 = [O_BA(1,2) + O_BA(2,1);
                1 + 2*O_BA(2,2) - tr_O_BA;                
                O_BA(2,3) + O_BA(3,2);];
    scal_BA2 = O_BA(3,1) - O_BA(1,3);
    norm2 = norm([vec_BA2;scal_BA2]);

    vec_BA3 = [O_BA(3,1) + O_BA(1,3);
               O_BA(2,3) + O_BA(3,2);
                1 + 2*O_BA(3,3) - tr_O_BA];
    scal_BA3 = O_BA(1,2) - O_BA(2,1);
    norm3 = norm([vec_BA3;scal_BA3]);

    vec_BA4 = [ O_BA(2,3) - O_BA(3,2);
                O_BA(3,1) - O_BA(1,3);
                O_BA(1,2) - O_BA(2,1)];
    scal_BA4 = 1 + tr_O_BA; 
    norm4 = norm([vec_BA4;scal_BA4]);

    [~,Idx] = max([norm1,norm2,norm3,norm4]);

    switch Idx
        case 1
            scal_BA = scal_BA1/norm1;
            vec_BA = vec_BA1/norm1;
        case 2
            scal_BA = scal_BA2/norm2;
            vec_BA = vec_BA2/norm2;
        case 3
            scal_BA = scal_BA3/norm3;
            vec_BA = vec_BA3/norm3;
        case 4
            scal_BA = scal_BA4/norm4;
            vec_BA = vec_BA4/norm4;
    end
end

