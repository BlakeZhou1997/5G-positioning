function parms_out = GMatrix2_Loosely_ECEF(R_b_e)

%% INPUTS
% R_b_e    ---> attitude DCM matrix(3x3) from body to ECEF frame

% returns the G matrix (5x4)


parms_out = [zeros(3),        zeros(3),      zeros(3),       zeros(3)      
              R_b_e,          zeros(3),      zeros(3),       zeros(3)       
              zeros(3)        R_b_e,          zeros(3),      zeros(3)      
              zeros(3),       zeros(3),        eye(3),       zeros(3)      
              zeros(3),       zeros(3),      zeros(3),        eye(3)   ];
%           
end