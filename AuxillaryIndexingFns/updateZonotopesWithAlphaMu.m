% scaling and reach set updating for implementation with lower and upper
% bound scalars on generators

function  [IC_z_generators_copy, alphas, mus, mu_cs, center_shift_x, scaledICs, R_zs, IC_z_copy, IC_c_copy] = ...
           updateZonotopesWithAlphaMu(alpha_c, alpha_x, mu_c, mu_x, IC_z_generators_copy, IC_c_copy, ...
           alphas, mus, mu_cs, center_shift_x, p, s, scaledICs, A_d, B_d, R_zs)
        % store scalars from the current iteration
        for j=1:p
            alphas(generatorIndex(p, s)+j) = alpha_x(j);
        end

        % store cumulative control effect on generators
        mus{s} = mu_x;
        mu_cs{s} = mu_c;
           
        % update cumulative center shift
        center_shift_x = center_shift_x + (A_d^(-1))^(s-1)*(IC_z_generators_copy*alpha_c);
%         center_shift_x = center_shift_x + (A_d^(-1))^(s-1)*(IC_z_generators_copy*alpha_c + sum(A_d^(-1)*B_d*mu_c,2));
        
        % update the center of the scaled initial set
        IC_c_copy = IC_c_copy + IC_z_generators_copy*alpha_c;
        
        % update the generator matrix to construct scaling of the initial set      
        IC_z_g =[];
        for i=1:p
            IC_z_g = horzcat(IC_z_g, alpha_x(i)*IC_z_generators_copy(:,i));
        end
        IC_z_copy_mat = horzcat(IC_c_copy, IC_z_g);
        IC_z_copy = zonotope(IC_z_copy_mat);
        scaledICs = horzcat(scaledICs, IC_z_copy);

        % construct a reach set
        reachSetGenerators = generatorsWithControls(A_d, B_d, IC_z_generators_copy, alpha_x, mu_x);
        reachSet = zonotope(horzcat(A_d*IC_c_copy + sum(B_d*mu_c,2), reachSetGenerators));
        R_zs = horzcat(R_zs, reachSet);
        
        % the next initial condition is the reach set from this
        % iteration 
        IC_c_copy = center(reachSet);
        IC_z_generators_copy = reachSetGenerators;
        IC_z_copy = reachSet;
           
end