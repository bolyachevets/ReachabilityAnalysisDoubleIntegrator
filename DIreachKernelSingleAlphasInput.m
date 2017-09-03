% Double Integrator Model
% approximation of finite horizon viability kernel
% -------------------------------------------------------------------------
% scaling is performed only once at the beginning of set evolution:
% within a single optimization problem both prescaling of the initial set 
% that guarantees constraint satisfaction and consequent controls are 
% computed simultaneously.
% mu variables represent hypothetical input dimension of the augmented state space. 
% Each mu_c (center) and mu_x (radius) belong to the corresponding
% generator. Similarly, alpha_c and alpha_x represent the center and radius
% of the rescaling of initial set.

% parameters of interest:
% - timeStep and tFinal
% - IC and CS (centered or skewed)
% n: state space dimension
% p: number of generators
% m: control space dimension

% TIME
%----------------------------------------------------------------------
tStart=0; %start time
tFinal=3; %final time
timeStep=0.1; %time step size for reachable set computation
% number of linearization points
number_steps = (tFinal-tStart)/timeStep;
%----------------------------------------------------------------------

% INITIAL CONDITIONS
%----------------------------------------------------------------------
% dimension of state space
n = 2;
IC = interval([-1;-1], [0;0]);
% dimension of control space
m = 1;

IC_z = zonotope(IC);

IC_z_generators = get(IC_z, 'Z');
IC_z_generators = IC_z_generators(:, 2:length(IC_z_generators));
IC_c = center(IC_z);

% add extra generators
%  p = n;
d_extra = 10;
p = d_extra + n;
extraScale = 1/d_extra;
extraG = extraScale*generatePts(d_extra, n);
IC_z_generators = horzcat(IC_z_generators, extraG);
IC_mat = horzcat(IC_c, IC_z_generators);
IC_z = zonotope(IC_mat);
%----------------------------------------------------------------------

% CONSTRAINTS
%----------------------------------------------------------------------
%CS = interval([-2;-2], [2;2]);
CS = IC;
%CS_z = zonotope(CS);
%----------------------------------------------------------------------

% MATRICES
%----------------------------------------------------------------------
A = [0   1;
     0   0];

B = [0; 1];

A_d = expm(A*timeStep);
integrandB = @(tau) expm(A*tau);
B_d = integral(integrandB, 0, timeStep, 'ArrayValued', true)*B;
%----------------------------------------------------------------------

% OPTIMIZATION
%----------------------------------------------------------------------
% makes sure that interval are not degenerate
epsilon = 0.01;

% Weight the inputs by the number of input dimensions to account for the
% extra degrees of freedom compared to the state scaling.  Could also
% downweight the inputs if you care more about the state space and/or
% separately weight the input dimensions if you care more about some
% inputs.
input_weights = ones(m, 1) / m;
     
        

cvx_begin quiet
            variable alpha_c(p)
            variable alpha_x(p)
            variable mu_c(m, p*number_steps )
            variable mu_x(m, p*number_steps )
            expression G_ij(n,p)
            
            maximize sum(alpha_x) + sum(input_weights .* sum(mu_x, 2));

            subject to
            % constraints on the magnitude of alpha_c, alpha_x can be
            % easily derived from the property that scaling factors must
            % exceed epsilon (if we want scaled generator intervals to be
            % non-empty)
            for j=1:p
                % assuming gamma_upper_bar >= gamma_lower_bar + epsilon and that gamma_lowe_bar = alpha_c - alpha_x
                % we can derive the following constraints on alphas:
                alpha_x(j) >= epsilon;
                % make sure that scaling factors are positive
                alpha_c(j) - alpha_x(j) >= 0;  
                
            %  (can only scale down) seems not necessary with a single optimization problem                
%                 alpha_x(j) <= 1;
%                 alpha_c(j) >=  alpha_x(j) - 1;
%                 alpha_c(j) <= 1 - alpha_x(j);
            end
           
           % individual control effects on generators are bounded by 1, 
           % where -1/+1 do not have the interpretation of a scaling factor, as
           % in alpha, but rather are min/max values for input
           for j=1:m
             for i=1:p*number_steps
%                 mu_c(j,i) - mu_x(j,i) >= - 1;
%                 mu_c(j,i) + mu_x(j,i) <= 1;
                % control is allowed to be 0
                mu_x(j,i) >= 0;
                % maximum radius is half the length of input interval
                % constraint: (1-(-1))/2
%                 mu_x(j,i) <= 1;
             end
           end
           % cumulative effect of control cannot exceed 1 in absolute value
           % (where 1 stands for the most extreme input value allowed, 
           % which in this case is symmetric around 0)
           for i=0:(number_steps-1)
               sum(mu_c(:, (p*i+1):p*(i+1)),2) - sum(mu_x(:, (p*i+1):p*(i+1)),2) >=  -1;
               sum(mu_c(:, (p*i+1):p*(i+1)),2) + sum(mu_x(:, (p*i+1):p*(i+1)),2) <= 1;
               % cumulative radius for a given control dimension is between
               % zero and half the distance between min and max input
               % values
%                sum(mu_x(:, (p*i+1):p*(i+1)),2) <= 1;
%                sum(mu_x(:, (p*i+1):p*(i+1)),2) >= 0;             
           end
          % scaled initial set is inside unscaled one, need this because
          % of extra generators added 
         IC_c + IC_z_generators*alpha_c ...
         + abs(IC_z_generators)*alpha_x  <= supremum(IC);
         IC_c + IC_z_generators*alpha_c ...
         - abs(IC_z_generators)*alpha_x >= infimum(IC);
         
        % prescaling generators
        for j=1:p
                G_ij(:,j) = alpha_x(j)*IC_z_generators(:,j);
        end
                   
         for t=1:number_steps
             (A_d^t)*(IC_c + IC_z_generators*alpha_c) ... 
                 + sum(computeInputEffect(A_d, B_d, mu_c, p, t), 2) ...
                 + sum(abs(A_d^t* G_ij + computeInputEffect(A_d, B_d, mu_x, p, t)),2) <= supremum(CS);
             (A_d^t)*(IC_c + IC_z_generators*alpha_c) ... 
                 + sum(computeInputEffect(A_d, B_d, mu_c, p, t), 2) ...
                 - sum(abs(A_d^t* G_ij + computeInputEffect(A_d, B_d, mu_x, p, t)),2) >= infimum(CS);
         end

 cvx_end
     
 alpha_c
 alpha_x

% ACCUMULATE SCALARS - APPLY TO INITIAL SET
%----------------------------------------------------------------------
IC_g_mat =[];
IC_z_mat =[];
for j=1:p
    tempGenerator = IC_z_generators(:,j);
    tempGenerator = alpha_x(j)*tempGenerator;
    IC_g_mat = horzcat(IC_g_mat, tempGenerator);
end
IC_z_mat = horzcat(IC_c + IC_z_generators*alpha_c, IC_g_mat);
IC_z_g = zonotope(IC_z_mat);
%----------------------------------------------------------------------

% PLOTS
%----------------------------------------------------------------------
figure;
hold on;

plot(IC, [1,2], 'y','lineWidth',2);
plot(CS, [1,2], 'g','lineWidth',2);
plot(IC_z_g, [1,2],'y','lineWidth',2);
%pause(1);
 

% evolve the safe set forward to illustrate the outcome (proof of concept)
Reach_generators = IC_g_mat;
for i=1:number_steps
    Reach_set_gen = [];
    controlEffect = computeInputEffect(A_d, B_d, mu_x, p, i);
    for j=1:p 
        tempGenerator = Reach_generators(:,j);
        tempGenerator = (A_d^i)*tempGenerator ...
            + controlEffect(:,j);
        Reach_set_gen = horzcat(Reach_set_gen, tempGenerator);
    end
    Reach_center = A_d^i*center(IC_z_g) + sum(computeInputEffect(A_d, B_d, mu_c, p, i),2);

    plot(zonotope(horzcat(Reach_center, Reach_set_gen)), [1,2], 'g', 'lineWidth', 2);
    %pause(1);
end