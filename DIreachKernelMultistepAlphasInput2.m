% scaling is performed only once at the beginning of set evolution.
% there is also one optimization problem that computes prescaling and
% consequent controls simultaneously
% mu tilda here is  the cumulative effect of rescaled control in each
% generator

% parameters of interest:
% - timeStep and tFinal
% - IC and CS (centered or skewed)
% n: state space dimension
% p: number of generators
% m: control space dimension

% TIME
%----------------------------------------------------------------------
% potentially will need to get rid of the similar fields in the options
tStart=0; %start time
tFinal=6; %final time
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
CS = interval([-2;-2], [2;2]);
%CS = IC;
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
            for j=1:p
                alpha_c(j) >=  alpha_x(j) - 1;
                alpha_c(j) <= 1 - alpha_x(j);
                % make sure that additional generators don't degenerate
                alpha_x(j) >= epsilon;
                alpha_x(j) <= 1;
            end
           
           % individual control effects on generators are bounded by 1
           for j=1:m
             for i=1:p*number_steps
                mu_c(j,i) >= mu_x(j,i) - 1;
                mu_c(j,i) <= 1 - mu_x(j,i);
                % make sure that additional generators don't degenerate
                mu_x(j,i) >= epsilon;
                mu_x(j,i) <= 1;
             end
           end
           % cumulative effect of control cannot exceed 1 in absolute value
           for i=0:(number_steps-1)
               sum(mu_c(:, (p*i+1):p*(i+1)),2) >= sum(mu_x(:, (p*i+1):p*(i+1))) - ones(n,1);
               sum(mu_c(:, (p*i+1):p*(i+1)),2) <= ones(n,1) - sum(mu_x(:, (p*i+1):p*(i+1)));
               sum(mu_x(:, (p*i+1):p*(i+1)),2) <= ones(n,1);
               sum(mu_x(:, (p*i+1):p*(i+1)),2) >= -ones(n,1);
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
        
        alpha_x
        alpha_c
        mu_x
        mu_c

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
pause(1);
 

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



% %reachSetOptimization(center(IC_z_g), IC_z_mat(:,2:end), tStart, tFinal, timeStep)
% % plot the center of the scaled initial set
% c = center(IC_z_g);
% x = c(1);
% y = c(2);
% plot(x, y, 'x');