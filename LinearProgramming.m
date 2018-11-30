function [ J_opt, u_opt_ind ] = LinearProgramming( P, G )
%LINEARPROGRAMMING Value iteration
%   Solve a stochastic shortest path problem by Linear Programming.
%
%   [J_opt, u_opt_ind] = LinearProgramming(P, G) computes the optimal cost
%   and the optimal control input for each state of the state space.
%
%   Input arguments:
%
%       P:
%           A (K x K x L)-matrix containing the transition probabilities
%           between all states in the state space for all control inputs.
%           The entry P(i, j, l) represents the transition probability
%           from state i to state j if control input l is applied.
%
%       G:
%           A (K x L)-matrix containing the stage costs of all states in
%           the state space for all control inputs. The entry G(i, l)
%           represents the cost if we are in state i and apply control
%           input l.
%
%   Output arguments:
%
%       J_opt:
%       	A (K x 1)-matrix containing the optimal cost-to-go for each
%       	element of the state space.
%
%       u_opt_ind:
%       	A (K x 1)-matrix containing the index of the optimal control
%       	input for each element of the state space.

% put your code here

% Note: we can use linprog to solve the linear programming problem

K = length(G(:,1));
U_optimal = ones(K,1);
% cost_to_go_candidates = zeros(length(P(1,1,:)),1);

% x = V
% f^T = -ones(1,K)
% A
% b = q(i,u) + cost_sum <- Note: This will be (5*K) x 1 because there are
% five inputs for each possible state and this five inequality constraaints
% for each state.
% x = linprog(f,A,b)

U = length(P(1,1,:)); % Number of inputs
b = zeros(K*U,1);
A = zeros(K*U,K);
A_ones = zeros(K*U,K);

idx = 0;
for current_state_row = 1:K
    for input_num = 1:U
        idx = idx+1;
        A_ones(idx,current_state_row) = 1;
    end
end

idx = 0;
for current_state_row = 1:K
    i = current_state_row;
        for u = 1:U
            idx = idx + 1;
            g_iu = G(i, u);
                for state_prob_col = 1:K
                    j = state_prob_col;
                    p_ij_u = P(i, j, u);
                    
                    if g_iu == Inf % Because the linprog solver throws an
                                   % error if the stage cost is = Inf
                        g_iu = 10^20;
                    end
                    b(idx) = g_iu;
                    A(idx,j) = -p_ij_u;
                end
        end
end
A = A_ones + A;
f = -ones(1,K);

x = linprog(f,A,b);

idx = 0;
for current_state_row = 1:K
    i = current_state_row;
        for u = 1:U
            idx = idx+1;
            g_iu = G(i,u);
            if g_iu == Inf % Because the linprog solver throws an
                           % error if the stage cost is = Inf
                g_iu = 10^20;
            end
            if abs(x(i) - (P(i, 1:end, u)*x + g_iu)) < 10^-5
                U_optimal(i) = u;
            end
%             disp([abs(x(i) - (P(i, 1:end, u)*x + g_iu)) < 10^-5, x(i), P(i, 1:end, u)*x + g_iu, u, U_optimal(i)])
        end
end

J_opt = x;
u_opt_ind = U_optimal;

end

