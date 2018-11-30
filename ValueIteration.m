function [ J_opt, u_opt_ind ] = ValueIteration( P, G )
%VALUEITERATION Value iteration
%   Solve a stochastic shortest path problem by Value Iteration.
%
%   [J_opt, u_opt_ind] = ValueIteration(P, G) computes the optimal cost and
%   the optimal control input for each state of the state space.
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

% Notes:
% J_k+1 (i) min_(u in U(i)) {g(i, u) + sum_(j=1)^n p_ij (u) J_k (i)}
tic
epsilon = 1e-5;
K = length(G(:,1));
J_kplus1 = ones(K, 1);
J_k = ones(K, 1);
U_optimal = ones(K,1);

cost_to_go_candidates = zeros(1,length(P(1,1,:)));
test = 0;

while test == 0
    for current_state_row = 1:K
            i = current_state_row;
                for u = 1:length(P(1,1,:))
                    cost_sum = 0;
                    g_iu = G(i, u);
                        for state_prob_col = 1:K
                            j = state_prob_col;
                            p_ij_u = P(i, j, u);
                            cost_sum = cost_sum + p_ij_u*J_k(j);
                        end
                    cost_to_go_candidates(u) = g_iu + cost_sum;
                end
            [J_kplus1(i), u_opt_idx] = min(cost_to_go_candidates);
            U_optimal(i) = u_opt_idx;
    end
    test = all(abs(J_k - J_kplus1) <= epsilon);
    J_k = J_kplus1;
end
    
    J_opt = J_k;
    u_opt_ind = U_optimal;
toc
end