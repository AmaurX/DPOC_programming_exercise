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

epsilon = 1e-100;
K = length(G(:,1));
J_kplus1 = ones(K, 1);
J_k = ones(K, 1);
U_optimal = ones(K,1);

%%
% state_label_mat = zeros(1,K^2);
%     for idx = 1:K^2
%         state_label_mat(idx) = mod(idx,K+1);
%     end
% state_label_mat = reshape(state_label_mat,[K,K])';
% state_label_mat = triu(state_label_mat) + triu(state_label_mat)' - eye(K,K);
%%

cost_to_go_candidates = zeros(1,length(P(1,1,:)));
input_candidates = zeros(1,length(P(1,1,:)));
iter_vec = zeros(K,1);

for current_state_row = 1:K
    test = 0;
        while test == 0
            iter_vec(current_state_row) = iter_vec(current_state_row) + 1;
            i = current_state_row;
                for u = 1:length(P(1,1,:))
                    g_iu = G(i, u);
                        for state_prob_col = 1:K
                            j = state_prob_col;
                            p_ij_u = P(i, j, u);
                            cost_sum = cost_sum + p_ij_u*J_k(j);
                        end
                    cost_to_go_candidates(u) = g_iu + cost_sum;
                    cost_sum = 0;
                    input_candidates(u) = u;
                end
            J_kplus1(i) = min(cost_to_go_candidates);
            U_optimal(i) = find(cost_to_go_candidates == J_kplus1(i),1,'last');
            test = (abs(J_k(i) - J_kplus1(i)) <= epsilon);
            J_k = J_kplus1;
        end                
end
    
    J_opt = J_k;
    u_opt_ind = U_optimal;
%     disp(iter_vec)
    
end