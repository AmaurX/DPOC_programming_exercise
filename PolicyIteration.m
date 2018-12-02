function [ J_opt, u_opt_ind ] = PolicyIteration( P, G )
%POLICYITERATION Value iteration
%   Solve a stochastic shortest path problem by Policy Iteration.
%
%   [J_opt, u_opt_ind] = PolicyIteration(P, G) computes the optimal cost and
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

K = size(P,1);
err = 1e-5;

%% Policy Iteration

% Initialize with random proper policy
curr_policy = 5*ones(K,1); % Go up for all states
next_policy = 5*ones(K,1);

J_curr = zeros(K,1);
J_candidate = zeros(size(G));

P_mat = zeros(K,K);
G_vec = zeros(K,1);

while true
    % Policy Evaluation 
    for curr_state = 1:K % for every state i
        for next_state = 1:K % for every state j
            P_mat(curr_state,next_state) = P(curr_state,next_state,curr_policy(curr_state));
            G_vec(curr_state) = G(curr_state,curr_policy(curr_state));
        end
    end
    J_curr = (eye(K) - P_mat)\G_vec;

    % Policy Improvement
    for u = 1:5
        for curr_state = 1:K % for every state i
                for next_state = 1:K % for every state j                
                    P_mat(curr_state,next_state) = P(curr_state,next_state,u);
                    G_vec(curr_state) = G(curr_state,u);
                end
        end
            J_candidate(:,u) = G_vec + P_mat*J_curr;
    end
    
    for curr_state = 1:K
        [~, next_policy(curr_state)] = min(J_candidate(curr_state,:));
    end
    
    if all(abs(next_policy - curr_policy) <= err)        
        break;
    end
    
    curr_policy = next_policy;
end

J_opt = J_curr;
u_opt_ind = curr_policy;

end

