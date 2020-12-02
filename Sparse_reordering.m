function [ordering_idx, original_idx] = ...
    Sparse_reordering(J_next, J_first, Scheme)
    switch Scheme
        case 0
            % TINNEY ZERO
            [ordering_idx,original_idx] = ...
                Sparse_methods.Tinney_zero_ordering(J_next, J_first); 
%         case 1
%             % TINNEY ONE
%             [ordering_idx,original_idx] = ...
%                 Sparse_methods.Tinney_one_ordering(J_next, J_first)            
%         case 2
%             % TINNEY TWO
%             [ordering_idx,original_idx] = ...
%                 Sparse_methods.Tinney_two_ordering(J_next, J_first)
    end
end