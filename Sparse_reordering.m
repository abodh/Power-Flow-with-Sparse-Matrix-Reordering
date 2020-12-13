function [ordering_idx, original_idx, A_ordered_next] = ...
    Sparse_reordering(J_next, J_first, Scheme)
    switch Scheme
        case 0
            % TINNEY ZERO
            [ordering_idx,original_idx, A_ordered_next] = ...
                Sparse_methods.Tinney_zero_ordering(J_next, J_first); 
        case 1
            % TINNEY ONE
            [ordering_idx,original_idx, A_ordered_next] = ...
                Sparse_methods.Tinney_one_ordering(J_next, J_first);  
            
%         case 2
%             % TINNEY TWO
%             [ordering_idx,original_idx] = ...
%                 Sparse_methods.Tinney_two_ordering(J_next, J_first);
        otherwise
            ordering_idx = 1:max(J_next(:,3));
            original_idx = ((sortrows([(1:length(ordering_idx))',...
                ordering_idx'],2))*[1;0])';
            A_ordered_next = J_next;
    end
end