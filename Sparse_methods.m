classdef Sparse_methods
    methods(Static)
        
        %%%%%%% TINNEY ZERO %%%%%%%%%%%%%
        function [ordering_idx,original_idx] = ...
                Tinney_zero_ordering(A_next, A_first)
            
            % calculate the number of rows/cols in the Sparse table
            N = max(A_next(:,3));
                        
            degree = zeros(1,N);
            for i=1:N
                for j=1:N
                    degree(i) = degree(i)+(sparse_table.retrieve...
                        (A_next,A_first,i,j)~=0);
                end
            end
            ordering_idx = (sortrows([(1:N)',degree'],2));
            % ordered one
            ordering_idx = ordering_idx(:,1);
            
            % original one
            original_idx = ((sortrows([(1:length(ordering_idx))',...
                ordering_idx],2))*[1;0])';
            
            % for test
            [A_ordered_next,A_ordered_first] = sparse_table.Blank_array();
            for i=1:max(A_next(:,3))
                for j=1:max(A_next(:,3))
                    [A_ordered_next,A_ordered_first] = ...
                        sparse_table.store(A_ordered_next,...
                        A_ordered_first,sparse_table.retrieve(A_next,...
                        A_first, ordering_idx(i),ordering_idx(j)),i,j);
                end
            end 
        end
    end
end