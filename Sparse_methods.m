classdef Sparse_methods
    methods(Static)
        
        %%%%%%%%%%%%%%%%%%%%%%% TINNEY ZERO %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [ordering_idx,original_idx,A_ordered_next] = ...
                Tinney_zero_ordering(A_next, A_first)
            
            % calculate the number of rows/cols in the Sparse table
            N = max(A_next(:,3));
                        
            % calculating the degree of each node
            degree = zeros(1,N);
            for i=1:N
                for j=1:N
                    degree(i) = degree(i)+(sparse_table.retrieve...
                        (A_next,A_first,i,j)~=0);
                end
            end
            
            ordering_idx = (sortrows([(1:N)',degree'],2));
            
            % ordered indices
            ordering_idx = ordering_idx(:,1);
            
            % original indices
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
        
        %%%%%%%%%%%%%%%%%%%%%%% TINNEY ONE %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ordering_idx,original_idx,A_ordered_next] = ...
                Tinney_one_ordering(A_next, A_first)
            % count the total number of rows/cols in the sparse table
            A_next_orig = A_next;
            A_first_orig = A_first;
            N = max(A_next(:,3));
            bus_number = (1:N);
            storage = zeros(1,N);
            for i = 1 : N
                
                % calculate 'N' as sparse table will be updated dynamically
                N1 = max(A_next(:,3));
                
                degree = zeros(1,N1);
                for k = 1 : N1
                    for j = 1 : N1
                        % step 1: calculate the degree of all vertices 
                        degree(k) = degree(k) + (sparse_table.retrieve...
                        (A_next,A_first,k,j)~=0);
                    end
                end
                pseudo_order = (sortrows([bus_number',...
                    (1:max(A_next(:,3)))',degree'],3));
                
                % step 2: choose the node with the lowest degree
                storage(i) = pseudo_order(1,1);
                
                B_next = A_next;
                B_first = A_first;
                
                for j = 1 : N1
                    for k = 1 : N1
                        % takes care of the fills that may occur after
                        % removing the vertex
                        [A_next,A_first] = sparse_table.store(A_next,...
                            A_first,sparse_table.retrieve(B_next,...
                            B_first,j,k) + sparse_table.retrieve(B_next,...
                            B_first,j, pseudo_order(1,2))*...
                            sparse_table.retrieve(B_next, B_first,...
                            pseudo_order(1,2),k),j,k);
                    end
                end
                
                % eliminate selected degree / choose natural order for tie
                pseudo_order_2 = find(sortrows(pseudo_order(:,2))...
                    ~= pseudo_order(1,2));
                
                % step 3: update the degrees
                [X_next,X_first]=sparse_table.Blank_array();
                for k=1:length(pseudo_order_2)
                    for j=1:length(pseudo_order_2)
                        [X_next,X_first] = sparse_table.store(X_next,...
                            X_first, (sparse_table.retrieve(A_next,...
                            A_first, pseudo_order_2(k),...
                            pseudo_order_2(j))),k,j);
                    end
                end
                
                % The original table will no longer have the removed node
                A_next = X_next;
                A_first = X_first;
                
                % remove the node from the bus
                bus_number(pseudo_order(1,2)) = [];
            end
            
            % ordered one
            ordering_idx = storage';
            
            % original one
            original_idx = ((sortrows([(1:length(ordering_idx))',...
                ordering_idx],2))*[1;0])'; 
            
            % for test
            [A_ordered_next,A_ordered_first] = sparse_table.Blank_array();
            for i=1:max(A_next_orig(:,3))
                for j=1:max(A_next_orig(:,3))
                    [A_ordered_next,A_ordered_first] = ...
                        sparse_table.store(A_ordered_next,...
                        A_ordered_first,sparse_table.retrieve(...
                        A_next_orig, A_first_orig, ordering_idx(i),...
                        ordering_idx(j)),i,j);
                end
            end             
        end        
    end
end