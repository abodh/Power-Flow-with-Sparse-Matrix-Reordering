function [Jacob_sparse_ordered, ordering_idx, original_idx] = ...
    Sparse_reordering(Jacob_sparse, Scheme)
    
    % calculate the number of rows/cols in the Jacobian
    N_Jacobian = max(Jacob_sparse(:,3));
    
    switch Scheme
        case 0
            %%%%%%% TINNEY ZERO %%%%%%%%%%%%%
            degree=zeros(1,N_Jacobian);
            for i=1:N_Jacobian
                for j=1:N_Jacobian
                    degree(i)=degree(i)+(sparse_table.retrieve...
                        (Jacob_sparse,i,j)~=0);
                end
            end
            ordering_idx = (sortrows([(1:N_Jacobian)',degree'],2));
            ordering_idx = ordering_idx(:,1);

            [Jacob_sparse_ordered,JB]=sparse_table.Blank_array();
            for i=1:max(Jacob_sparse(:,3))
                for j=1:max(Jacob_sparse(:,3))
                    [Jacob_sparse_ordered,JB]=...
                        sparse_table.store(Jacob_sparse_ordered,...
                        JB,sparse_table.retrieve(Jacob_sparse,...
                        ordering_idx(i),ordering_idx(j)),i,j);
                end
            end
            original_idx = ((sortrows([(1:length(ordering_idx))',...
                ordering_idx],2))*[1;0])';
    end
end