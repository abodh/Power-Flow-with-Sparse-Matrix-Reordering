%% check
A = [1 0 -2 0 0; 2 8 0 1 0; 0 0 3 0 -2; 0 -3 2 0 0; 1 2 0 0 -4];
B = A'
[col,row,V] = find(B)
B(row(1):end,col(1))



%%
clear;
clc;
B = [1 0 -2 0 0; 2 8 0 1 0; 0 0 3 0 -2; 0 -3 2 0 0; 1 2 0 0 -4];
A = B';
storage_table = zeros(nnz(A),6);

%{
* B = real sparse matrix  
A = transpose of sparse matrix B
Transpose is done here as some functions in matlab access the elements in
column order

* storage table is the table of the form:
Index | Value | NROW | NCOL | NIR | NIC

* Index = index of the stored non-zero element
    value :  value of the non-zero element in A
     NROW :  row index of non-zero element in A
     NCOL :  column index of non-zero element in A
      NIR :  the index of the next non-zero element in the row
      NIC : the index of the next non-zero element in the column 
%}

storage_table(:,1) = 1:nnz(A);
storage_table(:,2) = nonzeros(A);
[storage_table(:,4), storage_table(:,3)] = find(A);

storage_table(:,5) = next_in_row(B, storage_table);
% storage_table(:,6) = next_in_column(B, storage_table)

function NIR = next_in_row(B, storage_table)
    idx=1;
    visited = [];
    for i = 1:length(B)
        [row_R , col_R] = find(B(i,:)); 
        row_R = row_R + i - 1;
            for k = 1 : (length(col_R)+1)
                if k > length(col_R)
                    pseudo_table_NIR(idx) = 0;
                    idx = idx + 1;
                    break
                end                
                idx_search = find(storage_table(:,3)==row_R(k) & storage_table(:,4)==col_R(k));
                if k > 1 && sum(visited == idx_search)== 0 
                    pseudo_table_NIR(idx) = idx_search;
                    idx = idx + 1;
                end      
                visited = [visited idx_search];
            end
    end
    NIR = pseudo_table_NIR;
end

% function NIC = next_in_column(B, storage_table)
%     pseudo_table_NIC = storage_table(:,6);
%     idx=1;
%     visited = [];
%     for i = 1:length(B)
%         [row_R , col_R] = find(B(i,:));
%         row_R = row_R + i - 1;
%         [row_C, col_C] = find(B(:,i));
%         row_C = row_C'; 
%         col_C = col_C';
%         col_C = col_C + i - 1;
%             for k = 1 : (length(col_R)+1)
%                 if k > length(col_R)
%                     pseudo_table_NIC(idx) = 0;
%                     idx = idx + 1;
%                     break
%                 end                
%                 idx_search = find(storage_table(:,3)==row_C(k) & storage_table(:,4)==col_C(k));
%                 if k > 1 & sum(visited == idx_search)== 0 
%                     pseudo_table_NIC(idx) = idx_search;
%                     idx = idx + 1;       
%                 end      
%                 visited = [visited idx_search];
%             end
%     end
%     NIC = pseudo_table_NIC;
% end