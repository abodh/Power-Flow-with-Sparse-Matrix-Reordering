classdef sparse_table
    methods(Static)
        
        % returns a blank array  
        function [A,B] = Blank_array()
            A = [];
            B = [];
        end
        
        % function to store the elements in the sparse tabular format
        function [A,B] = store(A,B,value,row,col)
            if value == 0
                % does not do anything when the value is 0 
                return

            elseif isempty(A) == 1
                % initialize the storage table for the first entry
                A=[1,value,row,col,0,0];
                B(row,:) = [row,1,0];
                B(col,:) = [col,0,1];                
                if row == col
                    % when diagonal
                    B(row,:) = [row,1,1];
                end 
                
                return
            
            elseif sparse_table.retrieve(A,B,row,col) ~= 0
                A((A(:,3) == row) & (A(:,4) == col),2) = value;
                
            else
                % identify the element just before the specified 
                % element in the row
                prev_in_row = max(A((A(:,3) == row)&(A(:,4) < col),4));
                if isempty(prev_in_row)
                    % if no element before the specified element 
                    % prev_in_row = 0
                    prev_in_row = 0;
                end

                % identify the element just after the specified 
                % element in the row
                next_in_row = min(A((A(:,3) == row)&(A(:,4) > col),4));
                if isempty(next_in_row)
                    % if no element before the specified element 
                    % next_in_row = 0
                    next_in_row = 0;
                end

                % identify the element just before the specified 
                % element in the col
                prev_in_col = max(A((A(:,3) < row) & (A(:,4) == col),3));
                if isempty(prev_in_col)
                    % if no element before the specified element 
                    % prev_in_col = 0
                    prev_in_col = 0;
                end

                % identify the element just after the specified 
                % element in the col
                next_in_col = min(A((A(:,3) > row) & (A(:,4) == col),3));
                if isempty(next_in_col)
                     % if no element before the specified element 
                     % next_in_col = 0
                    next_in_col=0;
                end

                % next in row in table A
                NIR = A((A(:,3) == row) & (A(:,4) == next_in_row),1);
                if isempty(NIR)
                    NIR=0;
                end

                % next in column in table A
                NIC = A((A(:,3) == next_in_col) & (A(:,4) == col),1);
                if isempty(NIC)
                    NIC=0;
                end

                % increase the counter for the index and insert the 
                % elements in A
                A=[A;size(A)*[1;0]+1,value,row,col,NIR,NIC];

                % if any prev in row or prev in col elements are 
                % identified then NIR and NIC for the previous element 
                % is updated accordingly
                A(A(:,3) == row & A(:,4) == prev_in_row,5) = size(A)*[1;0];
                A(A(:,3) == prev_in_col & A(:,4) == col,6) = size(A)*[1;0];

                % the snippet below fills the FIR and FIC
                if prev_in_row == 0
                    % if there is no prev in row then this is the 
                    % first element
                    B(row,1) = row;
                    B(row,2) = size(A)*[1;0];
                end

                if prev_in_col == 0
                    % if there is no prev in col then this is the 
                    % first element
                    B(col,1) = col;
                    B(col,3) = size(A)*[1;0];
                end
            end
        end

        % function to retrieve the elements from the sparse tabular format
        function value = retrieve(A,B,row,col)
            % check for empty or out of bounds
            if isempty(A) || max(row,col) > size(B)*[1;0] 
                % return 0 if the storage table is empty or out of bound
                value = 0;
                return
            end
            FIR = B(row,2);
            FIC = B(col,3);
            if FIR == 0 || FIC == 0
                % return 0 if either of FIR or FIC is 0
                value = 0;
                return
            end
            
            % infinite loop until one of the conditions are met
            while 1                
                if A(FIR,4)== col
                    % ensures that the FIR element is extracted
                    value = A(FIR,2);
                    return
                    
                elseif A(FIC,3) == row
                    % ensures that the FIC element is extracted
                    value = A(FIC,2);
                    return
                    
                elseif A(FIR,5) == 0 || A(FIC,6) == 0
                    % ensures that we stop if we went ahead from last one
                    value = 0;
                    return
                    
                else
                    % if nothing works move ahead by keeping NIx to FIx
                    FIR = A(FIR,5); % A(,5) is NIR
                    FIC = A(FIC,6); % A(,6) is NIC
                end
            end
        end
    end
end