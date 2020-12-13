%{
Newton Rhapson Power Flow with Sparse Re-ordering
Author: Abodh Poudyal
Last updated: Dec 1, 2020
%}
clear;
clc;
format short % to display less significant digits in the result 

%% Initializations

% select the system on which the power flow should be performed
% currently works for 'IEEE14' and 'IEEE30'
System = "IEEE-30";

% iterate unless power mismatch < 0.01 (tolerance)
tolerance = 0.01;

% base power
base_MW = 100;

%{
    Sparse Ordering Scheme
    0 : Tinney zero
    1 : Tinney one
    2 : Tinney two (not completed yet)
%}

Scheme = 1;

if Scheme < 2
    fprintf("Tinney %d ordering scheme is being used ... \n", Scheme)
else 
    fprintf("No ordering is being done \n")
end

%% 1. Reading bus and branch data in common data format

% external function to extract the data from IEEE common data format
bus_path = strcat(System,'bus_data/bus_data.txt');
branch_path = strcat(System,'bus_data/branch_data.txt');
[bus_imp, branch_imp, bus_data, branch_data] = ...
    data_extraction(bus_path, branch_path,System);

%{ 
to reduce the computational complexity, we will only compute 
for existing branches
%}

% define some important variables

% from which bus
from = branch_data(:,1);
% to which bus
to = branch_data(:,2);

% extract voltage data
V_flat = bus_data.data(:,11);  

% flat start means |V| = 1.0 pu and delta = 0
% exact values for PV and slack bus whereas flat start for the rest
V_flat(find(V_flat == 0)) = 1;
delta_flat = zeros(length(V_flat),1);

% number of buses in the entire system
n_bus = length(bus_data.data(:,3));

% number of branches
n_branch = length(branch_imp);

%% 2. Calculating the Y-bus matrix

%{
* storage table is the table of the form:
Index | Value | NROW | NCOL | NIR | NIC

* For A_next: 
    Index : index of the stored non-zero element
    value : value of the non-zero element in A
     NROW : row index of non-zero element in A
     NCOL : column index of non-zero element in A
      NIR : the index of the next non-zero element in the row
      NIC : the index of the next non-zero element in the column 

* For A_first: 
    Index : index of the stored non-zero element
      FIR : the index of the first non-zero element in the row
      FIC : the index of the first non-zero element in the column 

Reference:
"Computational Methods for Electric Power Systems" by Mariesa L. Crow for
in-depth explanation"
 
%}

% Y-bus sparse storage
[Y_next, Y_first] = Ybus_sparse(n_bus, n_branch, branch_imp, bus_imp,...
                                from, to);
                            
%% 3,4,5. Calculating Jacobian Matrix, LU factorization, and NR power flow

Q_lim_status = 1;
while (Q_lim_status)
    
    % scheduled power
    Ps = (bus_data.data(:,8) - bus_data.data(:,6))/base_MW;
    Qs = (bus_data.data(:,9) - bus_data.data(:,7))/base_MW;

    % number of pq buses
    n_pq = length(find(bus_data.data(:,3) == 0));

    % number of PV buses 
    n_pv = length(find(bus_data.data(:,3) == 2));

    % stores an array of PQ bus IDs
    pq_bus_id = find(bus_data.data(:,3) == 0);

    % stores an array of PV bus IDs 
    pv_bus_id = find(bus_data.data(:,3) == 2);

    % Newton Rhapson Power Flow 
    [Volt, Angle, error_avg, L_next, U_next, J_next, J_ordered_next] = ...
        NewtonRhapson(tolerance, n_bus, n_pv, n_pq, pq_bus_id, V_flat, ...
        delta_flat, Y_next, Y_first, Ps, Qs, Scheme);
    
    V_final = Volt(:,end);
    Angle_final = Angle(:,end);
    
    plot_states(Volt, Angle)
    
    % Q-limit check
    [Q_lim_status, bus_data] = Qlim(V_final, Angle_final, bus_data, ...
        Y_next,Y_first, base_MW, pv_bus_id, n_bus);
    if (Q_lim_status)
        V_flat = V_final;
        delta_flat = Angle_final;
    end
end

%{
 calculated as -> nnz in Q(L+U) - nnz in J - J (subtracting common non 
 zeros that is already in J) 
%}
Fills = (size(L_next)+size(U_next)-size(J_next))*[1;0] - max(J_next(:,3));

%% Results
fprintf("\n")
fprintf("---------------------------------------------\n")
fprintf("--------------- SUMMARY ---------------------\n")
fprintf("---------------------------------------------\n")
fprintf("System = %s bus system \n",System)
fprintf("Total non-zeroes before fills = %d \n",size(J_next)*[1;0]);
fprintf("Total fills = %d \n",Fills);
fprintf("Total non-zeroes after fills = %d \n",Fills + size(J_next)*[1;0]);
fprintf("%% of non-zeroes = %.2f%%",((Fills + size(J_next)*[1;0])/...
    (max(L_next(:,3))^2)*100));
fprintf("\n")

% fprintf("---------------------------------------------\n")
% fprintf("--------------- POWER FLOW ------------------\n")
% fprintf("---------------------------------------------\n")
% 
% fprintf("Voltage  Angle Voltage Angle \n")
% fprintf("---------------------------------------------\n")
% collect = [bus_data.data(:,4),bus_data.data(:,5),V_final,Angle_final*180/pi];
% collect

figure('color', [1,1,1])
aH=axes;
scatter(aH, J_ordered_next(:,3),J_ordered_next(:,4),40,'filled');
aH.YDir = 'reverse';
aH.XLim = [0 max(J_ordered_next(:,3))+1];
aH.YLim = [0 max(J_ordered_next(:,3))+1];
% aH.PlotBoxAspectRatio = [1 1 1];
title ('Jacobian')
grid minor

figure('color', [1,1,1])
aH=axes;
scatter(aH,[L_next(:,3);U_next(:,3)],[L_next(:,4);U_next(:,4)],40,'filled');
aH.YDir = 'reverse';
aH.XLim = [0 max(L_next(:,3))+1];
aH.YLim = [0 max(L_next(:,3))+1];
% aH.PlotBoxAspectRatio = [1 1 1];
title ('Q')
grid minor

% %% 6. Fast decoupled power flow
% [Volt_FD, Angle_FD, error_avg_FD] = ...
%     FastDecoupledPF(tolerance, from, to, n_branch, n_bus, n_pv, n_pq, ...
%     pq_bus_id, V_flat, delta_flat, Y_first,Y_next, Ps, Qs, branch_imp,...
%     bus_imp); 

%% plots of the result

function plot_states(Volt, Angle)
    % NRLF Voltage
    figure('color', [1,1,1])
    str = "bus 1";
    for i = 1: length(Volt)
        plot(Volt(i,:), 'Linewidth', 1.5)
        hold on
        if i > 1
            str = [str , strcat('bus',' ', num2str(i))];
        end
    end

    ylabel('Voltage (pu)')
    xlabel('Number of iteration')
    grid on
    set(gca,'XTick',(1:1:10))
    set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
    lgd = legend (str, 'NumColumns', 4);
    lgd.FontSize = 9;
    hold off

    % NRLF Angle
    figure('color', [1,1,1])
    for i = 1: length(Angle)
        plot(Angle(i,:), 'Linewidth', 1.5)
        hold on
    end
    ylabel('Angle (rad)')
    xlabel('Number of iteration')
    grid on
    set(gca,'XTick',(1:1:10))
    set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
    lgd = legend (str, 'NumColumns', 3);
    lgd.FontSize = 9;
    hold off   

    % % FDLF Voltage
    % figure('color', [1,1,1])
    % for i = 1: length(Volt_FD)
    %     plot(Volt_FD(i,:), 'Linewidth', 1.5)
    %     hold on
    % end
    % ylabel('Voltage (pu)')
    % xlabel('Number of iteration')
    % grid on
    % set(gca,'XTick',(1:1:10))
    % set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
    % lgd = legend (str, 'NumColumns', 4);
    % lgd.FontSize = 9;
    % hold off
    % 
    % % FDLF Angle
    % figure('color', [1,1,1])
    % for i = 1: length(Angle_FD)
    %     plot(Angle_FD(i,:), 'Linewidth', 1.5)
    %     hold on
    % %     ylim([0:1])
    % end
    % ylabel('Angle (rad)')
    % xlabel('Number of iteration')
    % grid on
    % set(gca,'XTick',(1:1:10))
    % set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
    % lgd = legend (str, 'NumColumns', 3);
    % lgd.FontSize = 9;
    % hold off
    % 
    % lgd = legend (str, 'NumColumns', 4);
    % % set(lgd,'position',poshL);      % Adjusting legend's position
    % % axis(hL,'off');                 % Turning its axis off
end