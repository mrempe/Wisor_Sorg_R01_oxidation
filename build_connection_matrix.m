function conn_mat = build_connection_matrix(N_PV,N_PC,N_SOM,PV_cell_index,PC_cell_index,SOM_cell_index,...
											PVPV_conn_prob,PVPC_conn_prob,PVSOM_conn_prob, ...
											PCPV_conn_prob,PCPC_conn_prob,PCSOM_conn_prob, ...
											SOMPV_conn_prob,SOMPC_conn_prob,SOMSOM_conn_prob)

%
% USAGE: conn_mat = build_connection_matrix(N_PV,N_PC,N_SOM,PV_cell_index,PC_cell_index,SOM_cell_index,...
%											PVPV_conn_prob,PVPC_conn_prob,PVSOM_conn_prob, ...
%											PCPV_conn_prob,PCPC_conn_prob,PCSOM_conn_prob, ...
%											SOMPV_conn_prob,SOMPC_conn_prob,SOMSOM_conn_prob)
%
% This function builds one large connection matrix for a network involving PV cells, PC cells, and SOM cells 
% based on the connection probabilities within each cell group and connections between cells in distinct groups. 
%
%
% INPUTS:  	N_PV:  # of Parvalbumin-positive cells
%			N_PC:  # of Pyramidal cells
%			N_SOM: # of Somatostatin cells
%			...

% the large matrix is built by concatenating 9 smaller matrices.  The colums of the matrix 
% refer to the presynaptic cells, and the rows refer to the post-synaptic cells

% --- PV->PV --------------------------------------------
conn_mat_PVPV  = rand(N_PV);
locs_low_PVPV  = find(conn_mat_PVPV<=PVPV_conn_prob);
locs_high_PVPV = find(conn_mat_PVPV>PVPV_conn_prob);
conn_mat_PVPV(locs_low_PVPV)  = 1;
conn_mat_PVPV(locs_high_PVPV) = 0;
% -------------------------------------------------------


% --- PV->PC --------------------------------------------
conn_mat_PVPC  = rand(N_PC,N_PV);
locs_low_PVPC  = find(conn_mat_PVPC<=PVPC_conn_prob);
locs_high_PVPC = find(conn_mat_PVPC>PVPC_conn_prob);
conn_mat_PVPC(locs_low_PVPC)  = 1;
conn_mat_PVPC(locs_high_PVPC) = 0;
% -------------------------------------------------------


% --- PV->SOM --------------------------------------------
conn_mat_PVSOM  = rand(N_SOM,N_PV);
locs_low_PVSOM  = find(conn_mat_PVSOM<=PVSOM_conn_prob);
locs_high_PVSOM = find(conn_mat_PVSOM>PVSOM_conn_prob);
conn_mat_PVSOM(locs_low_PVSOM)  = 1;
conn_mat_PVSOM(locs_high_PVSOM) = 0;
% -------------------------------------------------------


% --- PC->PV --------------------------------------------
conn_mat_PCPV  = rand(N_PV,N_PC);
locs_low_PCPV  = find(conn_mat_PCPV<=PCPV_conn_prob);
locs_high_PCPV = find(conn_mat_PCPV>PCPV_conn_prob);
conn_mat_PCPV(locs_low_PCPV)  = 1;
conn_mat_PCPV(locs_high_PCPV) = 0;
% -------------------------------------------------------


% --- PC->PC --------------------------------------------
conn_mat_PCPC  = rand(N_PC,N_PC);
locs_low_PCPC  = find(conn_mat_PCPC<=PCPC_conn_prob);
locs_high_PCPC = find(conn_mat_PCPC>PCPC_conn_prob);
conn_mat_PCPC(locs_low_PCPC)  = 1;
conn_mat_PCPC(locs_high_PCPC) = 0;
% -------------------------------------------------------

% --- PC->SOM --------------------------------------------
conn_mat_PCSOM  = rand(N_SOM,N_PC);
locs_low_PCSOM  = find(conn_mat_PCSOM<=PCSOM_conn_prob);
locs_high_PCSOM = find(conn_mat_PCSOM>PCSOM_conn_prob);
conn_mat_PCSOM(locs_low_PCSOM)  = 1;
conn_mat_PCSOM(locs_high_PCSOM) = 0;
% -------------------------------------------------------


% --- SOM->PV --------------------------------------------
conn_mat_SOMPV  = rand(N_PV,N_SOM);
locs_low_SOMPV  = find(conn_mat_SOMPV<=SOMPV_conn_prob);
locs_high_SOMPV = find(conn_mat_SOMPV>SOMPV_conn_prob);
conn_mat_SOMPV(locs_low_SOMPV)  = 1;
conn_mat_SOMPV(locs_high_SOMPV) = 0;
% -------------------------------------------------------


% --- SOM->PC --------------------------------------------
conn_mat_SOMPC  = rand(N_PC,N_SOM);
locs_low_SOMPC  = find(conn_mat_SOMPC<=SOMPC_conn_prob);
locs_high_SOMPC = find(conn_mat_SOMPC>SOMPC_conn_prob);
conn_mat_SOMPC(locs_low_SOMPC)  = 1;
conn_mat_SOMPC(locs_high_SOMPC) = 0;
% -------------------------------------------------------


% --- SOM->SOM --------------------------------------------
conn_mat_SOMSOM  = rand(N_SOM,N_SOM);
locs_low_SOMSOM  = find(conn_mat_SOMSOM<=SOMSOM_conn_prob);
locs_high_SOMSOM = find(conn_mat_SOMSOM>SOMSOM_conn_prob);
conn_mat_SOMSOM(locs_low_SOMSOM)  = 1;
conn_mat_SOMSOM(locs_high_SOMSOM) = 0;
% -------------------------------------------------------

conn_mat = [conn_mat_PVPV  conn_mat_PCPV  conn_mat_SOMPV; ...
			conn_mat_PVPC  conn_mat_PCPC  conn_mat_SOMPC; ...
			conn_mat_PVSOM conn_mat_PCSOM conn_mat_SOMSOM];