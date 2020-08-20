%GAP
%analyzes variance among genes within same-strain bacteria
%Noga Aharony, Kishony Lab, June 11, 2020


%% initialize
clear all; close all;
%initialize variables
initialize_gap

% 
% % run programs in commandline
% 
commandline_step

%get a pangneome file of the proteins
%get_protein_pangenome

% %% run functions on info
% 
create_alignment_biomaps

summarize_alignment_reports

process_plasfinder

process_resfinder

create_clusters

add_res_seq

% truncate_alignment_report
% 
% plot_clusters
% % % 
% observe_alignment
% % 
% 





