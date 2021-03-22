% input for proxy identification
ident_input_proxy    = {'proxy',VAR_Est, DATASET_VAR,{'DGS1'},3,[1979 7],[2007 1]};

% input for sign restriction identification
ident_input_sign     = {'sign', Para, VAR_Ident_short...
    ,[ 1           4           1          %  VAR1
       1           4           1          %  VAR2
       0           0           0]};

% input for max share identification
ident_input_maxshare     = {'maxshare', VAR_Ident_short, 1, 8, 1};