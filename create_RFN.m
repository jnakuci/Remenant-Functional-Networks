function [fc_masked,bp_binary,bp_weighted]  = create_RFN(fc,bp,dn)

%% Input
% FC: Funcitonal connectivity 
% MC: Network used to mask the FC
% dn: desired density of the MC

%% Output: 
% FC_masked: FC matrix with biophyiscal network

%%
bp_weighted = threshold_proportional(bp,dn);
bp_binary = bp_weighted; 
bp_binary(bp_binary == 0) = 1;
bp_binary(bp_binary ~= 1) = 0; 
fc_masked = fc .* bp_binary;


end


