
function [perc_overlap_all_bp,bp_overlap] = percent_overlap(fc,sc,ec,gc,rc,gc_dist,rc_dist,dn)

v(:,1) = unwrapMat(threshold_proportional(sc, dn));
v(:,2) = unwrapMat(threshold_proportional(ec,dn));
v(:,3) = unwrapMat(threshold_proportional(gc,dn));
v(:,4) = unwrapMat(threshold_proportional(rc, dn));
v(:,5) = unwrapMat(threshold_proportional(gc_dist, dn));
v(:,6) = unwrapMat(threshold_proportional(rc_dist, dn));

% Determine the edges across biophysical networks
v(v~=0) = 1; 
v3 = sum(v,2);
v3(v3~=0) = 1;


% Edges in FC
v1(:,1) = unwrapMat(threshold_proportional(fc,dn));
v1(v1~=0) = 1;

%% Percent Overal Across All Biophysical Networks
perc_overlap_all_bp = 100 .* (1 - sum((v1-v3) > 0) ./sum(v1));


%% Percent Overal For Individual Biophysical Networks
idx1 = find(v1 == 1);

for j = 1:size(v,2)
    idx2 = find(v(:,j) == 1);
    bp_overlap(j) = 100.8*(sum(ismember(idx1,idx2))/numel(idx2));
end
