function obj = A_BH_ApplySpotMerging(obj)
mergeRadXY = 2;
mergeRadZ = 1;
rxysq = double(mergeRadXY) ^ 2;
rzsq = double(mergeRadZ) ^ 2;
mergeRad3 = sqrt(rxysq + rxysq + rzsq);
image_numbers = unique(obj.spotTable.image);
obj.spotTable{:, 'uid'} = zeros(size(obj.spotTable,1), 1);
obj.spotTable{:, 'likely_dup'} = uint32(zeros(size(obj.spotTable,1), 1)); %Reset all to 0
obj.spotTable{:, 'is_a_duplicate'} = uint32(zeros(size(obj.spotTable,1), 1)); %Reset all to 0
for im = 1:length(image_numbers)
    spotCount = sum(obj.spotTable.image==image_numbers(im));
    spots_in_image = find(obj.spotTable.image==image_numbers(im));
    obj.spotTable{spots_in_image, 'uid'} = [1:1:spotCount]';
    obj.spotTable{spots_in_image, 'likely_dup'} = uint32(zeros(spotCount, 1)); %Reset all to 0
    obj.spotTable{spots_in_image, 'is_a_duplicate'} = uint32(zeros(spotCount, 1)); %Reset all to 0
    %Matrix
    callmtx = NaN(spotCount, 3);
    callmtx(:,1) = obj.spotTable{spots_in_image, 'xinit'};
    callmtx(:,2) = obj.spotTable{spots_in_image, 'yinit'};
    callmtx(:,3) = obj.spotTable{spots_in_image, 'zinit'};
    cand_table = RNACoords.findMatchCandidates(callmtx, callmtx, mergeRad3, mergeRadZ);
    %Remove any self-matches
    %dist3 row col dist2 distz
    cand_table = cand_table((cand_table(:,1) > 0), :);
    if (~isempty(cand_table))
        %Find clusters
        gg = digraph(cand_table(:,2), cand_table(:,3));
        gbins = conncomp(gg);
        binnedSpots = size(gbins, 2);
        clusterCount = max(gbins, [], 'all', 'omitnan');
        for c = 1:clusterCount
            in_cluster = (gbins == c);
            nodeCount = nnz(in_cluster);
            if nodeCount > 1
                %find the biggest member to be "center" (using max
                %total intensity)
                totExpInt = immultiply(obj.spotTable{spots_in_image(1:binnedSpots), 'TotExpInt'}', in_cluster);
                [~, idx] = max(totExpInt, [], 'all', 'omitnan');
                nodeId = obj.spotTable{spots_in_image(idx), 'uid'};
                obj.spotTable{spots_in_image(in_cluster), 'likely_dup'} = nodeId;
                obj.spotTable{setdiff(spots_in_image(in_cluster),spots_in_image(idx)), 'is_a_duplicate'} = 1;
            end
        end
    end

end
end