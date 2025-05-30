function is_main_spot = A_BHJH_ApplySpotMergingGeneric(spotTable,xvar,yvar,zvar,int_var)
mergeRadXY = 2;
mergeRadZ = 1;
rxysq = double(mergeRadXY) ^ 2;
rzsq = double(mergeRadZ) ^ 2;
mergeRad3 = sqrt(rxysq + rxysq + rzsq);
image_numbers = unique(spotTable.image);
spotTable{:, 'uid'} = zeros(size(spotTable,1), 1);
spotTable{:, 'likely_dup'} = uint32(zeros(size(spotTable,1), 1)); %Reset all to 0
spotTable{:, 'is_a_duplicate'} = uint32(zeros(size(spotTable,1), 1)); %Reset all to 0
for im = 1:length(image_numbers)
    cells_with_spots = setdiff(unique(spotTable.cell(find(double(spotTable.image==image_numbers(im))))),0);
    for ic = 1:length(cells_with_spots)
        spotCount = sum(double(spotTable.image==image_numbers(im)).*double(spotTable.cell==cells_with_spots(ic)));
        spots_in_image_and_cell = find(double(spotTable.image==image_numbers(im)).*double(spotTable.cell==cells_with_spots(ic)));
        spotTable{spots_in_image_and_cell, 'uid'} = [1:1:spotCount]';
        spotTable{spots_in_image_and_cell, 'likely_dup'} = uint32(zeros(spotCount, 1)); %Reset all to 0
        spotTable{spots_in_image_and_cell, 'is_a_duplicate'} = uint32(zeros(spotCount, 1)); %Reset all to 0
        %Matrix
        callmtx = NaN(spotCount, 3);
        callmtx(:,1) = spotTable{spots_in_image_and_cell, xvar};
        callmtx(:,2) = spotTable{spots_in_image_and_cell, yvar};
        callmtx(:,3) = spotTable{spots_in_image_and_cell, zvar};
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
                    totExpInt = immultiply(spotTable{spots_in_image_and_cell(1:binnedSpots), int_var}', in_cluster);
                    [~, idx] = max(totExpInt, [], 'all', 'omitnan');
                    nodeId = spotTable{spots_in_image_and_cell(idx), 'uid'};
                    spotTable{spots_in_image_and_cell(in_cluster), 'likely_dup'} = nodeId;
                    spotTable{setdiff(spots_in_image_and_cell(in_cluster),spots_in_image_and_cell(idx)), 'is_a_duplicate'} = 1;
                end
            end
        end
    end
end
is_main_spot = double(~spotTable.is_a_duplicate);
end