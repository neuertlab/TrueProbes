<<<<<<< HEAD
<<<<<<< HEAD
function statsStruct = takeStats(data, mask, statsStruct, localXY, localZ)
    if nargin < 4; localXY = 0; end
    if nargin < 5; localZ = 0; end

%Apply mask if nonempty
    if ~isempty(mask) & (nnz(~mask) > 0)
        data = double(data);
        %fprintf('DEBUG -- Expected NaNs: %d\n', nnz(~mask));
        data(~mask) = NaN;
        %fprintf('DEBUG -- Added NaNs: %d\n', nnz(isnan(data)));
    end

    if nnz(isfinite(data)) < 1
        statsStruct.histo_y = [];
        statsStruct.histo_x = [];
        statsStruct.min = 0;
        statsStruct.max = 0;
        statsStruct.median = NaN;
        statsStruct.mad = NaN;
        statsStruct.mean = NaN;
        statsStruct.stdev = NaN;
        return;
    end

%Histo
    data_all = uint16(data(isfinite(data)));
    statsStruct.min = min(data_all, [], 'all', 'omitnan');
    statsStruct.max = max(data_all, [], 'all', 'omitnan');

    edges = [0:1:statsStruct.max];
    [bins, edges] = histcounts(data_all, edges);
    statsStruct.histo_y = uint32(bins);
    statsStruct.histo_x = uint16(edges(1:size(bins,2)));
    clear edges bins

    data_all = double(data_all);

    statsStruct.median = median(data_all, 'all', 'omitnan');
    statsStruct.mad = mad(data_all, 1, 'all');
    statsStruct.mean = mean(data_all, 'all', 'omitnan');
    statsStruct.stdev = std(data_all, 0, 'all', 'omitnan');
    clear data_all

%Local blocks, if applicable
    if (localXY > 0) & (localZ > 0)
        X = size(data, 2);
        Y = size(data, 1);
        Z = size(data, 3);
        zBlocks = ceil(Z ./ localZ);
        xBlocks = ceil(X ./ localXY);
        yBlocks = ceil(Y ./ localXY);
        totalBlocks = zBlocks * xBlocks * yBlocks;

        [fieldNames, fieldTypes] = getLocalRegionStatsFields();
        table_size = [totalBlocks size(fieldNames,2)];
        statsStruct.localBlocks = table('Size', table_size, 'VariableTypes',fieldTypes, 'VariableNames',fieldNames);

        i = 1;
        z0 = 1;
        for zb = 1:zBlocks
            z1 = min((z0 + localZ) - 1, Z);
            x0 = 1;
            for xb = 1:xBlocks
                x1 = min((x0 + localXY) - 1, X);
                y0 = 1;
                for yb = 1:yBlocks
                    y1 = min((y0 + localXY) - 1, Y);

                    %Do actual block.
                    statsStruct.localBlocks{i, 'x'} = x0;
                    statsStruct.localBlocks{i, 'y'} = y0;
                    statsStruct.localBlocks{i, 'z'} = z0;
                    statsStruct.localBlocks{i, 'width'} = x1 - x0 + 1;
                    statsStruct.localBlocks{i, 'height'} = y1 - y0 + 1;
                    statsStruct.localBlocks{i, 'depth'} = z1 - z0 + 1;

                    %Don't need mask again because it's already been
                    %applied.
                    localSStruct = ...
                        takeStats(data(y0:y1,x0:x1,z0:z1), [], struct(), 0, 0);

                    statsStruct.localBlocks{i, 'min'} = localSStruct.min;
                    statsStruct.localBlocks{i, 'max'} = localSStruct.max;
                    statsStruct.localBlocks{i, 'median'} = single(localSStruct.median);
                    statsStruct.localBlocks{i, 'mad'} = single(localSStruct.mad);
                    statsStruct.localBlocks{i, 'mean'} = single(localSStruct.mean);
                    statsStruct.localBlocks{i, 'stdev'} = single(localSStruct.stdev);

                    ccx = cell(1,1);
                    ccx{1} = localSStruct.histo_x;
                    ccy = cell(1,1);
                    ccy{1} = localSStruct.histo_y;
                    statsStruct.localBlocks{i, 'histo_x'} = ccx;
                    statsStruct.localBlocks{i, 'histo_y'} = ccy;

                    y0 = y1 + 1;
                    i = i + 1;
                end
                x0 = x1 + 1;
            end
            z0 = z1+1;
        end
    end

=======
function statsStruct = takeStats(data, mask, statsStruct, localXY, localZ)
    if nargin < 4; localXY = 0; end
    if nargin < 5; localZ = 0; end

%Apply mask if nonempty
    if ~isempty(mask) & (nnz(~mask) > 0)
        data = double(data);
        %fprintf('DEBUG -- Expected NaNs: %d\n', nnz(~mask));
        data(~mask) = NaN;
        %fprintf('DEBUG -- Added NaNs: %d\n', nnz(isnan(data)));
    end

    if nnz(isfinite(data)) < 1
        statsStruct.histo_y = [];
        statsStruct.histo_x = [];
        statsStruct.min = 0;
        statsStruct.max = 0;
        statsStruct.median = NaN;
        statsStruct.mad = NaN;
        statsStruct.mean = NaN;
        statsStruct.stdev = NaN;
        return;
    end

%Histo
    data_all = uint16(data(isfinite(data)));
    statsStruct.min = min(data_all, [], 'all', 'omitnan');
    statsStruct.max = max(data_all, [], 'all', 'omitnan');

    edges = [0:1:statsStruct.max];
    [bins, edges] = histcounts(data_all, edges);
    statsStruct.histo_y = uint32(bins);
    statsStruct.histo_x = uint16(edges(1:size(bins,2)));
    clear edges bins

    data_all = double(data_all);

    statsStruct.median = median(data_all, 'all', 'omitnan');
    statsStruct.mad = mad(data_all, 1, 'all');
    statsStruct.mean = mean(data_all, 'all', 'omitnan');
    statsStruct.stdev = std(data_all, 0, 'all', 'omitnan');
    clear data_all

%Local blocks, if applicable
    if (localXY > 0) & (localZ > 0)
        X = size(data, 2);
        Y = size(data, 1);
        Z = size(data, 3);
        zBlocks = ceil(Z ./ localZ);
        xBlocks = ceil(X ./ localXY);
        yBlocks = ceil(Y ./ localXY);
        totalBlocks = zBlocks * xBlocks * yBlocks;

        [fieldNames, fieldTypes] = getLocalRegionStatsFields();
        table_size = [totalBlocks size(fieldNames,2)];
        statsStruct.localBlocks = table('Size', table_size, 'VariableTypes',fieldTypes, 'VariableNames',fieldNames);

        i = 1;
        z0 = 1;
        for zb = 1:zBlocks
            z1 = min((z0 + localZ) - 1, Z);
            x0 = 1;
            for xb = 1:xBlocks
                x1 = min((x0 + localXY) - 1, X);
                y0 = 1;
                for yb = 1:yBlocks
                    y1 = min((y0 + localXY) - 1, Y);

                    %Do actual block.
                    statsStruct.localBlocks{i, 'x'} = x0;
                    statsStruct.localBlocks{i, 'y'} = y0;
                    statsStruct.localBlocks{i, 'z'} = z0;
                    statsStruct.localBlocks{i, 'width'} = x1 - x0 + 1;
                    statsStruct.localBlocks{i, 'height'} = y1 - y0 + 1;
                    statsStruct.localBlocks{i, 'depth'} = z1 - z0 + 1;

                    %Don't need mask again because it's already been
                    %applied.
                    localSStruct = ...
                        takeStats(data(y0:y1,x0:x1,z0:z1), [], struct(), 0, 0);

                    statsStruct.localBlocks{i, 'min'} = localSStruct.min;
                    statsStruct.localBlocks{i, 'max'} = localSStruct.max;
                    statsStruct.localBlocks{i, 'median'} = single(localSStruct.median);
                    statsStruct.localBlocks{i, 'mad'} = single(localSStruct.mad);
                    statsStruct.localBlocks{i, 'mean'} = single(localSStruct.mean);
                    statsStruct.localBlocks{i, 'stdev'} = single(localSStruct.stdev);

                    ccx = cell(1,1);
                    ccx{1} = localSStruct.histo_x;
                    ccy = cell(1,1);
                    ccy{1} = localSStruct.histo_y;
                    statsStruct.localBlocks{i, 'histo_x'} = ccx;
                    statsStruct.localBlocks{i, 'histo_y'} = ccy;

                    y0 = y1 + 1;
                    i = i + 1;
                end
                x0 = x1 + 1;
            end
            z0 = z1+1;
        end
    end

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function statsStruct = takeStats(data, mask, statsStruct, localXY, localZ)
    if nargin < 4; localXY = 0; end
    if nargin < 5; localZ = 0; end

%Apply mask if nonempty
    if ~isempty(mask) & (nnz(~mask) > 0)
        data = double(data);
        %fprintf('DEBUG -- Expected NaNs: %d\n', nnz(~mask));
        data(~mask) = NaN;
        %fprintf('DEBUG -- Added NaNs: %d\n', nnz(isnan(data)));
    end

    if nnz(isfinite(data)) < 1
        statsStruct.histo_y = [];
        statsStruct.histo_x = [];
        statsStruct.min = 0;
        statsStruct.max = 0;
        statsStruct.median = NaN;
        statsStruct.mad = NaN;
        statsStruct.mean = NaN;
        statsStruct.stdev = NaN;
        return;
    end

%Histo
    data_all = uint16(data(isfinite(data)));
    statsStruct.min = min(data_all, [], 'all', 'omitnan');
    statsStruct.max = max(data_all, [], 'all', 'omitnan');

    edges = [0:1:statsStruct.max];
    [bins, edges] = histcounts(data_all, edges);
    statsStruct.histo_y = uint32(bins);
    statsStruct.histo_x = uint16(edges(1:size(bins,2)));
    clear edges bins

    data_all = double(data_all);

    statsStruct.median = median(data_all, 'all', 'omitnan');
    statsStruct.mad = mad(data_all, 1, 'all');
    statsStruct.mean = mean(data_all, 'all', 'omitnan');
    statsStruct.stdev = std(data_all, 0, 'all', 'omitnan');
    clear data_all

%Local blocks, if applicable
    if (localXY > 0) & (localZ > 0)
        X = size(data, 2);
        Y = size(data, 1);
        Z = size(data, 3);
        zBlocks = ceil(Z ./ localZ);
        xBlocks = ceil(X ./ localXY);
        yBlocks = ceil(Y ./ localXY);
        totalBlocks = zBlocks * xBlocks * yBlocks;

        [fieldNames, fieldTypes] = getLocalRegionStatsFields();
        table_size = [totalBlocks size(fieldNames,2)];
        statsStruct.localBlocks = table('Size', table_size, 'VariableTypes',fieldTypes, 'VariableNames',fieldNames);

        i = 1;
        z0 = 1;
        for zb = 1:zBlocks
            z1 = min((z0 + localZ) - 1, Z);
            x0 = 1;
            for xb = 1:xBlocks
                x1 = min((x0 + localXY) - 1, X);
                y0 = 1;
                for yb = 1:yBlocks
                    y1 = min((y0 + localXY) - 1, Y);

                    %Do actual block.
                    statsStruct.localBlocks{i, 'x'} = x0;
                    statsStruct.localBlocks{i, 'y'} = y0;
                    statsStruct.localBlocks{i, 'z'} = z0;
                    statsStruct.localBlocks{i, 'width'} = x1 - x0 + 1;
                    statsStruct.localBlocks{i, 'height'} = y1 - y0 + 1;
                    statsStruct.localBlocks{i, 'depth'} = z1 - z0 + 1;

                    %Don't need mask again because it's already been
                    %applied.
                    localSStruct = ...
                        takeStats(data(y0:y1,x0:x1,z0:z1), [], struct(), 0, 0);

                    statsStruct.localBlocks{i, 'min'} = localSStruct.min;
                    statsStruct.localBlocks{i, 'max'} = localSStruct.max;
                    statsStruct.localBlocks{i, 'median'} = single(localSStruct.median);
                    statsStruct.localBlocks{i, 'mad'} = single(localSStruct.mad);
                    statsStruct.localBlocks{i, 'mean'} = single(localSStruct.mean);
                    statsStruct.localBlocks{i, 'stdev'} = single(localSStruct.stdev);

                    ccx = cell(1,1);
                    ccx{1} = localSStruct.histo_x;
                    ccy = cell(1,1);
                    ccy{1} = localSStruct.histo_y;
                    statsStruct.localBlocks{i, 'histo_x'} = ccx;
                    statsStruct.localBlocks{i, 'histo_y'} = ccy;

                    y0 = y1 + 1;
                    i = i + 1;
                end
                x0 = x1 + 1;
            end
            z0 = z1+1;
        end
    end

>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
end