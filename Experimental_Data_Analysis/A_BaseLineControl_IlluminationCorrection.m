
rescale_func = @(x,minVal,maxVal) (x-minVal)/(maxVal-minVal); 
render_function = @(CM,R,G,B) cat(3,...
    CM(1,1)*R+CM(1,2)*G+CM(1,3)*B,...
    CM(2,1)*R+CM(2,2)*G+CM(2,3)*B,...
    CM(3,1)*R+CM(3,2)*G+CM(3,3)*B);
CM = 255*[0 0 0;1 0 0;0 0 0];
x_horizontal_center = @(Img) [1 size(Img,2)];
y_horizontal_center = @(Img) [size(Img,1)/2 size(Img,1)/2];
y_vertical_center = @(Img) [1 size(Img,2)];
x_vertical_center = @(Img) [size(Img,1)/2 size(Img,1)/2];
y_diagonal = @(Img) [1 size(Img,2)];
x_diagonal = @(Img) [1 size(Img,2)];
y_antidiagonal = @(Img)  [size(Img,2) 1];
x_antidiagonal = @(Img) [1 size(Img,2)];
profile_function_horizontal_center = @(Img) CATnWrapper(arrayfun(@(z) improfile(Img(:,:,z),x_horizontal_center(Img),y_horizontal_center(Img)),1:size(Img,3),'Un',0),2);
profile_function_vertical_center = @(Img) CATnWrapper(arrayfun(@(z) improfile(Img(:,:,z),x_vertical_center(Img),y_vertical_center(Img)),1:size(Img,3),'Un',0),2);
profile_function_diagonal = @(Img) CATnWrapper(arrayfun(@(z) improfile(Img(:,:,z),x_diagonal(Img),y_diagonal(Img)),1:size(Img,3),'Un',0),2);
profile_function_antidiagonal = @(Img) CATnWrapper(arrayfun(@(z) improfile(Img(:,:,z),x_antidiagonal(Img),y_antidiagonal(Img)),1:size(Img,3),'Un',0),2);

%PaintSHOP 01_16_22 has very large background cell, weird illumination
%Oligostan-HT 01_16_ Image 3, has very bright spots in nucleus but not
%labeled as nascent for some reason, vs nuclear well if all the spots are many nascent,
imgPath0 = 'JA_2025_02_18_CorrectionBaseline_Simulated_1_MMStack_Pos0.ome.tif';
[BaseLineControl_img_CY5,~] =  LoadTif(imgPath0,3,3,1);
Smoothed_BaseLineControl = zeros(size(BaseLineControl_img_CY5{3}));
BaseLineControl_ImageCorrectionFactor = zeros(size(BaseLineControl_img_CY5{3}));
for zed = 1:size(BaseLineControl_img_CY5{3},3)
        Smoothed_BaseLineControl(:,:,zed) = smoothdata2(BaseLineControl_img_CY5{3}(:,:,zed));
        BaseLineControl_ImageCorrectionFactor(:,:,zed) = double(Smoothed_BaseLineControl(:,:,zed))/double(max(Smoothed_BaseLineControl(:,:,zed),[],'all'));
end
BaseLineControl_ImageCorrectionFactor(:,:,zed) = double(Smoothed_BaseLineControl(:,:,zed))/double(max(Smoothed_BaseLineControl(:,:,zed),[],'all'));

plot(squeeze(max(BaseLineControl_img_CY5{3},[],[1 2])))
figure(200)
imagesc(BaseLineControl_img{3}(:,:,34),[190 655])
title('BaseLine Control Image','FontSize',15,'FontWeight','Bold','FontName','Apotos')
ax = gca;
ax.XTick = [];
ax.YTick = [];

plot(smoothdata2(profile_function_horizontal_center(BaseLineControl_img{3})))
plot(smoothdata2(profile_function_vertical_center(BaseLineControl_img{3})))
plot(smoothdata2(profile_function_diagonal(BaseLineControl_img{3})))
plot(smoothdata2(profile_function_antidiagonal(BaseLineControl_img{3})))





fn = 100;
figure(fn);set(gcf, 'Units','normalized','outerposition', [0 0 1 1]);
tFigReplicaImages{fn} = tiledlayout('horizontal');
tFigReplicaImages{fn}.TileSpacing = 'tight';
tFigReplicaImages{fn}.Padding = 'tight';
tFigReplicaImages{fn}.GridSize = [1 2];
ax = nexttile(tFigReplicaImages{fn},1);
imagesc(BaseLineControl_img{3}(:,:,34),[190 655]); hold on
colorbar
plot([1 2048],[2048 1],'r');
title('BaseLine Control Image','FontSize',15,'FontWeight','Bold','FontName','Apotos')
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax = nexttile(tFigReplicaImages{fn},2);
control_profile = profile_function_antidiagonal(BaseLineControl_img{3});
xlim([1 2048])
plot(1:2048,control_profile(:,34),'r');hold on
plot(1:2048,smoothdata(control_profile(:,34)),'k');hold on
xlabel('x position along anti-diagonal','FontWeight','Bold','FontSize',tFon,'FontName',fName);
ylabel('Intensity','FontWeight','Bold','FontSize',tFon,'FontName',fName);
title('BaseLine Control CY5 Intensity Profile','FontWeight','Bold','FontSize',tFon,'FontName',fName);
%plot(smoothdata2(profile_function_horizontal_center(BaseLineControl_img{3})))
legend({'raw trace','smoothed trace'},'FontSize',tFon)
%only change fit 



fn = 100;
figure(fn);set(gcf, 'Units','normalized','outerposition', [0 0 1 1]);
tFigReplicaImages{fn} = tiledlayout('horizontal');
tFigReplicaImages{fn}.TileSpacing = 'tight';
tFigReplicaImages{fn}.Padding = 'tight';
tFigReplicaImages{fn}.GridSize = [1 2];
ax = nexttile(tFigReplicaImages{fn},1);
imagesc(squeeze(chn{3,1}(:,:,34)),[min(squeeze(chn{3,1}(:,:,34)),[],'all') max(squeeze(chn{3,1}(:,:,34)),[],'all')])
title('With Dead Pixels','FontWeight','Bold','FontSize',tFon,'FontName',fName);
ax.XTick = [];
ax.YTick = [];
ax = nexttile(tFigReplicaImages{fn},2);
imagesc(squeeze(img_clean(:,:,34)),[min(squeeze(img_clean(:,:,34)),[],'all') max(squeeze(img_clean(:,:,34)),[],'all')])
title('Without Dead Pixels','FontWeight','Bold','FontSize',tFon,'FontName',fName);
%fit light source background drift profile 



[Xsurf,Ysurf] = meshgrid(1:2048,1:2048);
mesh(Xsurf,Ysurf,BaseLineControl_ImageCorrectionFactor(:,:,34),BaseLineControl_ImageCorrectionFactor(:,:,34))
 imagesc(BaseLineControl_ImageCorrectionFactor(:,:,34));colorbar;
 clim([0.75 1])
gaussThresh = 0.125;
X = 2048;Y=2048;Z=67;

sim_spot_func = @(spotTable,s) RNASpot.generateSimSpotFromFit_Table(spotTable, s, spotTable{s,'spotZFits'}, ceil(max(spotTable.xFWHM(s), spotTable.yFWHM(s)) * 2));
sim_mask_func = @(spotTable,s) double(sim_spot_func(spotTable,s) >= gaussThresh*max(sim_spot_func(spotTable,s), [], 'all', 'omitnan'));
in_spot_mask_func = @(spotTable,s) padarray(sim_mask_func(spotTable,s),...
    max(floor([0.5*(length(max(floor(spotTable.yabsloc(s)-spotTable.yFWHM(s)-5),1):min(ceil(spotTable.yabsloc(s)+spotTable.yFWHM(s)+5),2048))-2*ceil(double(max(spotTable.xFWHM(s), spotTable.yFWHM(s))) * 2)-1)...
             0.5*(length(max(floor(spotTable.xabsloc(s)-spotTable.xFWHM(s)-5),1):min(ceil(spotTable.xabsloc(s)+spotTable.xFWHM(s)+5),2048))-2*ceil(double(max(spotTable.xFWHM(s), spotTable.yFWHM(s))) * 2)-1) ...
             0.5*(length(max(floor(spotTable.zinit(s)-3),1):min(ceil(spotTable.zinit(s)+3),67))-length(spotTable{s,'spotZFits'}))]),0),0,'both');
in_spot_position_x_func = @(spotTable,s) min(max(ind2sub_Wrapper3D(in_spot_mask_func(spotTable,s),1).x+spotTable.xabsloc(s)-floor(0.5*(length(max(floor(spotTable.xabsloc(s)-spotTable.xFWHM(s)-5),1):min(ceil(spotTable.xabsloc(s)+spotTable.xFWHM(s)+5),2048))-1))-1,1),2048);
in_spot_position_y_func = @(spotTable,s) min(max(ind2sub_Wrapper3D(in_spot_mask_func(spotTable,s),1).y+spotTable.yabsloc(s)-floor(0.5*(length(max(floor(spotTable.yabsloc(s)-spotTable.yFWHM(s)-5),1):min(ceil(spotTable.yabsloc(s)+spotTable.yFWHM(s)+5),2048))-1))-1,1),2048);
in_spot_position_z_func = @(spotTable,s) min(max(ind2sub_Wrapper3D(in_spot_mask_func(spotTable,s),1).z+double(spotTable.zinit(s))-floor(0.5*(length(max(floor(spotTable.zinit(s)-3),1):min(ceil(spotTable.zinit(s)+3),67))-1))-1,1),67);
out_spot_position_x_func = @(spotTable,s) min(max(ind2sub_Wrapper3D(in_spot_mask_func(spotTable,s),0).x+spotTable.xabsloc(s)-floor(0.5*(length(max(floor(spotTable.xabsloc(s)-spotTable.xFWHM(s)-5),1):min(ceil(spotTable.xabsloc(s)+spotTable.xFWHM(s)+5),2048))-1))-1,1),2048);
out_spot_position_y_func = @(spotTable,s) min(max(ind2sub_Wrapper3D(in_spot_mask_func(spotTable,s),0).y+spotTable.yabsloc(s)-floor(0.5*(length(max(floor(spotTable.yabsloc(s)-spotTable.yFWHM(s)-5),1):min(ceil(spotTable.yabsloc(s)+spotTable.yFWHM(s)+5),2048))-1))-1,1),2048);
out_spot_position_z_func = @(spotTable,s) min(max(ind2sub_Wrapper3D(in_spot_mask_func(spotTable,s),0).z+double(spotTable.zinit(s))-floor(0.5*(length(max(floor(spotTable.zinit(s)-3),1):min(ceil(spotTable.zinit(s)+3),67))-1))-1,1),67);
bkgd_baseline_int_correction = @(spotTable,s) mean(BaseLineControl_ImageCorrectionFactor(sub2ind(size(BaseLineControl_ImageCorrectionFactor),out_spot_position_y_func(spotTable,s),out_spot_position_x_func(spotTable,s),out_spot_position_z_func(spotTable,s))));
spot_baseline_int_correction = @(spotTable,s) mean(BaseLineControl_ImageCorrectionFactor(sub2ind(size(BaseLineControl_ImageCorrectionFactor),in_spot_position_y_func(spotTable,s),in_spot_position_x_func(spotTable,s),in_spot_position_z_func(spotTable,s))));
sim_cellmask_func = @(imgTable,c) repmat(bwareafilt(imgTable.cell_mask{c},1),[1 1 min(imgTable.box{c}.z_top + 5,Z)-max(imgTable.box{c}.z_bottom - 5,1)+1]);
in_cell_mask_func = @(imgTable,c) padarray(sim_cellmask_func(imgTable,c),max( ...
                                                            floor(0.5*[min(imgTable.box{c}.bottom + 5,Y)-max(imgTable.box{c}.top - 5, 1)-min(imgTable.box{c}.bottom,Y)+max(imgTable.box{c}.top, 1) ...
                                                             min(imgTable.box{c}.right + 5,X)-max(imgTable.box{c}.left - 5,1)-min(imgTable.box{c}.right,X)+max(imgTable.box{c}.left,1) ... 
                                                             min(imgTable.box{c}.z_top + 5,Z)-max(imgTable.box{c}.z_bottom - 5,1)-min(imgTable.box{c}.z_top,Z)+max(imgTable.box{c}.z_bottom,1)]),0),0,'both');
in_cell_position_x_func = @(imgTable,c) min(max(ind2sub_Wrapper3D(in_cell_mask_func(imgTable,c),1).x+imgTable.box{c}.left-1,1),X);
in_cell_position_y_func = @(imgTable,c) min(max(ind2sub_Wrapper3D(in_cell_mask_func(imgTable,c),1).y+imgTable.box{c}.top-1,1),Y);
in_cell_position_z_func = @(imgTable,c) min(max(ind2sub_Wrapper3D(in_cell_mask_func(imgTable,c),1).z+imgTable.box{c}.z_bottom-1,1),Z);
out_cell_position_x_func = @(imgTable,c) min(max(ind2sub_Wrapper3D(in_cell_mask_func(imgTable,c),0).x+imgTable.box{c}.left-1,1),X);
out_cell_position_y_func = @(imgTable,c) min(max(ind2sub_Wrapper3D(in_cell_mask_func(imgTable,c),0).y+imgTable.box{c}.top-1,1),Y);
out_cell_position_z_func = @(imgTable,c) min(max(ind2sub_Wrapper3D(in_cell_mask_func(imgTable,c),0).z+imgTable.box{c}.z_bottom-1,1),Z);
imgbkgd_baseline_int_correction = @(imgTable,c) mean(BaseLineControl_ImageCorrectionFactor(sub2ind(size(BaseLineControl_ImageCorrectionFactor),out_cell_position_y_func(imgTable,c),out_cell_position_x_func(imgTable,c),out_cell_position_z_func(imgTable,c))));
cellbkgd_baseline_int_correction = @(imgTable,c) mean(BaseLineControl_ImageCorrectionFactor(sub2ind(size(BaseLineControl_ImageCorrectionFactor),in_cell_position_y_func(imgTable,c),in_cell_position_x_func(imgTable,c),in_cell_position_z_func(imgTable,c))));

