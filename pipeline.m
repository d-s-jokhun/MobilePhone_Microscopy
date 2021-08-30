

% pix_size=0.043;  %in um
pix_size=0.215;  %in um


% sample = 'M';
% filenames = dir (['*',sample,'*','.tif']);
filenames = dir (['*.tif']);



'loading images'
XYZ={};
parfor reader_count=1:size(filenames,1);
    reader_count
    filename = filenames(reader_count).name;
    
    Reader = bfGetReader (filename);
    OmeMeta = Reader.getMetadataStore();
    Num_of_Pixels_Z = OmeMeta.getPixelsSizeZ(0).getValue();
    iSeries =1;
    Reader.setSeries(iSeries - 1);
    iT=1;
    iCh=1
    XYZ_temp =uint16([]);
    for iZ=1:Num_of_Pixels_Z;
        iPlane = Reader.getIndex(iZ-1, iCh-1, iT-1) + 1;     %%% The last '1-1' is for timepoint 0 (the 1st timepoint)
        XYZ_temp(:,:,iZ)= bfGetPlane(Reader, iPlane);
    end
    XYZ{1,reader_count}=XYZ_temp;
    
    
end
'Files imported'


%%


%%
'Manual Selection'

Selected_proper_2d_Raw=cell(size(XYZ));
for file_count=1:size(XYZ,2)
    ['file ',num2str(file_count),' of ',num2str(size(XYZ,2))]
    [Selected_proper_2d_Raw{1,file_count}]= manual_sel(XYZ{1,file_count});
end
Selected_proper_2d_Raw=Selected_proper_2d_Raw(~cellfun('isempty',Selected_proper_2d_Raw));

%%




%%
tic
'Basic Measurements'


Selected_proper_2d_Raw=Selected_Mobile_BJ_23g;
MetaData.Filename='Mobile BJ_23g';
MetaData.Voxel_Size_X=pix_size;
MetaData.Voxel_Size_Y=pix_size;

Basic_Measurements=cell(size(Selected_proper_2d_Raw));
Boundary_RadOfCurva_K=cell(size(Selected_proper_2d_Raw));
Norm_Pix_Intensities=cell(size(Selected_proper_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_proper_2d_Raw));
num_of_nuc=zeros(size(Selected_proper_2d_Raw));
total_num_of_nuc=zeros(size(Selected_proper_2d_Raw,1),1);


for condition_count=1:size(Selected_proper_2d_Raw,1)
    parfor file_count=1:size(Selected_proper_2d_Raw,2)
        if size(Selected_proper_2d_Raw{condition_count,file_count},1)*size(Selected_proper_2d_Raw{condition_count,file_count},2)>0
            
            
            
            ['File ',num2str(file_count),' of ',num2str(size(Selected_proper_2d_Raw,2)),' from condition ',num2str(condition_count),' of ',num2str(size(Selected_proper_2d_Raw,1))]
            
            [Basic_Measurements{condition_count,file_count},Boundary_RadOfCurva_K{condition_count,file_count},Norm_Pix_Intensities{condition_count,file_count}]=basic_mea(Selected_proper_2d_Raw{condition_count,file_count},MetaData);
            
        end
    end
    
    %Reorganizing Boundary_RadOfCurva_K into Result_Boundary_RadOfCurva_K
    %adding nuc_label to the Result_Basic_Meansurements sheet
    %adding nuc_label to the Result_Norm_Pix_Intensities and reorganizing it.
    'Reorganizing output'
    
    for file_count1=1:size(Selected_proper_2d_Raw,2)
        number=size(Norm_Pix_Intensities{condition_count,file_count1},1);
        if number~=0
            number=number-2;
        end
        num_of_pix_in_largest_nuc(condition_count,file_count1)=number;
        num_of_nuc(condition_count,file_count1)=size(Norm_Pix_Intensities{condition_count,file_count1},2);
    end
    
    total_num_of_nuc(condition_count,1)=sum(num_of_nuc(condition_count,:));
    
    Result_Norm_Pix_Intensities=cell(max(num_of_pix_in_largest_nuc(condition_count,:)),total_num_of_nuc(condition_count,1));
    Result_Boundary_RadOfCurva_K=cell(size(Selected_proper_2d_Raw,1),1);
    lower_bound=1;
    for file_count=1:size(Selected_proper_2d_Raw,2)
        if size(Selected_proper_2d_Raw{condition_count,file_count},1)*size(Selected_proper_2d_Raw{condition_count,file_count},2)>0
            Result_Boundary_RadOfCurva_K{condition_count,1}(1:size(Boundary_RadOfCurva_K{condition_count,file_count},1),file_count)=Boundary_RadOfCurva_K{condition_count,file_count}(:,1);
            for nuc_count=1:size (Selected_proper_2d_Raw{condition_count,file_count},3)
                Norm_Pix_Intensities{condition_count,file_count}(2,nuc_count)=num2cell(1);%num2cell(mean(nonzeros(Selected_proper_2d_AutoLabelled{condition_count,file_count}(:,:,nuc_count))));
                Basic_Measurements{condition_count,file_count}{nuc_count,2}=1;%mean(nonzeros(Selected_proper_2d_AutoLabelled{condition_count,file_count}(:,:,nuc_count)));
            end
            
            upper_bound=lower_bound-1+size(Norm_Pix_Intensities{condition_count,file_count},2);
            Result_Norm_Pix_Intensities(1:size(Norm_Pix_Intensities{condition_count,file_count},1),lower_bound:upper_bound)=Norm_Pix_Intensities{condition_count,file_count};
            lower_bound=upper_bound+1;
        end
    end
    
    
end

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Prob of Int Distri being Normal' 'Entropy'...
    'Relative concavity' 'Average curvature x Perimeter' 'SD in curvature(um^-1)' 'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter' 'Num. of times curvature changes polarity'...
    'Max positive curvature (um^-1)' 'Num of peaks with curvature >= 80% of Max positive' 'Avg prominance of peaks with curvature >= 80% of Max positive(um)' 'Avg width of peaks with curvature >= 80% of Max positive(um)'...
    'Max negative curvature (um^-1)' 'Num of peaks with curvature >= 80% of Max negative' 'Avg prominance of peaks with curvature >= 80% of Max negative(um)' 'Avg width of peaks with curvature >= 80% of Max negative(um)'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements{:});


Result_Combined_Norm_Pix_Intensities={};
Result_Combined_Norm_Pix_Intensities=Result_Norm_Pix_Intensities(3:end,:);
Result_Combined_Norm_Pix_Intensities=vertcat(Result_Combined_Norm_Pix_Intensities{:});

Result_Basic_Measurements=Result_Basic_Measurements(:,1:18);


clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities...
    Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number Result_Boundary_RadOfCurva_K...
    Result_Norm_Pix_Intensities Result_Combined_Norm_Pix_Intensities
toc










'Spatial Correlation'
tic

cell_num=0;

result_norm_corr_maps_Cir={};
result_norm_corr_maps_Rec={};
result_MeanCorr_vs_r={};
result_VarInCorr_vs_r={};
result_MaxDiffInCorr_vs_r={};
result_LenghtAt0Corr=[];
Corr_along_X={};
Corr_along_Y={};
result_lengthscale_along_X=[];
result_lengthscale_along_Y=[];



for file_count=1:size(Selected_proper_2d_Raw,2);
    file_count
    
    excluded=[];  % files to exclude.  e.g, if you want to exclude file_count 7,10 and 12, write [7 10 12]. Othervise write [].
    condition = file_count~=excluded;
    
    if prod(condition)==1
        
        raw=sum(Selected_proper_2d_Raw{1,file_count},3);
        nuc_bw=raw>0;
        
        stats_orient=regionprops(nuc_bw,'PixelIdxList','Area','Orientation');
        if size(stats_orient,1)>1
            area_check=[];
            template=zeros(size(nuc_bw));
            for area_count=1:size(stats_orient,1)
                area_check(area_count)=stats_orient(area_count).Area;
            end
            [~,id]=max(area_check);
            template(stats_orient(id).PixelIdxList)=1;
            raw=raw.*double(template);
            nuc_bw=raw>0;
            stats_orient=regionprops(nuc_bw,'PixelIdxList','Area','Orientation');
        end
        
        
        nuc=imrotate(raw,-stats_orient.Orientation,'bilinear');
        nuc_bw=nuc>0;
        
        %                 nuc_bw=imclose(nuc_bw,strel('disk',1));
        %                 nuc_edge=imdilate(edge(nuc_bw),strel('disk',1));
        %                        imtool(sum(XYZ{1,1},3)+(nuc_edge*max(max(sum(XYZ{1,1},3)))),[]);
        
        stats_select = regionprops(nuc_bw,'Area','PixelIdxList');  %to check the num of objects detected
        if size(stats_select,1)>1
            ['More than 1 object in file ',num2str(file_count)]
            imtool(nuc_bw,[])
            [~,max_area_idx] = max([stats_select.Area]);
            
            nuc_bw=zeros(size(nuc_bw));
            nuc_bw(stats_select(max_area_idx).PixelIdxList)=1;
            nuc=nuc.*double(nuc_bw);
        end
        
        
        cell_num=cell_num+1;
        %                     cell_num
        
        stats = regionprops(nuc_bw,'Centroid','ConvexHull','BoundingBox');
        
        convex_row_coor=stats.ConvexHull(:,2);
        convex_col_coor=stats.ConvexHull(:,1);
        convex_row_dist_frm_cen=convex_row_coor-stats.Centroid(1,2);
        convex_col_dist_frm_cen=convex_col_coor-stats.Centroid(1,1);
        convex_dist_frm_cen= sqrt((convex_row_dist_frm_cen.^2)+(convex_col_dist_frm_cen.^2));
        radius=floor(min(convex_dist_frm_cen))-5; % -10 because the edge is usally just scattered light from the actual nucleus
        
        %Creating a circle inside the nucleus
        rows = size(nuc_bw,1);  % circle will be in a matrix of same size as the original image
        cols = size(nuc_bw,2);
        center = stats.Centroid;  % circle will have the same centroid as the original image
        [xMat,yMat] = meshgrid(1:cols,1:rows);
        distFromCenter = sqrt((xMat-center(1)).^2 + (yMat-center(2)).^2);
        CirMat = distFromCenter<=radius;
        
        cir_edge=imdilate(edge(CirMat),strel('disk',1));
        %                                   imtool(nuc+(cir_edge*max(max(nuc))),[]);
        
        cropped_cir = double(nuc) .* double(CirMat);
        %                     imtool(cropped_cir,[])
        
        RecMat=zeros(size(nuc));
        RecMat(ceil(stats.BoundingBox(2))+1:(floor(stats.BoundingBox(2)))+stats.BoundingBox(4)-1,ceil(stats.BoundingBox(1))+1:(floor(stats.BoundingBox(1)))+stats.BoundingBox(3)-1)=1;
        for count_erode=1:stats.BoundingBox(4)/2
            if sum(sum(RecMat))>sum(sum(RecMat.*nuc_bw))
                RecMat=imerode(RecMat,strel('square',3));
            else
                break
            end
        end
        RecMat=imerode(RecMat,strel('square',9));
        
        rec_edge=imdilate(edge(RecMat),strel('disk',1));
        %                                   imtool(nuc+(rec_edge*max(max(nuc))),[]);
        
        cropped_rec = double(nuc) .* double(RecMat);
        %                     imtool(cropped_rec,[])
        
        %         cropped_cir = histeq(cropped_cir);
        %         cropped_cir = adapthisteq(cropped_cir);
        mean_int_cir = mean (nonzeros(cropped_cir));
        pre_CirImg_0_mean=cropped_cir - mean_int_cir;
        pre_CirImg_0_mean=pre_CirImg_0_mean/max(max(pre_CirImg_0_mean));
        CirImg_0_mean = pre_CirImg_0_mean.*double(CirMat);
        CirImg_to_be_analysed = imcrop(CirImg_0_mean,[(stats.Centroid(1,1)-radius+1) (stats.Centroid(1,2)-radius+1) (2*radius)-1 (2*radius)-1]);   %[xmin ymin width height]  %leaving 5 empty pixels on each side
        
        %         cropped_rec = histeq(cropped_rec);
        %         cropped_rec = adapthisteq(cropped_rec);
        mean_int_rec = mean (nonzeros(cropped_rec));
        pre_RecImg_0_mean=cropped_rec - mean_int_rec;
        pre_RecImg_0_mean=pre_RecImg_0_mean/max(max(pre_RecImg_0_mean));
        RecImg_0_mean = pre_RecImg_0_mean.*double(RecMat);
        RecCrop_stats=regionprops(RecMat,'BoundingBox');
        RecImg_to_be_analysed = imcrop(RecImg_0_mean,RecCrop_stats.BoundingBox);   %[xmin ymin width height]  %leaving 5 empty pixels on each side
        
        
        norm2DXCorrCir=normxcorr2(CirImg_to_be_analysed,CirImg_to_be_analysed);
        norm2DXCorrRec=normxcorr2(RecImg_to_be_analysed,RecImg_to_be_analysed);
        
        
        %Creating a circle to remove abberant values from the map(correlation for positions without overlap of the cropped circular nucleus
        rowsCorrMap = size(norm2DXCorrCir,1);  % circle will be in a matrix of same size as the original image
        colsCorrMap = size(norm2DXCorrCir,2);
        centerCorrMap = [((colsCorrMap-1)/2)+1 ((rowsCorrMap-1)/2)+1];
        [xMatCorrMap,yMatCorrMap] = meshgrid(1:colsCorrMap,1:rowsCorrMap);
        distFromCenterCorrMap = sqrt((xMatCorrMap-centerCorrMap(1)).^2 + (yMatCorrMap-centerCorrMap(2)).^2);
        circleMatCorrMap = distFromCenterCorrMap<=((radius*2)-1);
        norm_corr_map_Cir_n= norm2DXCorrCir.*double(circleMatCorrMap);
        
        norm_corr_map_Rec_n = norm2DXCorrRec;
        
        result_norm_corr_maps_Cir{cell_num}=norm_corr_map_Cir_n;
        result_norm_corr_maps_Rec{cell_num}=norm_corr_map_Rec_n;
        
        result_MeanCorr_vs_r {1,(cell_num*2)-1}=0;
        result_MeanCorr_vs_r {1,(cell_num*2)}=1;
        result_VarInCorr_vs_r {1,(cell_num*2)-1}=0;
        result_VarInCorr_vs_r {1,(cell_num*2)}=0;
        result_MaxDiffInCorr_vs_r {1,(cell_num*2)-1}=0;
        result_MaxDiffInCorr_vs_r {1,(cell_num*2)}=0;
        
        Corr_along_X {1,cell_num}(1:size(norm_corr_map_Rec_n,2),1)=(-(size(norm_corr_map_Rec_n,2)-1)/2:(size(norm_corr_map_Rec_n,2)-1)/2)*pix_size;
        Corr_along_X {1,cell_num}(1:size(norm_corr_map_Rec_n,2),2)=norm_corr_map_Rec_n(((size(norm_corr_map_Rec_n,1)-1)/2)+1,:);
        Corr_along_Y {1,cell_num}(1:size(norm_corr_map_Rec_n,1),1)=(-(size(norm_corr_map_Rec_n,1)-1)/2:(size(norm_corr_map_Rec_n,1)-1)/2)*pix_size;
        Corr_along_Y {1,cell_num}(1:size(norm_corr_map_Rec_n,1),2)=norm_corr_map_Rec_n(:,((size(norm_corr_map_Rec_n,2)-1)/2)+1);
        
        
        %%finding frequency along X
        neg_template_max=Corr_along_X {1,cell_num}(:,1)<-0.25;
        neg_template_min=Corr_along_X {1,cell_num}(:,1)>(min(Corr_along_X {1,cell_num}(:,1))+0.5);
        neg_template=neg_template_max.*neg_template_min;
        pos_template_max=Corr_along_X {1,cell_num}(:,1)<(max(Corr_along_X {1,cell_num}(:,1))-0.5);
        pos_template_min=Corr_along_X {1,cell_num}(:,1)>0.25;
        pos_template=pos_template_max.*pos_template_min;
        
        neg_segmentX=[];
        pos_segmentX=[];
        for count_selection=1:size(Corr_along_X{1,cell_num},1)
            if neg_template(count_selection,1)==1
                neg_segmentX(end+1,1)=Corr_along_X{1,cell_num}(count_selection,1);
                neg_segmentX(end,2)=Corr_along_X{1,cell_num}(count_selection,2);
            end
            if pos_template(count_selection,1)==1
                pos_segmentX(end+1,1)=Corr_along_X{1,cell_num}(count_selection,1);
                pos_segmentX(end,2)=Corr_along_X{1,cell_num}(count_selection,2);
            end
        end
        
        fitobject_negX = fit(neg_segmentX(:,1),neg_segmentX(:,2),'smoothingspline','SmoothingParam',0.1);
        fitobject_posX = fit(pos_segmentX(:,1),pos_segmentX(:,2),'smoothingspline','SmoothingParam',0.1);
        
        %                     figure('Name','neg segment of autocorr along X')
        %                     plot(fitobject_negX,neg_segmentX(:,1),neg_segmentX(:,2))
        %                     figure('Name','pos segment of autocorr along X')
        %                     plot(fitobject_posX,pos_segmentX(:,1),pos_segmentX(:,2))
        
        residualX_NegativeSide=neg_segmentX(:,2)-fitobject_negX(neg_segmentX(:,1));
        residualX_NegativeSide=residualX_NegativeSide-mean(residualX_NegativeSide);
        residualX_PositiveSide=pos_segmentX(:,2)-fitobject_posX(pos_segmentX(:,1));
        residualX_PositiveSide=residualX_PositiveSide-mean(residualX_PositiveSide);
        autocorrX_neg=xcorr(residualX_NegativeSide,residualX_NegativeSide);
        autocorrX_pos=xcorr(residualX_PositiveSide,residualX_PositiveSide);
        
        [~,locsX]=findpeaks(-autocorrX_neg(((length(autocorrX_neg)-1)/2)+2:end));
        lengthscale_negX=2*locsX(1)*pix_size;
        result_lengthscale_along_X(cell_num,1)=lengthscale_negX;
        [~,locsX]=findpeaks(-autocorrX_pos(((length(autocorrX_pos)-1)/2)+2:end));
        lengthscale_posX=2*locsX(1)*pix_size;
        result_lengthscale_along_X(cell_num,2)=lengthscale_posX;
        result_lengthscale_along_X(cell_num,3)=(lengthscale_negX+lengthscale_posX)/2;
        
        
        %%finding frequency along Y
        neg_template_max=Corr_along_Y {1,cell_num}(:,1)<-0.25;
        neg_template_min=Corr_along_Y {1,cell_num}(:,1)>(min(Corr_along_Y {1,cell_num}(:,1))+0.5);
        neg_template=neg_template_max.*neg_template_min;
        pos_template_max=Corr_along_Y {1,cell_num}(:,1)<(max(Corr_along_Y {1,cell_num}(:,1))-0.5);
        pos_template_min=Corr_along_Y {1,cell_num}(:,1)>0.25;
        pos_template=pos_template_max.*pos_template_min;
        
        neg_segmentY=[];
        pos_segmentY=[];
        for count_selection=1:size(Corr_along_Y{1,cell_num},1)
            if neg_template(count_selection,1)==1
                neg_segmentY(end+1,1)=Corr_along_Y{1,cell_num}(count_selection,1);
                neg_segmentY(end,2)=Corr_along_Y{1,cell_num}(count_selection,2);
            end
            if pos_template(count_selection,1)==1
                pos_segmentY(end+1,1)=Corr_along_Y{1,cell_num}(count_selection,1);
                pos_segmentY(end,2)=Corr_along_Y{1,cell_num}(count_selection,2);
            end
        end
        
        fitobject_negY = fit(neg_segmentY(:,1),neg_segmentY(:,2),'smoothingspline','SmoothingParam',0.1);
        fitobject_posY = fit(pos_segmentY(:,1),pos_segmentY(:,2),'smoothingspline','SmoothingParam',0.1);
        
        %                     figure('Name','neg segment of autocorr along Y')
        %                     plot(fitobject_negY,neg_segmentY(:,1),neg_segmentY(:,2))
        %                     figure('Name','pos segment of autocorr along Y')
        %                     plot(fitobject_posY,pos_segmentY(:,1),pos_segmentY(:,2))
        
        residualY_NegativeSide=neg_segmentY(:,2)-fitobject_negY(neg_segmentY(:,1));
        residualY_NegativeSide=residualY_NegativeSide-mean(residualY_NegativeSide);
        residualY_PositiveSide=pos_segmentY(:,2)-fitobject_posY(pos_segmentY(:,1));
        residualY_PositiveSide=residualY_PositiveSide-mean(residualY_PositiveSide);
        autocorrY_neg=xcorr(residualY_NegativeSide,residualY_NegativeSide);
        autocorrY_pos=xcorr(residualY_PositiveSide,residualY_PositiveSide);
        
        [~,locsY]=findpeaks(-autocorrY_neg(((length(autocorrY_neg)-1)/2)+2:end));
        lengthscale_negY=2*locsY(1)*pix_size;
        result_lengthscale_along_Y(cell_num,1)=lengthscale_negY;
        [~,locsY]=findpeaks(-autocorrY_pos(((length(autocorrY_pos)-1)/2)+2:end));
        lengthscale_posY=2*locsY(1)*pix_size;
        result_lengthscale_along_Y(cell_num,2)=lengthscale_posY;
        result_lengthscale_along_Y(cell_num,3)=(lengthscale_negY+lengthscale_posY)/2;
        
        
        
        sign_check=0; %used to determine when the mean correlation for each cell reaches 0
        
        for r=1:(size(norm_corr_map_Cir_n,2)-1)/2
            %Creating a circle to find the mean correlation along the edge of incremental radii
            rowsCorrMapR = size(norm_corr_map_Cir_n,1);  % circle will be in a matrix of same size as the original image
            colsCorrMapR = size(norm_corr_map_Cir_n,2);
            centerCorrMapR = [((colsCorrMapR-1)/2)+1 ((rowsCorrMapR-1)/2)+1];
            [xMatCorrMapR,yMatCorrMapR] = meshgrid(1:colsCorrMapR,1:rowsCorrMapR);
            distFromCenterCorrMapR = sqrt((xMatCorrMapR-centerCorrMapR(1)).^2 + (yMatCorrMapR-centerCorrMapR(2)).^2);
            circleMatCorrMapR = distFromCenterCorrMapR<=r;
            cir_edge=edge(circleMatCorrMapR);
            
            CorrMap_at_r = norm_corr_map_Cir_n.*double(cir_edge);
            CorrVals_at_r = nonzeros(CorrMap_at_r);  %takes only values from the edge but misses the zeros from that edge
            if size(CorrVals_at_r,1) < sum(sum(cir_edge))
                CorrVals_at_r (sum(sum(cir_edge)),1)=0;  %adds the missing zeros to the list
            end
            
            result_MeanCorr_vs_r {r+1,(cell_num*2)-1}=r*pix_size;
            result_MeanCorr_vs_r {r+1,(cell_num*2)}=mean(CorrVals_at_r);
            result_VarInCorr_vs_r {r+1,(cell_num*2)-1}=r*pix_size;
            result_VarInCorr_vs_r {r+1,(cell_num*2)}=var(CorrVals_at_r);
            result_MaxDiffInCorr_vs_r {r+1,(cell_num*2)-1}=r*pix_size;
            result_MaxDiffInCorr_vs_r {r+1,(cell_num*2)}=max(CorrVals_at_r)-min(CorrVals_at_r);
            
            %getting the 0 correlation lengthscale for each cell
            if sign_check==0
                if mean(CorrVals_at_r)<0
                    sign_check=1;
                    result_LenghtAt0Corr(cell_num,1)=r*pix_size;
                end
            end
        end
        
        
        %                     imtool(sum(XYZ{1,file_count},3),[])
        
    end
    
end


result_Corr_along_X={};
result_Corr_along_Y={};
for cell_count=1:size(Corr_along_X,2)
    for row_count1=1:size(Corr_along_X{1,cell_count},1)
        result_Corr_along_X{row_count1,(cell_count*2)-1}=Corr_along_X{1,cell_count}(row_count1,1);
        result_Corr_along_X{row_count1,(cell_count*2)}=Corr_along_X{1,cell_count}(row_count1,2);
    end
    for row_count2=1:size(Corr_along_Y{1,cell_count},1)
        result_Corr_along_Y{row_count2,(cell_count*2)-1}=Corr_along_Y{1,cell_count}(row_count2,1);
        result_Corr_along_Y{row_count2,(cell_count*2)}=Corr_along_Y{1,cell_count}(row_count2,2);
    end
end

result_lengthscale_along_X_by_lengthscale_along_Y=result_lengthscale_along_X./result_lengthscale_along_Y;



correlation_summary={'LenghtAt0Corr','lengthscale_along_X','lengthscale_along_Y'};
for count=1:size(result_LenghtAt0Corr,1)
    correlation_summary{count+1,1}=result_LenghtAt0Corr(count,1);
    correlation_summary{count+1,2}=result_lengthscale_along_X(count,3);
    correlation_summary{count+1,3}=result_lengthscale_along_Y(count,3);
end

Result_Basic_Measurements=horzcat(Result_Basic_Measurements,correlation_summary);

Basic_Measurements=cell2mat(Result_Basic_Measurements(2:end,5:end));


clearvars -except *Basic_Measurements* MetaData pix_size result_Corr_along_X result_Corr_along_Y result_LenghtAt0Corr result_lengthscale_along_X result_lengthscale_along_Y result_MeanCorr_vs_r result_norm_corr_maps_Cir result_norm_corr_maps_Rec Selected*






'Spatial correlation complete'





% figure ('Name','Cir Corr Map')
% colormap jet
% surf(result_norm_corr_maps_Cir{cell_to_disp})
% figure ('Name','Rec Corr Map')
% colormap jet
% surf(result_norm_corr_maps_Rec{cell_to_disp})


% %% finding mean from all cells
% for count_r=1:size(result_MeanCorr_vs_r,1)
%     list_MeanCorr_vs_r=[];
%     for count_n=1:size(result_MeanCorr_vs_r,2)/2
%         if isempty(result_MeanCorr_vs_r{count_r,count_n*2})==0
%             list_MeanCorr_vs_r(1,size(list_MeanCorr_vs_r,2)+1)=result_MeanCorr_vs_r{count_r,count_n*2};
%         end
%     end
%     result_AllCellsAvg_MeanCorr_vs_r(count_r,1)=(count_r-1)*pix_size;  %distance in um
%     result_AllCellsAvg_MeanCorr_vs_r(count_r,2)=mean(list_MeanCorr_vs_r);  %average correlation
%     result_AllCellsAvg_MeanCorr_vs_r(count_r,3)=std(list_MeanCorr_vs_r);
% end






%%





%%
'Check data'

grp1=1:80;
grp2=81:157;
% grp3=54:82;


figure
subgroups=[];
subgroups(grp1,1)=1;
subgroups(grp2,1)=2;
% subgroups(grp3,1)=3;
draw=0;
for count=41:60%size(All,2)
    draw=draw+1;
subplot(4,5,draw)
boxplot(All(:,count),subgroups(:,1),'Notch','marker','Widths',0.8)
end


%%


'PCA'

mean_all=mean(All,'omitNaN');
SD_all=std(All,'omitNaN');
standardized_all=(All-mean_all)./SD_all;


[coeff,score,latent,~,explained] = pca(standardized_all);


% Plot the total variance explained by each individual principal component along with the cumulative total that has been explained.
figure
plot([cumsum(latent(1:size(latent,1)))/sum(latent) latent(1:size(latent,1))/sum(latent)]*100,'.',...
    'MarkerSize',18)
xlabel('# of principal component');
ylabel('% of variance of dataset explained');
legend('Cumulative','Individual')
title('Individual and Cumulative Variance Explained by PCA')
grid on



% visualize first two principal components
figure
scatter(score(grp1,1),score(grp1,2),'r','filled','MarkerFaceAlpha',0.6)
hold on
scatter(score(grp2,1),score(grp2,2),'b','filled','MarkerFaceAlpha',0.6)
% scatter(score(grp3,1),score(grp3,2),'g','filled','MarkerFaceAlpha',0.6)
title('1st and 2nd PCs')
axis equal
xlabel('First Principal Component')
ylabel('Second Principal Component')


figure
scatter(score(grp1,3),score(grp1,4),'r','filled','MarkerFaceAlpha',0.6)
hold on
scatter(score(grp2,3),score(grp2,4),'b','filled','MarkerFaceAlpha',0.6)
% scatter(score(grp3,3),score(grp3,4),'g','filled','MarkerFaceAlpha',0.6)
title('3rd and 4th PCs')
axis equal
xlabel('Third Principal Component')
ylabel('Fourth Principal Component')


figure
scatter3(score(grp1,1),score(grp1,2),score(grp1,3),'r','filled','MarkerFaceAlpha',0.6)
hold on
scatter3(score(grp2,1),score(grp2,2),score(grp2,3),'b','filled','MarkerFaceAlpha',0.6)
% scatter3(score(grp3,1),score(grp3,2),score(grp3,3),'g','filled','MarkerFaceAlpha',0.6)
title('1st, 2nd and 3rd PCs')
axis equal
xlabel('First Principal Component')
ylabel('Second Principal Component')
zlabel('Third Principal Component')




% The first column of coeff contains the contributions from each signal to
% the first principal component. The values go from -1 to +1, with a larger
% absolute value denoting greater significance.




%%

