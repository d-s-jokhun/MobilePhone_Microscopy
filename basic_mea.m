

%%% written by D.S.JOKHUN on 04/09/2018



function [Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_mea(segmented_raw_2d,MetaData)

 MetaData.Num_of_Pixels_X=size(segmented_raw_2d,2);
 MetaData.Num_of_Pixels_Y=size(segmented_raw_2d,1);



Filename=cell(size(segmented_raw_2d,3),1);
Nuc_label=zeros(size(segmented_raw_2d,3),1);
Centroid=zeros(size(segmented_raw_2d,3),2);
Pro_area=zeros(size(segmented_raw_2d,3),1);
Perimeter=zeros(size(segmented_raw_2d,3),1);
AR=zeros(size(segmented_raw_2d,3),1);
Shape_factor=zeros(size(segmented_raw_2d,3),1);
PDI=zeros(size(segmented_raw_2d,3),1);;
Centre_mismatch=zeros(size(segmented_raw_2d,3),1);
I80_by_I20=zeros(size(segmented_raw_2d,3),1);
nHigh_by_nLow=zeros(size(segmented_raw_2d,3),1);
Mean_Norm_Int=zeros(size(segmented_raw_2d,3),1);
Median_Norm_Int=zeros(size(segmented_raw_2d,3),1);
SD_Norm_Int=zeros(size(segmented_raw_2d,3),1);
p_of_IntDistri_being_Normal=zeros(size(segmented_raw_2d,3),1);
Entropy=zeros(size(segmented_raw_2d,3),1);

Relative_concavity=zeros(size(segmented_raw_2d,3),1);
AvgCurva_X_Peri=zeros(size(segmented_raw_2d,3),1);
SDinCurva=zeros(size(segmented_raw_2d,3),1);
AvgPosCurva_X_Peri=zeros(size(segmented_raw_2d,3),1);
AvgNegCurva_X_Peri=zeros(size(segmented_raw_2d,3),1);
n_CurvaChangesSign=zeros(size(segmented_raw_2d,3),1);
MaxPosCurva=zeros(size(segmented_raw_2d,3),1);
nPeaks_Within80PercOf_MaxPos=zeros(size(segmented_raw_2d,3),1);
AvgProminance_of_PeaksWithin80PercOf_MaxPos=zeros(size(segmented_raw_2d,3),1);
AvgWidth_of_PeaksWithin80PercOf_MaxPos=zeros(size(segmented_raw_2d,3),1);
MaxNegCurva=zeros(size(segmented_raw_2d,3),1);
nPeaks_Within80PercOf_MaxNeg=zeros(size(segmented_raw_2d,3),1);
AvgProminance_of_PeaksWithin80PercOf_MaxNeg=zeros(size(segmented_raw_2d,3),1);
AvgWidth_of_PeaksWithin80PercOf_MaxNeg=zeros(size(segmented_raw_2d,3),1);


Boundary_RadOfCurva_K=cell(size(segmented_raw_2d,3),1);


for nuc_count=1:size(segmented_raw_2d,3);
    
    Filename{nuc_count,1}=MetaData.Filename;
    
    raw=segmented_raw_2d(:,:,nuc_count);
    bw=raw>0;
    stats=regionprops(bw,raw,'PixelIdxList','PixelValues','Area','MajorAxisLength','MinorAxisLength','Perimeter', 'Centroid','WeightedCentroid', 'ConvexArea');
    
    if size(stats,1)>1
        area_check=[];
        template=zeros(size(raw));
        for area_count=1:size(stats,1)
            area_check(area_count)=stats(area_count).Area;
        end
        [~,id]=max(area_check);
        template(stats(id).PixelIdxList)=1;
        raw=raw.*uint16(template);
    end
        
    bw=raw>0;
    stats=regionprops(bw,raw,'PixelIdxList','PixelValues','Area','MajorAxisLength','MinorAxisLength','Perimeter', 'Centroid','WeightedCentroid', 'ConvexArea');
    
    
    Norm_Pix_Intensities{1,nuc_count}=MetaData.Filename;
    Norm_Pix_Intensities(3:2+size(stats.PixelValues,1),nuc_count)=num2cell(mat2gray(stats.PixelValues));
    
    Mean_Norm_Int(nuc_count,1)=mean(mat2gray(stats.PixelValues));
    Median_Norm_Int(nuc_count,1)=median(mat2gray(stats.PixelValues));
    SD_Norm_Int(nuc_count,1)=std(mat2gray(stats.PixelValues));
    [~,p_of_IntDistri_being_Normal(nuc_count,1)] = kstest(zscore(mat2gray(stats.PixelValues)));
    
    Centroid(nuc_count,1)=stats.Centroid(1)*MetaData.Voxel_Size_X;
    Centroid(nuc_count,2)=stats.Centroid(2)*MetaData.Voxel_Size_Y;
    Pro_area(nuc_count,1)=(stats.Area * (MetaData.Voxel_Size_X*MetaData.Voxel_Size_Y));
    Perimeter(nuc_count,1)=(stats.Perimeter * MetaData.Voxel_Size_X);
    AR(nuc_count,1)=(stats.MajorAxisLength/stats.MinorAxisLength);
    Shape_factor(nuc_count,1)=((stats.Perimeter^2)/(4*pi*stats.Area));
    
    
    %% PDI
    distance_frm_cen_squared=zeros(MetaData.Num_of_Pixels_Y, MetaData.Num_of_Pixels_X);
    for countX=1:MetaData.Num_of_Pixels_X
        for countY=1:MetaData.Num_of_Pixels_Y
            distance_frm_cen_squared(countY,countX)=((sqrt(((countX-stats.Centroid(1))^2)...
                +(((MetaData.Num_of_Pixels_Y-countY+1)-stats.Centroid(2))^2)))...
                *MetaData.Voxel_Size_X).^2;
        end
    end
    actual_IntMoment_2=distance_frm_cen_squared.*mat2gray(raw);
    total_NormInt=sum(sum(mat2gray(raw)));
    uniform_IntMoment_2=distance_frm_cen_squared.*(double(bw)*(total_NormInt/stats.Area));
    PDI(nuc_count,1)=sum(sum(actual_IntMoment_2))/sum(sum(uniform_IntMoment_2));  %peripheral distributon index
    
    %% Centre mismatch
    Centre_mismatch(nuc_count,1)=(sqrt(((stats.WeightedCentroid(1)-stats.Centroid(1))^2)...
        +((stats.WeightedCentroid(2)-stats.Centroid(2))^2)))...
        *MetaData.Voxel_Size_X;
    
    
    %% Intensity Histogram Analysis
    I_exclude_percentiles=prctile(single(nonzeros(raw)),[0.1,99.9]);   %elimimating extreme values from the image (e.g saturated pixels etc.)
    aft_excl_extremes=raw.*uint16(raw>=I_exclude_percentiles(1)).*uint16(raw<=I_exclude_percentiles(2));
    
    I_percentiles=prctile(single(nonzeros(aft_excl_extremes)),[20,80]);
    I80_by_I20(nuc_count,1)=I_percentiles(2)/I_percentiles(1);
    
    normalize_aft_excl_extremes=mat2gray(nonzeros(aft_excl_extremes));
    nHigh_by_nLow(nuc_count,1)=sum(normalize_aft_excl_extremes>=0.8)/sum(normalize_aft_excl_extremes<=0.2);
    
    
    Entropy(nuc_count,1)=entropy(normalize_aft_excl_extremes);
        
    
    
    %% global curvature analysis
    Relative_concavity(nuc_count,1)=(stats.ConvexArea-stats.Area)/stats.ConvexArea;
    
    %% local curvature analysis
% % %     boundary=bwboundaries(bw);
% % %     boundary_xy=[];
% % %     boundary_xy(:,1)=boundary{1,1}(1:end-1,2);   %The last point is a repeat of the 1st point.
% % %     boundary_xy(:,2)=boundary{1,1}(1:end-1,1);
% % %     smooth_span=2/MetaData.Voxel_Size_X; %smooth span of 2um
% % %     cond = mod(smooth_span,2)<1;  % =1 if remainder is less than 1 (rounding down will bring it to an even num) and =0 if remainder is more than 1 (rounding down will bring it to an odd number).
% % %     smooth_span = floor(smooth_span);
% % %     smooth_span = smooth_span+cond;  % if cond was 1, 1 will be added to floor(span). if cond was 0, floor(span) is already odd and nothing is added.
% % %     %padding the x-y series for proper smoothing and proper subsequent fitting
% % %     boundary_xy((2*smooth_span)+1:(2*smooth_span)+size(boundary_xy,1),:)=boundary_xy(:,:); %sifting the series by some rows
% % %     boundary_xy(1:(2*smooth_span),:)=boundary_xy(end-(2*smooth_span)+1:end,:); %padding the top
% % %     boundary_xy(end+1:end+(2*smooth_span),:)=boundary_xy((2*smooth_span)+1:(2*smooth_span)+(2*smooth_span),:); %padding the bottom
% % %     boundary_xy(:,1)=smooth(boundary_xy(:,1),smooth_span,'lowess');
% % %     boundary_xy(:,2)=smooth(boundary_xy(:,2),smooth_span,'lowess');
% % %     
% % %     [MetaData.Filename,' nuc No.',num2str(nuc_count),' of ',num2str(size(segmented_raw_2d,3)),' -boundary analysis']
% % %     fit_span=1/MetaData.Voxel_Size_X; %%tangent 1 will be taken over 1um stretch before point P and tangent 2 will be take 1um stretch after point P.
% % %     cond = mod(fit_span,2)>=1;  % =1 if remainder is more or equal to 1 (rounding down will bring it to an odd num) and =0 if remainder is less than 1 (rounding down will bring it to an even number).
% % %     fit_span = floor(fit_span);
% % %     fit_span = fit_span+cond;  % if cond was 1, 1 will be added to floor(span). if cond was 0, floor(span) is already even and nothing is added.
% % %     
% % %     for point_count=(2*smooth_span)+1:size(boundary_xy,1)-(2*smooth_span);
% % % %         point_count
% % %         line_1=fit(boundary_xy(point_count-fit_span:point_count,1),boundary_xy(point_count-fit_span:point_count,2),'poly1');
% % %         coeffvalues_1=coeffvalues(line_1);
% % %         syms x y   %declaring system of equations
% % %         if abs(coeffvalues_1(1,1))>10^10
% % %             norm_eqn_1 = y==boundary_xy(point_count-(fit_span/2),2);
% % %             norm_grad_1=0;
% % %         end
% % %         if abs(coeffvalues_1(1,1))<10^-10
% % %             norm_eqn_1 = x==boundary_xy(point_count-(fit_span/2),1);
% % %             norm_grad_1=Inf;
% % %         end
% % %         if abs(coeffvalues_1(1,1)) <=10^10
% % %             if abs(coeffvalues_1(1,1)) >=10^-10
% % %                 norm_grad_1=-1/(coeffvalues_1(1,1));
% % %                 c1=(boundary_xy(point_count-(fit_span/2),2)) - (norm_grad_1*(boundary_xy(point_count-(fit_span/2),1))); %c=y-mx
% % %                 norm_eqn_1 = y==(norm_grad_1*x)+c1;
% % %             end
% % %         end
% % %         
% % %         line_2=fit(boundary_xy(point_count:point_count+fit_span,1),boundary_xy(point_count:point_count+fit_span,2),'poly1');
% % %         coeffvalues_2=coeffvalues(line_2);
% % %         if abs(coeffvalues_2(1,1))>10^10
% % %             norm_eqn_2 = y==boundary_xy(point_count+(fit_span/2),2);
% % %             norm_grad_2=0;
% % %         end
% % %         if abs(coeffvalues_2(1,1))<10^-10
% % %             norm_eqn_2 = x==boundary_xy(point_count+(fit_span/2),1);
% % %             norm_grad_2=Inf;
% % %         end
% % %         if abs(coeffvalues_2(1,1)) <=10^10
% % %             if abs(coeffvalues_2(1,1)) >=10^-10
% % %                 norm_grad_2=-1/(coeffvalues_2(1,1));
% % %                 c2=(boundary_xy(point_count+(fit_span/2),2)) - (norm_grad_2*(boundary_xy(point_count+(fit_span/2),1))); %c=y-mx
% % %                 norm_eqn_2 = y==(norm_grad_2*x)+c2;
% % %             end
% % %         end
% % %         
% % %         if round(norm_grad_1,10) == round(norm_grad_2,10)
% % %             radius_of_curvature=Inf;
% % %         else
% % %             sol = solve([norm_eqn_1, norm_eqn_2], [x, y]);
% % %             centre_of_curvature=double([sol.x,sol.y]);
% % %             radius_of_curvature=double((sqrt(((boundary_xy(point_count,1)-sol.x)^2)+((boundary_xy(point_count,2)-sol.y)^2))))*MetaData.Voxel_Size_X;
% % %                         
% % %             cen_of_curva_to_pt=[boundary_xy(point_count,1)-centre_of_curvature(1,1),boundary_xy(point_count,2)-centre_of_curvature(1,2)];   %vector to move from centre of curvature to point P on the boundary
% % %             unit_vec=cen_of_curva_to_pt./norm(cen_of_curva_to_pt);  % unit vector to move from centre of curvature to point P on the boundary
% % %             query_point=[boundary_xy(point_count,1)+(2*unit_vec(1,1)),boundary_xy(point_count,2)+(2*unit_vec(1,2))];
% % %             if inpolygon(query_point(1,1),query_point(1,2),boundary_xy(smooth_span+1:length(boundary_xy)-smooth_span,1),boundary_xy(smooth_span+1:length(boundary_xy)-smooth_span,2))==1    %If point P+2*(unit vector from the centre of curvature to point P) lies inside the cell, this statement will be true
% % %                 radius_of_curvature=-radius_of_curvature;
% % %             end
% % %         end
% % %         boundary_xy(point_count,3)=radius_of_curvature;
% % %         boundary_xy(point_count,4)=1/radius_of_curvature;
% % %     end
    
    Boundary_RadOfCurva_K_header={'x(pix)' 'y(pix)' 'Radius of curvature(um)' 'curvature,K=1/r(um^-1)'};
    pre_Boundary_RadOfCurva_K=[0 0 0 0];%boundary_xy((2*smooth_span)+1:end-(2*smooth_span),1:4);
    Boundary_RadOfCurva_K{nuc_count,1}=vertcat(Boundary_RadOfCurva_K_header,num2cell(pre_Boundary_RadOfCurva_K));
    
    AvgCurva_X_Peri(nuc_count,1)=0;%mean(pre_Boundary_RadOfCurva_K(:,4))*Perimeter(nuc_count,1);
    SDinCurva(nuc_count,1)=0;%std(pre_Boundary_RadOfCurva_K(:,4));
    AvgPosCurva_X_Peri(nuc_count,1)=0;% (mean(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)>0)))) *Perimeter(nuc_count,1);
    AvgNegCurva_X_Peri(nuc_count,1)=0;% (mean(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)<0)))) *Perimeter(nuc_count,1);
    
% % %     SignTest=diff((nonzeros(pre_Boundary_RadOfCurva_K(:,4)))>0);
    n_CurvaChangesSign(nuc_count,1)=0;%numel(nonzeros(SignTest));
    
    MaxPosCurva(nuc_count,1)= 0;%max(findpeaks(pre_Boundary_RadOfCurva_K(:,4)));
% % %     [pks,~,w,p] =findpeaks(pre_Boundary_RadOfCurva_K(:,4),'MinPeakDistance',round(3/MetaData.Voxel_Size_X),'MinPeakHeight',(0.8*MaxPosCurva(nuc_count,1)));   %If 2 peaks are reparated by less than 3um along the boundary, they are counted as 1
    nPeaks_Within80PercOf_MaxPos(nuc_count,1)=0;% numel(pks); 
    AvgProminance_of_PeaksWithin80PercOf_MaxPos(nuc_count,1)=0;%mean(p);
    AvgWidth_of_PeaksWithin80PercOf_MaxPos(nuc_count,1)=0;%mean(w);
    MaxNegCurva(nuc_count,1)=0;% -1*(max(findpeaks(-1*(pre_Boundary_RadOfCurva_K(:,4)))));
% % %     [pks,~,w,p] =findpeaks(-1*(pre_Boundary_RadOfCurva_K(:,4)),'MinPeakDistance',round(3/MetaData.Voxel_Size_X),'MinPeakHeight',abs(0.8*MaxNegCurva(nuc_count,1)));
    nPeaks_Within80PercOf_MaxNeg(nuc_count,1)=0;% numel(pks);
    AvgProminance_of_PeaksWithin80PercOf_MaxNeg(nuc_count,1)=0;%mean(p);
    AvgWidth_of_PeaksWithin80PercOf_MaxNeg(nuc_count,1)=0;%mean(w);
    
    
    

%     figure ('Name',['Curvature,K of ',MetaData.Filename])
%     colormap jet
%     patch(pre_Boundary_RadOfCurva_K(:,1),pre_Boundary_RadOfCurva_K(:,2),pre_Boundary_RadOfCurva_K(:,4),'EdgeColor','interp','FaceColor','none','LineWidth',5)
%     
%     figure ('Name',['Curvature,K of ',MetaData.Filename])
%     plot(pre_Boundary_RadOfCurva_K(:,4))
    
    
    
    
    
end


Basic_Measurements=horzcat(Filename,num2cell(Nuc_label),num2cell(Centroid),num2cell(Pro_area),num2cell(Perimeter),num2cell(AR),num2cell(Shape_factor),num2cell(PDI),num2cell(Centre_mismatch),num2cell(I80_by_I20),num2cell(nHigh_by_nLow),num2cell(Mean_Norm_Int),num2cell(Median_Norm_Int),num2cell(SD_Norm_Int),num2cell(p_of_IntDistri_being_Normal),num2cell(Entropy),...
    num2cell(Relative_concavity),num2cell(AvgCurva_X_Peri),num2cell(SDinCurva),num2cell(AvgPosCurva_X_Peri),num2cell(AvgNegCurva_X_Peri),num2cell(n_CurvaChangesSign)...
    ,num2cell(MaxPosCurva),num2cell(nPeaks_Within80PercOf_MaxPos),num2cell(AvgProminance_of_PeaksWithin80PercOf_MaxPos),num2cell(AvgWidth_of_PeaksWithin80PercOf_MaxPos)...
    ,num2cell(MaxNegCurva),num2cell(nPeaks_Within80PercOf_MaxNeg),num2cell(AvgProminance_of_PeaksWithin80PercOf_MaxNeg),num2cell(AvgWidth_of_PeaksWithin80PercOf_MaxNeg));


end









