
%%% written by D.S.JOKHUN on 20/06/2018


function [Selected_2d_Raw]=manual_obj_selection(Proper_2d_segments)

Selected_2d_Raw=[];


        
        figure ('Name','Select or not? ');
        subplot(1,2,1)
        imshow(Proper_2d_segments,[])
        subplot(1,2,2)
        selection=bwselect(Proper_2d_segments>0);
        close all
        
        
        if sum(sum(selection))>0
            Selected_2d_Raw=Proper_2d_segments;
        end


end
