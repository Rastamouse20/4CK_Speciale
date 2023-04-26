%% Function made for plotting
function [] = plot2DElements(fignum,nodalCoordinates,elementConnectivity,displacement,displScale,elementValue,density)

    figure(fignum); clf;
    plottitle = 'Undeformed structure';
    
    if (~isempty(displacement))
        if (isempty(displScale))
            displScale = 0*max(max(nodalCoordinates(:,2)),max(nodalCoordinates(:,3)))/max(abs(displacement));
        end
        nodalCoordinates(:,2) = nodalCoordinates(:,2) + displScale*displacement(1:2:end);
        nodalCoordinates(:,3) = nodalCoordinates(:,3) + displScale*displacement(2:2:end);
        plottitle = ['Deformed structure - scale = ' num2str(displScale)];
    end
        
    if (~isempty(elementValue)) 
        if (size(elementConnectivity,1) < 100)
            patch('Vertices',nodalCoordinates(:,2:3),...
                  'Faces',elementConnectivity(:,2:5),...
                  'FaceVertexCData',elementValue,...
                  'FaceColor','flat');
        else
            patch('Vertices',nodalCoordinates(:,2:3),...
                  'Faces',elementConnectivity(:,2:5),...
                  'FaceVertexCData',elementValue,...
                  'FaceColor','flat','LineStyle','none');
        end
        axis off; colorbar;
    else
        patch('Vertices',nodalCoordinates(:,2:3),...
              'Faces',elementConnectivity(:,2:5),...
              'FaceColor','white');        
    end    
    axis tight; axis equal; title(plottitle);
    if (density); caxis([0 1]); colormap(flipud(gray)); end
    drawnow;
end