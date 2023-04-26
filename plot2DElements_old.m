function [] = plot2DElements_old(fignum,nodalCoordinates,elementConnectivity,displacement,elementStress)

    % % % Self-explanatory input

    figure(fignum); clf;
    plottitle = 'Undeformed structure';
    
    if (~isempty(displacement))
        displScale = 0.1*max(max(nodalCoordinates(:,2)),max(nodalCoordinates(:,3)))/max(abs(displacement));
        nodalCoordinates(:,2) = nodalCoordinates(:,2) + displScale*displacement(1:2:end);
        nodalCoordinates(:,3) = nodalCoordinates(:,3) + displScale*displacement(2:2:end);
        plottitle = ['Deformed structure - scale = ' num2str(displScale)];
    end
        
    if (~isempty(elementStress)) 
        elementColour = (elementStress - min(elementStress))./(max(elementStress)-min(elementStress));
        if (size(elementConnectivity,1) < 100)
            patch('Vertices',nodalCoordinates(:,2:3),...
                  'Faces',elementConnectivity(:,2:5),...
                  'FaceVertexCData',elementStress,...
                  'FaceColor','flat');
        else
            patch('Vertices',nodalCoordinates(:,2:3),...
                  'Faces',elementConnectivity(:,2:5),...
                  'FaceVertexCData',elementStress,...
                  'FaceColor','flat','LineStyle','none');
        end
        axis off; colorbar;
    else
        patch('Vertices',nodalCoordinates(:,2:3),...
              'Faces',elementConnectivity(:,2:5),...
              'FaceColor','white');        
    end    
    axis tight; axis equal; 
    title(plottitle);
    drawnow;
end