% Visualizes rocket flight using InertialExport simulation extension for OpenRocket

clear; clc; close;
ORDat = ORDataImport(0.01);
quatList = [ORDat.qW ORDat.qX ORDat.qY ORDat.qZ];
posList = [ORDat.relPosX ORDat.relPosY ORDat.relPosZ];

%%
% for i = 1:2:height(ORDat)
%     q = quaternion(quatList(i,:));
%     pos = posList(i,:);
%     poseplot(q,pos,'ScaleFactor',0.2);
% 
%     timeStr = sprintf('Time: %.2f s', ORDat.t(i));
%     title(timeStr, 'FontSize', 14); 
% 
%     drawnow;
% end

%%
    
% v = VideoWriter('flight_of_fuckin_bumblebee.mp4','MPEG-4');
% v.FrameRate = 60;
% open(v);

    q = quaternion(quatList(1,:));
    patch = poseplot(q);
    patch.ScaleFactor = 15;

    xlabel("X")
    ylabel("Y")
    zlabel("Z")
    
    for i = 1:4:height(ORDat)
        q = quaternion(quatList(i,:));
        pos = posList(i,:);
        set(patch, Orientation=q, Position=pos); hold on
        plot3(pos(1), pos(2), pos(3), '.b', 'MarkerSize', 2); hold on
    
        xlim([-100, 100]);
        ylim([-100, 100]);
        zlim([0, 400]);
    
        set(gca,'ZDir','normal')  
        title(sprintf("t = %0.2f", ORDat.t(i)))
        drawnow
        % frame = getframe(gcf);
        % writeVideo(v, frame);
    end
% close(v);