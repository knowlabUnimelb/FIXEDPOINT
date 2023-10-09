% angle = linspace(0,2*pi,360);
% x = cos(angle);
% y = sin(angle);
% hold off
% plot(x,y,'.r', 'MarkerSize', 69)
% axis('equal')

% saturation = [16 14 10]';        % levels of saturation (for output file)
% 
% colors = nan(numel(saturation), 3); % Preallocate matrix
% for i = 1:numel(saturation)
%     colors(i,:) =  munsell2rgb('10PB', 5, saturation(i)); % RGB values for hue 5R, brightness 5, saturation defined above
% end


fig = figure('WindowStyle', 'docked');

hold on

pos = [4.1 4 2 2];
rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'r')
axis equal

pos = [6.2 4 2 2];
rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'b')
axis equal

pos = [2 4 2 2];
rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'b')
axis equal