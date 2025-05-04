function plot_curves(paramFunc, wRange,Delta,color,nArrows, arrowScale,arrowColor)
% Plot parametric curves and display tangent and normal vector arrows
% Input parameters:
% paramFunc  - Function handle, accepts scalar parameter w, returns [x, y] coordinates
% wRange     - 1x2 array, parameter range [w_start, w_end]
% Delta      - Sign of the Jacobian determinant
% color      - Optional, curve color (default: blue)
% nArrows    - Optional, number of gradient vectors (default: 20)
% arrowScale - Optional, scaling factor for vector length (default: 0.1)
% arrowColor - Optional, color of gradient vectors (default: red)

% Handle optional parameters
if nargin<4
    color='b';
end
if nargin < 5
    nArrows = 20;
end
if nargin < 6
    arrowScale = 0.1;
end
if nargin<7
    arrowColor='r';
end

% Generate curve data
wCurve = linspace(wRange(1), wRange(2), 1000);
xCurve = zeros(size(wCurve));
yCurve = zeros(size(wCurve));

for i = 1:length(wCurve)
    xy = paramFunc(wCurve(i));
    xCurve(i) = xy(1);
    yCurve(i) = xy(2);
end

% Generate curve data
plot(xCurve, yCurve,'Color',color,'LineWidth',2);
hold on;

% Generate sample points and plot arrows
wSamples = linspace(wRange(1), wRange(2), nArrows);
h = 1e-5; % Step size for numerical differentiation

for i = 2:length(wSamples)
    w = wSamples(i);
    
    % Get current point coordinates
    xy = paramFunc(w);
    x = xy(1);
    y = xy(2);
    
    % Numerical derivative calculation
    xy_plus = paramFunc(w + h);
    xy_minus = paramFunc(w - h);
    dx = (xy_plus(1) - xy_minus(1)) / (2*h);
    dy = (xy_plus(2) - xy_minus(2)) / (2*h);
    
    % Handle zero vector case
    if norm([dx, dy]) < 1e-10
        continue;
    end
    
    % Calculate normal vector
    
    normal = Delta(i)*[dy, -dx];
    normal_unit = normal / norm(normal);
    
    % Plot normal vector
    quiver(x, y, normal_unit(1)*arrowScale, normal_unit(2)*arrowScale,...
        0, 'Color', arrowColor, 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
end
end