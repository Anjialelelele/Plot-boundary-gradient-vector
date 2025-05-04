function plotGradientVector(r, u_range, v_range,Delta,sur_color,n_iso_u, n_iso_v,n_color, n_scale,n_vectorStyle)
% Plot parametric surface and display gradient vectors on it
% Input parameters:
% r             - Function handle, accepts scalar parameters u and v, returns [x,y,z] coordinates
% u_range       - 1x2 array, parameter range [u_start, u_end], first free variable
% v_range       - 1x2 array, parameter range [v_start, v_end], second free variable
% Delta         - Sign of the Jacobian determinant,it is a one-dimensional array for two parameters and a two-dimensional array for three
% sur_color     - Optional, surface color (default: [0.5 0.5 1])
% n_iso_u       - Optional, number of u parameters for gradient vectors, together with n_iso_v determines the number of gradient vectors (default: 5)
% n_iso_v       - Optional, number of v parameters for gradient vectors (default: 5)
% n_color       - Optional, gradient vector color (default: red)
% n_scale       - Optional, scaling factor for vector length (default: 0.2)
% n_vectorStyle - Optional, gradient vector line style (default: '-')

    %Handle optional parameters
    if nargin<5
        sur_color=[0.5 0.5 1];
    end
    if nargin<6
        n_iso_u=5;
    end
    if nargin<7
        n_iso_v=5;
    end
    if nargin<8
        n_color="r";
    end
    if nargin<9
        n_scale=0.2;
    end
    if nargin<10
        n_vectorStyle='-';
    end

    % Set step size for numerical differentiation
    delta = 1e-6;
    
    % Generate surface mesh
    u = linspace(u_range(1), u_range(2), 100);
    v = linspace(v_range(1), v_range(2), 100);
    [U, V] = meshgrid(u, v);
    [X, Y, Z] = r(U, V);
    
    % Plot surface
    surf(X, Y, Z,'FaceColor',sur_color, 'EdgeColor', 'none');
    colormap(gray);
    hold on;
    camlight;
    lighting gouraud;
    
    % Calculate and plot cross product direction arrows
    u_values = linspace(u_range(1), u_range(2), n_iso_u);
    v_values = linspace(v_range(1), v_range(2), n_iso_v);
    for i = 1:length(u_values)
        for j = 1:length(v_values)
            u0 = u_values(i);
            v0 = v_values(j);
            
            up = min(u0 + delta, u_range(2));
            um = max(u0 - delta, u_range(1));
            [xp, yp, zp] = r(up, v0);
            [xm, ym, zm] = r(um, v0);
            du = [xp - xm, yp - ym, zp - zm] / (up - um);
            
            vp = min(v0 + delta, v_range(2));
            vm = max(v0 - delta, v_range(1));
            [xp, yp, zp] = r(u0, vp);
            [xm, ym, zm] = r(u0, vm);
            dv = [xp - xm, yp - ym, zp - zm] / (vp - vm);
            
            % Calculate cross product and plot arrows
            cross_vec = Delta(i,j)*cross(du, dv);
            cross_unit=cross_vec/norm(cross_vec);
            [x, y, z] = r(u0, v0);

            quiver3(x, y, z, cross_unit(1), cross_unit(2), cross_unit(3), n_scale, 'Color', n_color, 'LineWidth', 2,'MaxHeadSize', 1,'LineStyle',n_vectorStyle);
        end
    end
    hold on;
end