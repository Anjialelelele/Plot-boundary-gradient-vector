function Delta = computeDelta(Fr, Fl, omega, lambda, q1, q2)
% Check if it's the two-parameter case
if nargin < 4 || isempty(lambda)
    isTwoParam = true;
else
    isTwoParam = false;
end

h = 1e-8;          % Step size for numerical differentiation
tol = 1e-12;       % Tolerance for zero detection

if isTwoParam
    % Process the two-parameter case
    n = numel(omega);
    omega = omega(:); q1 = q1(:); q2 = q2(:);
    Delta = zeros(n, 1);
    
    for i = 1:n
        omg = omega(i);
        q1i = q1(i);
        q2i = q2(i);
        
        % Compute initial partial derivatives
        dFr_dq1 = (Fr(omg, q1i+h, q2i) - Fr(omg, q1i, q2i)) / h;
        dFr_dq2 = (Fr(omg, q1i, q2i+h) - Fr(omg, q1i, q2i)) / h;
        dFl_dq1 = (Fl(omg, q1i+h, q2i) - Fl(omg, q1i, q2i)) / h;
        dFl_dq2 = (Fl(omg, q1i, q2i+h) - Fl(omg, q1i, q2i)) / h;
        
        J1 = dFr_dq1 * dFl_dq2 - dFr_dq2 * dFl_dq1;
        
        if abs(J1) < tol
            % Compute alternative partial derivatives
            dFr_domega = (Fr(omg+h, q1i, q2i) - Fr(omg, q1i, q2i)) / h;
            dFl_domega = (Fl(omg+h, q1i, q2i) - Fl(omg, q1i, q2i)) / h;
            J2 = dFr_dq2 * dFl_domega - dFr_domega * dFl_dq2;
            Delta(i) = sign(J2);
        else
            Delta(i) = sign(J1);
        end
    end
    Delta = reshape(Delta, size(omega));
else
    % Process the three-parameter case
    [rows, cols] = size(omega);
    assert(isequal(size(lambda), [rows, cols]) && ...
           isequal(size(q1), [rows, cols]) && ...
           isequal(size(q2), [rows, cols]), ...
           'Input arrays must have consistent dimensions.');
    
    Delta = zeros(rows, cols);
    
    for i = 1:rows
        for j = 1:cols
            omg = omega(i,j);
            lam = lambda(i,j);
            q1ij = q1(i,j);
            q2ij = q2(i,j);
            
            dFr_dq1 = (Fr(omg, lam, q1ij+h, q2ij) - Fr(omg, lam, q1ij, q2ij)) / h;
            dFr_dq2 = (Fr(omg, lam, q1ij, q2ij+h) - Fr(omg, lam, q1ij, q2ij)) / h;
            dFl_dq1 = (Fl(omg, lam, q1ij+h, q2ij) - Fl(omg, lam, q1ij, q2ij)) / h;
            dFl_dq2 = (Fl(omg, lam, q1ij, q2ij+h) - Fl(omg, lam, q1ij, q2ij)) / h;
            
            J1 = dFr_dq1 * dFl_dq2 - dFr_dq2 * dFl_dq1;
            
            if abs(J1) < tol
                dFr_domega = (Fr(omg+h, lam, q1ij, q2ij) - Fr(omg, lam, q1ij, q2ij)) / h;
                dFl_domega = (Fl(omg+h, lam, q1ij, q2ij) - Fl(omg, lam, q1ij, q2ij)) / h;
                J2 = dFr_dq2 * dFl_domega - dFr_domega * dFl_dq2;
                Delta(i,j) = sign(J2);
            else
                Delta(i,j) = sign(J1);
            end
        end
    end
end
end