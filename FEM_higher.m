function [S, A, Al] = FEM_higher(M, p, bc, ICache)
%[S, A] = FEM_HIGHER(M, p) Compute the stiffness (S) and mass (A) matrices
%for a order p Finite Element discretization of Laplace-Beltrami operator
%on a mesh M. Uses Neumann boundary conditions.
%
%[S, A, Al] = FEM_HIGHER(M, p) Also compute the lumped mass matrix. Throws
%a warning if the lumped mass results singular.
%
%[--] = FEM_HIGHER(M, p, bc) Also allows to choose the boundary conditions.
%Permitted choices are:
%   - 'Neumann': Uses Neumann boundary conditions.
%   - 'Dirichlet': Uses Dirichlet boundary conditions.
%
%[--] = FEM_HIGHER(M, p, bc, ICache) Since the computation of the I#
%matrices is costly, the function caches them on disk after the first time 
%they are computed. The optional parameter ICache tells to the function the
%directory where such matrices can be found. Notice the matrices are also
%saved in that directory in a single .mat file named "OrderP.mat".
%
%
%
%Dependencies:          None, the function only depends on onther functions
%                       in the same file.
%Required MAT-Files:    None, the function caches the files it needs at the
%                       first execution with a certain order.
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli@di.uniroma1.it
%Last Revision: 15 June 2020


    % Default arguments
    if nargin < 3
        bc = "Neumann";
    end
    if nargin < 4
        ICache = "Cache/FEM/";
    end
    
    % Arguments validity
    if lower(bc) ~= "neumann" && lower(bc) ~= "dirichlet"
        errmsg = strcat("Accepted boundary conditions are 'Neumann' ", ...
                        sprintf("and 'Dirichlet'. Input is %s.", bc));
        error(errmsg);
    end

    % Compute the I matrices
    [I1, I2, I3, I4] = calc_FEM_test_funcs(p, ICache);
    q = size(I1, 1);
    
    % Compute the boundary edges
    E = [M.TRIV(:, [1 2]);
         M.TRIV(:, [2 3]);
         M.TRIV(:, [3 1])];
    E = sort(E, 2);
    [Eu, Euidx, ~] = unique(E, 'rows');
    E = E(setdiff(1:length(E), Euidx), :);
    Boundary = Eu(~ismember(Eu, E, 'rows'), :);
    
    % Linear FEM ==> use cotangent method
    if p == 1
        if nargout > 2
            [S, A, Al] = cotangent_method(M, bc, Boundary);
        else
            [S, A] = cotangent_method(M, bc, Boundary);
        end
        return;
    end
    
    % Compute the adjacency matrix and the matrices entries and types
    Ad = calc_indices_adj_order_p(M, p);
    Nedges = nnz(triu(Ad));
    OnEdges = p - 1;
    Internal = (p + 1) * (p + 2) / 2 - (3 * p);
    Ntot = M.n + OnEdges * Nedges + Internal * M.m;
    
    % Compute TRItot
    TRItot = zeros(M.m, q);
    TRItot(:, 1:3) = M.TRIV;
    % On edges, consistency is achieved assigning increasing values to
    % virtual nodes if the vertices are in ascending order, otherwise the
    % assigned values are decreasing
    ExtraBoundary = [];
    for ii = 1:OnEdges
        i = ii + 3;
        TRItot(:, i) = Ad(sub2ind(size(Ad), M.TRIV(:, 1), M.TRIV(:, 2))) + M.n + ...
                       (OnEdges - ii) .* (M.TRIV(:, 1) > M.TRIV(:, 2)) + ...
                       (ii - 1) .* (~(M.TRIV(:, 1) > M.TRIV(:, 2)));
        % Virtual nodes on boundary are identified as the virtual nodes on
        % boundary edges
        ExtraBoundary = [ExtraBoundary;
                         TRItot(ismember(sort(M.TRIV(:, [1 2]), 2), Boundary, 'rows'), i)];
        
        i = i + OnEdges;
        TRItot(:, i) = Ad(sub2ind(size(Ad), M.TRIV(:, 2), M.TRIV(:, 3))) + M.n + ...
                       (OnEdges - ii) .* (M.TRIV(:, 2) > M.TRIV(:, 3)) + ...
                       (ii - 1) .* (~(M.TRIV(:, 2) > M.TRIV(:, 3)));
        ExtraBoundary = [ExtraBoundary;
                         TRItot(ismember(sort(M.TRIV(:, [2 3]), 2), Boundary, 'rows'), i)];
        
        i = i + OnEdges;
        TRItot(:, i) = Ad(sub2ind(size(Ad), M.TRIV(:, 3), M.TRIV(:, 1))) + M.n + ...
                       (OnEdges - ii) .* (M.TRIV(:, 3) > M.TRIV(:, 1)) + ...
                       (ii - 1) .* (~(M.TRIV(:, 3) > M.TRIV(:, 1)));
        ExtraBoundary = [ExtraBoundary;
                         TRItot(ismember(sort(M.TRIV(:, [3 1]), 2), Boundary, 'rows'), i)];
    end
    % Consistency on internal nodes is easy and is achieved by simply
    % assigning arbitrary values (still ensuring the values are different)
    for i = 1:Internal
        TRItot(:, q - Internal + i) = (1:M.m)' + M.n + OnEdges * Nedges + (i - 1) * M.m;
    end
    

    % Compute the edges
    P1 = M.VERT(M.TRIV(:, 2), :) - M.VERT(M.TRIV(:, 1), :);
    P2 = M.VERT(M.TRIV(:, 3), :) - M.VERT(M.TRIV(:, 1), :);
    % Compute the edges' dot products
    P11 = dot(P1, P1, 2);
    P22 = dot(P2, P2, 2);
    P12 = dot(P1, P2, 2);
    % Compute the norms of the edges' cross products
    pre = cross(P1, P2, 2);
    pre = sqrt(sum(pre.^2, 2));
    % Compute the local mini-stiffness matrices
    Alocal = reshape(P11, 1, 1, M.m) .* repmat(I2, 1, 1, M.m) + ...
             reshape(P22, 1, 1, M.m) .* repmat(I1, 1, 1, M.m) - ...
             reshape(P12, 1, 1, M.m) .* repmat(I3, 1, 1, M.m);
    Alocal = Alocal ./ reshape(pre, 1, 1, M.m);
    % Compute the local mini-mass matrices
    Blocal = repmat(I4, 1, 1, M.m) .* reshape(pre, 1, 1, M.m);
    % Compute the indices to transfer local mini-mass/stiffness matrices to
    % global mass/stiffness
    AI = reshape(TRItot(:, repmat(1:q, 1, q))', q*q*M.m, 1);
    AJ = reshape(TRItot(:, reshape(repmat(1:q, q, 1), 1, q*q))', q*q*M.m, 1);
    BI = reshape(TRItot(:, repmat(1:q, 1, q))', q*q*M.m, 1);
    BJ = reshape(TRItot(:, reshape(repmat(1:q, q, 1), 1, q*q))', q*q*M.m, 1);
    
    % Compute the boundary vertices
    BoundaryVerts = [Boundary(:); ExtraBoundary];
    InternalVerts = setdiff(1:Ntot, BoundaryVerts);
    
    % S = Stiffness
    if lower(bc) == "neumann"
        S = sparse(AI, AJ, reshape(Alocal, 1, q*q*M.m), Ntot, Ntot);
    else
        S_full = sparse(AI, AJ, reshape(Alocal, 1, q*q*M.m), Ntot, Ntot);
        S = sparse(Ntot, Ntot);
        S(BoundaryVerts, BoundaryVerts) = eye(length(BoundaryVerts));
        S(InternalVerts, InternalVerts) = S_full(InternalVerts, InternalVerts);
    end
    
    % A = Mass
    if lower(bc) == "neumann"
        A = sparse(BI, BJ, reshape(Blocal, 1, q*q*M.m), Ntot, Ntot);
    else
        A_full = sparse(BI, BJ, reshape(Blocal, 1, q*q*M.m), Ntot, Ntot);
        A = sparse(Ntot, Ntot);
        A(InternalVerts, InternalVerts) = A_full(InternalVerts, InternalVerts);
    end
    
    if nargout > 2
        % Al = Lumped Mass
        Al = sparse(1:Ntot, 1:Ntot, sum(A), Ntot, Ntot, Ntot);
        % Dirichlet lumped mass is for sure singular. Don't even check
        if nnz(Al) < size(Al, 1) && lower(bc) == "neumann"
            warning("The lumped mass is singular (%d non-zeros out of %d)", ...
                    nnz(Al), size(Al, 1));
        end
    end

end






function Ad = calc_indices_adj_order_p(M, p)

indicesI = [M.TRIV(:,1); M.TRIV(:,2); M.TRIV(:,3); M.TRIV(:,3); M.TRIV(:,2); M.TRIV(:,1)];
indicesJ = [M.TRIV(:,2); M.TRIV(:,3); M.TRIV(:,1); M.TRIV(:,2); M.TRIV(:,1); M.TRIV(:,3)];
Ad = sparse(indicesI, indicesJ, ones(M.m,6), M.n, M.n);
[indicesI, indicesJ, ~] = find(triu(Ad));
Nedges = nnz(triu(Ad));
Ad = sparse(indicesI, indicesJ, 1:p-1:(p-1)*Nedges, M.n, M.n);
Ad = Ad + Ad';

end


function [W, Sc, Sl] = cotangent_method(M, bc, Boundary)
    % Stiffness (p.s.d.)
    angles = zeros(M.m,3);
    for i=1:3
        a = mod(i-1,3)+1;
        b = mod(i,3)+1;
        c = mod(i+1,3)+1;
        ab = M.VERT(M.TRIV(:,b),:) - M.VERT(M.TRIV(:,a),:);
        ac = M.VERT(M.TRIV(:,c),:) - M.VERT(M.TRIV(:,a),:);
        %normalize edges
        ab = ab ./ (sqrt(sum(ab.^2,2))*[1 1 1]);
        ac = ac ./ (sqrt(sum(ac.^2,2))*[1 1 1]);
        % normalize the vectors
        % compute cotan of angles
        angles(:,a) = cot(acos(sum(ab.*ac,2)));
        %cotan can also be computed by x/sqrt(1-x^2)
    end

    indicesI = [M.TRIV(:,1);M.TRIV(:,2);M.TRIV(:,3);M.TRIV(:,3);M.TRIV(:,2);M.TRIV(:,1)];
    indicesJ = [M.TRIV(:,2);M.TRIV(:,3);M.TRIV(:,1);M.TRIV(:,2);M.TRIV(:,1);M.TRIV(:,3)];
    values   = [angles(:,3);angles(:,1);angles(:,2);angles(:,1);angles(:,3);angles(:,2)]*0.5;
    W = sparse(indicesI, indicesJ, -values, M.n, M.n);
    W = W-sparse(1:M.n,1:M.n,sum(W));

    % Mass
    areas = mesh.proc.tri_areas(M);
    indicesI = [M.TRIV(:,1);M.TRIV(:,2);M.TRIV(:,3);M.TRIV(:,3);M.TRIV(:,2);M.TRIV(:,1)];
    indicesJ = [M.TRIV(:,2);M.TRIV(:,3);M.TRIV(:,1);M.TRIV(:,2);M.TRIV(:,1);M.TRIV(:,3)];
    values   = [areas(:); areas(:); areas(:); areas(:); areas(:); areas(:)]./12;
    Sc = sparse(indicesI, indicesJ, values, M.n, M.n);
    Sc = Sc+sparse(1:M.n, 1:M.n, sum(Sc));
    
    if lower(bc) == "dirichlet"
        Boundary = Boundary(:);
        Inner = setdiff(1:M.n, Boundary);
        W_full = W;
        W = sparse(M.n, M.n);
        W(Boundary, Boundary) = eye(length(Boundary));
        W(Inner, Inner) = W_full(Inner, Inner);
        Sc_full = Sc;
        Sc = sparse(M.n, M.n);
        Sc(Inner, Inner) = Sc_full(Inner, Inner);
    end
    
    % Lumped mass
    if nargout > 2
        Sl = spdiags(sum(Sc,2), 0, M.n, M.n);
    end
    
    W = (W + W')/2;
    Sc = (Sc + Sc')/2;
end

function [I1, I2, I3, I4] = calc_FEM_test_funcs(Order, ICache)

    [I1, I2, I3, I4, found] = check_cache(Order, ICache);
    if found
        return;
    end

    syms u v;
    U(1) = u^0;
    V(1) = v^0;
    for i = 1:Order
        U(i + 1) = u * U(i);
        V(i + 1) = v * V(i);
    end
    
    UV = V' * U;
    UV = fliplr(UV);
    Vars = [];
    for i = 0:Order
        Vars = [Vars; diag(UV, Order - i)];
    end
    Vars = Vars';
    
    M = [];
    M = [subs(subs(Vars, u, 0), v, 0);
         subs(subs(Vars, u, 1), v, 0);
         subs(subs(Vars, u, 0), v, 1)];
    Space = linspace(0, 1, Order + 1);
    for i = 2:Order
        M = [M; subs(subs(Vars, u, Space(i)), v, 0)];
    end
    for i = 2:Order
        M = [M; subs(subs(Vars, u, 1 - Space(i)), v, Space(i))];
    end
    for i = 2:Order
        M = [M; subs(subs(Vars, u, 0), v, 1 - Space(i))];
    end
    for i = 2:Order
        for j = 2:Order - i + 1
            M = [M; subs(subs(Vars, u, Space(i)), v, Space(j))];
        end
    end
    
    % Compute coefficients for each test function
    C = zeros(length(Vars));
    for i = 1:length(Vars)
        rhs = zeros(length(Vars), 1);
        rhs(i) = 1;
        C(:, i) = M \ rhs;
    end
    H = C' * Vars';
    
    
    I1 = zeros(length(Vars));
    I2 = zeros(length(Vars));
    I3 = zeros(length(Vars));
    I4 = zeros(length(Vars));
    for i = 1:length(Vars)
        for j = i:length(Vars)
            f = H(i) * H(j);
            I4(i, j) = int(int(f, v, 0, 1 - u), u, 0, 1);
            f = diff(H(i), u) * diff(H(j), u);
            I1(i, j) = int(int(f, v, 0, 1 - u), u, 0, 1);
            f = diff(H(i), v) * diff(H(j), v);
            I2(i, j) = int(int(f, v, 0, 1 - u), u, 0, 1);
            f = diff(H(i), v) * diff(H(j), u) + diff(H(i), u) * diff(H(j), v);
            I3(i, j) = int(int(f, v, 0, 1 - u), u, 0, 1);
        end 
    end
    
    
    I1 = I1 + I1';
    I1(1:length(Vars)+1:end) = diag(I1) ./ 2;
    
    I2 = I2 + I2';
    I2(1:length(Vars)+1:end) = diag(I2) ./ 2;
    
    I3 = I3 + I3';
    I3(1:length(Vars)+1:end) = diag(I3) ./ 2;
    
    I4 = I4 + I4';
    I4(1:length(Vars)+1:end) = diag(I4) ./ 2;
    
    cachename = sprintf("%s/Order%d.mat", ICache, Order);
    save(cachename, 'I1', 'I2', 'I3', 'I4');
    
end


function [I1, I2, I3, I4, found] = check_cache(p, cachedir)
    cachename = sprintf("%s/Order%d.mat", cachedir, p);
    if exist(cachename, 'file')
        Cache = load(cachename);
        I1 = Cache.I1;
        I2 = Cache.I2;
        I3 = Cache.I3;
        I4 = Cache.I4;
        found = true;
    else
        I1 = [];
        I2 = [];
        I3 = [];
        I4 = [];
        found = false;
    end
end


