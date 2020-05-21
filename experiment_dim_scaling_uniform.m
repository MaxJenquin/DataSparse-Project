%% Numerical Experiment: fixed dataset scale, uniform grid, dimensionality scaling
% This experiment probes the scaling of different SKI paradigms (Toeplitz
% structure exploiting, RS-factorization exploiting, HODLR-factorization
% exploiting, Toeplitz-HODLR-factorization exploiting and default SKI) for
% a kernel with axis-multiplicative structure as problem dimensionality 
% grows, for a fixed-datascale problem.
 
% This experiment is complementary to the experiment "fixed dimensionality,
% uniform grid, dataset scaling", which scales dataset size
% with a fixed problem dimensionality instead of
% scaling dimensionality with a fixed data-per-axis size. That experiment is
% located in "experiment_data_scaling_uniform.m"
 
% In order to include the Toeplitz methodologies we use a stationary
% (Matern 3/2) kernel as the per-axis kernel factor. It is also necessary
% that inputs are on a grid with uniform spacing (in each dimension,
% although different axes need not have the same spacing), so this
% experiment is performed on a connected domain [-1,1]^3.
 
%% Initialization
 
% problem dimensionalities
min_num_dims = 2; % 1D case means trivial Kronecker product, which changes some helper functions
max_num_dims = 5;
 
% dataset size parameters:
% each dimension will have n grid locations, 
% where n is in [data_low, data_high]
data_low = 10; data_high = 15;
 
% kernel parameters: Matern
maternparam = [1.2, 2.3];
% noise variance:
noise = 1e-2;
 
% conjugate gradients error tolerance
cg_tol = 1e-12;
 
% RS factorization tolerance
rs_tol = 1e-12;
 
% containers for output information:
% timing
% default/toeplitz:
% (:,1)-total time needed for solve without linsolve
% (:,2)-total time needed for linsolve via eigendecomp
% (:,3)-total time needed for linsolve via CG
default_times = zeros(max_num_dims, 3);  
toeplitz_times = zeros(max_num_dims, 3);
% the rest skip the eigendecomp, so no linear solve shortcuts (CG only)
% therefore: 
% (:,1)-total time needed for solve without linsolve
% (:,2)-total time neede for linsolve via CG
rskelf_times = zeros(max_num_dims, 2);
hodlr_times = zeros(max_num_dims, 2);
toeplitz_hodlr_times = zeros(max_num_dims, 2);
% problem sizes
problem_sizes = zeros(max_num_dims, 1);
% linear CG iterations, for reference
default_iters = zeros(max_num_dims, 1);
toeplitz_iters = zeros(max_num_dims, 1);
rskelf_iters = zeros(max_num_dims, 1);
hodlr_iters = zeros(max_num_dims, 1);
toeplitz_hodlr_iters = zeros(max_num_dims, 1);
 
%% Main Loop, over problem dimensionality
 
dimstep = 0;
for num_dims=min_num_dims:max_num_dims
    dimstep = dimstep + 1;
    
    %%%%%%%%%%%
    % SECTION % Domain Generation
    %%%%%%%%%%%
    
    domain = cell(1, num_dims);
    size_counter = zeros(1,num_dims);
    for axis=1:num_dims
        size_counter(axis) = randsample(data_low:data_high,1);
        domain{axis} = linspace(-1,1,size_counter(axis))';
    end
    problem_sizes(dimstep) = prod(size_counter);
    
    % Additional setup: generate a right-hand-side for mat-vecs and linear
    % solves. The values are arbitrary.
    rhs = -2.5 + 5*rand(problem_sizes(dimstep), 1);
    
    % generate component kernel matrices for Default, Toeplitz, HODLR, and
    % HODLR-Toeplitz methodologies
    K_cells = cell(1,num_dims);
    for axis=1:num_dims
        K_cells{axis} = Kmatern(domain{axis}, domain{axis}, maternparam);
    end
    
    fprintf('Beginning dimension-%i experiments with %i points\n', [num_dims, problem_sizes(dimstep)])
    
    %%%%%%%%%%%
    % SECTION % Default SKI, Toeplitz SKI
    %%%%%%%%%%%
    
    fprintf('Toeplitz SKI and default SKI\n')
    
    eig_decomp = cell(1, num_dims);
    vecs_decomp = cell(1, num_dims);
    vecs_decomp_trps = cell(1, num_dims);
    
    tstart = tic;
    % compute eigendecomposition and logdets - used for both
    full_eigs = 1;
    for axis=1:num_dims
        [vecs_decomp{axis}, foo] = eig(K_cells{axis});
        eig_decomp{axis} = diag(foo);
        vecs_decomp_trps{axis} = vecs_decomp{axis}';
        full_eigs = kron(full_eigs, eig_decomp{axis});
    end
    % compute log-determinant via eigendecomposition
    logdet_via_eigs = sum(reallog(full_eigs+noise^2));
    % time for logdet via eigendecomposition
    ttemp = toc(tstart);
    default_times(dimstep, 1) = ttemp;
    toeplitz_times(dimstep, 1) = ttemp;
    
    % DEFAULT SKI
    tstart = tic;
    % compute linear solve via eigendecomposition
    inv_diag = 1./(full_eigs+noise^2);
    % combining the next lines into one is finicky, hence separation
    trps_mult = kron_default_mult(vecs_decomp_trps, 0, rhs);
    linsolve_via_eigs = kron_default_mult(vecs_decomp, 0, inv_diag.*trps_mult);
    ttemp = toc(tstart);
    default_times(dimstep, 2) = ttemp;
    
    tstart = tic;
    % compute linear solve via conjugate gradients
    [linsolve_via_cg, default_iters(dimstep)] = my_lin_cg(@kron_default_mult,rhs,K_cells,noise^2,cg_tol);
    ttemp = toc(tstart);
    default_times(dimstep, 3) = ttemp;
    
    % clear reused variables for timing fairness
    clear inv_diag trps_mult linsolve_via_eigs linsolve_via_cg
    
    % TOEPLITZ SKI
    tstart = tic;
    % compute linear solve via eigendecomposition
    inv_diag = 1./(full_eigs+noise^2);
    % combining the next lines into one is finicky, hence separation
    trps_mult = kron_toeplitz_mult(vecs_decomp_trps, 0, rhs);
    linsolve_via_eigs = kron_toeplitz_mult(vecs_decomp, 0, inv_diag.*trps_mult);
    ttemp = toc(tstart);
    toeplitz_times(dimstep, 2) = ttemp;
    
    tstart = tic;
    % compute linear solve via conjugate gradients
    [linsolve_via_cg, toeplitz_iters(dimstep)] = my_lin_cg(@kron_toeplitz_mult,rhs,K_cells,noise^2,cg_tol);
    ttemp = toc(tstart);
    toeplitz_times(dimstep, 3) = ttemp;
    
    clear eig_decomp vecs_decomp vecs_decomp_trps full_eigs foo logdet_via_eigs
    clear inv_diag trps_mult linsolve_via_eigs linsolve_via_cg
    
    %%%%%%%%%%%
    % SECTION % RS SKI
    %%%%%%%%%%%
    
    fprintf('RS SKI\n')
    
    rs_factors = cell(1,num_dims);
    logdets = zeros(1,num_dims);
    
    tstart = tic;
    for axis=1:num_dims
        Afunwrapper = @(i,j)Afun_matern(i,j,domain{axis}',maternparam);
        rs_factors{axis} = rskelf(Afunwrapper, domain{axis}', 64, rs_tol);
        logdets(axis) = real(rskelf_logdet(rs_factors{axis}));
    end
    % compute logdet and correct using matrix determinant lemma
    logdet_via_rs = sum(logdets.*size_counter);  % logdet of base kernel matrix, no noise
    logdet_via_rs = logdet_via_rs + detlemma_correction_rskelf(rs_factors,noise^2);
    ttemp = toc(tstart);
    rskelf_times(dimstep, 1) = ttemp;
    % compute inverse via conjugate gradients
    tstart = tic;
    [linsolve_via_cg, rskelf_iters(dimstep)] = my_lin_cg(@kron_rskelf_mult,rhs,rs_factors,noise^2,cg_tol);
    ttemp = toc(tstart);
    rskelf_times(dimstep, 2) = ttemp;
    
    clear rs_factors logdets Afunwrapper logdet_via_rs linsolve_via_cg
    
    %%%%%%%%%%%
    % SECTION % HODLR SKI, Toeplitz-HODLR SKI
    %%%%%%%%%%%
    
    fprintf('HODLR SKI and Toeplitz-HODLR SKI\n')
    
    % HODLR SKI
    hr_factors = cell(1, num_dims);
    logdets = zeros(1, num_dims);
    
    tstart = tic;
    for axis=1:num_dims
        hr_factors{axis} = hodlr(K_cells{axis});
        logdets(axis) = reallog(det(hr_factors{axis}));
    end
    % compute logdet and correct using matrix determinant lemma
    logdet_via_hr = sum(logdets.*size_counter);
    logdet_via_hr = logdet_via_hr + detlemma_correction_hodlr(hr_factors,noise^2);
    ttemp = toc(tstart);
    hodlr_times(dimstep, 1) = ttemp;
    % compute inverse via conjugate gradients
    tstart = tic;
    [linsolve_via_cg, hodlr_iters(dimstep)] = my_lin_cg(@kron_hodlr_mult,rhs,hr_factors,noise^2,cg_tol);
    ttemp = toc(tstart);
    hodlr_times(dimstep, 2) = ttemp;
    
    clear hr_factors logdets logdet_via_hr linsolve_via_cg
    
    % Toeplitz-HODLR SKI
    t_hr_factors = cell(1,num_dims);
    t_logdets = zeros(1, num_dims);
    
    tstart = tic;
    for axis=1:num_dims
        t_hr_factors{axis} = hodlr('toeplitz',K_cells{axis}(:,1),K_cells{axis}(1,:));
        t_logdets(axis) = reallog(det(t_hr_factors{axis}));
    end
    % compute logdet and correct using matrix determinant lemma
    logdet_via_t_hr = sum(t_logdets.*size_counter);
    logdet_via_t_hr = logdet_via_t_hr + detlemma_correction_hodlr(t_hr_factors,noise^2);
    ttemp = toc(tstart);
    toeplitz_hodlr_times(dimstep, 1) = ttemp;
    % compute inverse via conjugate gradients
    tstart = tic;
    [linsolve_via_cg, toeplitz_hodlr_iters(dimstep)] = my_lin_cg(@kron_hodlr_mult,rhs,t_hr_factors,noise^2,cg_tol);
    ttemp = toc(tstart);
    toeplitz_hodlr_times(dimstep, 2) = ttemp;
    
    clear t_hr_factors t_logdets logdet_via_t_hr linsolve_via_cg
end
 
% save data for plotting later
save('data/dim_scaling_uniform.mat','default_times','toeplitz_times',...
    'rskelf_times','hodlr_times','toeplitz_hodlr_times','problem_sizes',...
    'default_iters','toeplitz_iters','rskelf_iters','hodlr_iters',...
    'toeplitz_hodlr_iters');
 
%% Helper functions
 
function K = Kmatern(x,y,maternparam)
    % x: column vector of 1D left arguments
    % y: column vector of 1D right arguments
    % maternparam: [magnitude scale, length-scale]
    % returns the 3/2 matern kernel applied to |x-y'|
    dr = abs(x-y');
    K = (maternparam(1)^2)*(1+sqrt(3)*dr/maternparam(2)).*exp(-sqrt(5)*dr/maternparam(2));
end
 
function K = Kmatern_rs(x,y,maternparam)
    % returns Kmatern but with inputs appropriate for the rskelf factorizer
    % in FLAM - specifically rskelf wants row vectors as domain inputs.
    dr = abs(x'-y);
    K = (maternparam(1)^2)*(1+sqrt(3)*dr/maternparam(2)).*exp(-sqrt(5)*dr/maternparam(2));
end
 
function A = Afun_matern(i,j,domain,maternparam)
    % i: indices of left arguments required
    % j: indices of right arguments required
    % domain: data locations
    % noise: model noise variance
    % maternparam: [magnitude scale, length-scale]
    % returns matern-covariance matrix samples for rskelf
    A = Kmatern_rs(domain(:,i),domain(:,j),maternparam);
end
 
function [x, iters] = my_lin_cg(multfun, rhs, Afactors, diag, tol)
    % Computes the linear solve (Kron(Afactors)+diag)x = rhs using CG
    % multfun: kronecker-multiplication subroutine function handle
    % rhs: column vector containing right hand side
    % Afactors: cell array containing kronecker factors
    % diag: diagonal addition to Kron(Afactors)
    % tol: error tolerance for CG iteration
    x = zeros(size(rhs));
    r = rhs;  % since our initial guess is the zero vector
    p = r;  
    iters = 0;
    while true
        iters = iters + 1;
        alpha = (r'*r)/(p'*multfun(Afactors, diag, p));
        x = x + alpha*p;
        new_r = r - alpha*multfun(Afactors, diag, p);
        if (norm(new_r)<tol) || (iters==10000)
            % returns current x
            break
        end
        beta = (new_r'*new_r)/(r'*r);
        p = new_r + beta*p;
        r = new_r;
    end
end
 
function y = kron_default_mult(factors, diag, x)
    % subroutine for computing (Kron(factors)+diag)x
    % factors: cell array containing kronecker factors
    % diag: diagonal correction
    % x: input vector
    num_factors = length(factors);
    sizes = zeros(1,num_factors);
    for fact_idx=1:num_factors
        sizes(fact_idx) = length(factors{fact_idx});
    end
    X = reshape(x,sizes);
    for fact_idx=1:num_factors
        X = refold(mtimes(factors{fact_idx},unfold(X, fact_idx, sizes)), fact_idx, sizes);
    end
    y = reshape(X, [prod(sizes), 1]);
    y = y + diag.*x;
end
 
function y = kron_toeplitz_mult(factors, diag, x)
    % subroutine for computing (Kron(factors)+diag)x
    % factors: cell array containing toeplitz kronecker factors
    % diag: diagonal correction
    % x: input vector
    num_factors = length(factors);
    sizes = zeros(1, num_factors);
    for fact_idx = 1:num_factors
        sizes(fact_idx) = length(factors{fact_idx});
    end
    X = reshape(x, sizes);
    for fact_idx=1:num_factors
        X = refold(toeplitzmult(factors{fact_idx}, unfold(X, fact_idx, sizes)), fact_idx, sizes);
    end
    y = reshape(X, [prod(sizes), 1]);
    y = y + diag.*x;
end
 
function y = kron_rskelf_mult(factors, diag, x)
    % subroutine for computing (Kron(factors)+diag)x
    % factors: cell array containing rskelf-decomposed kronecker factors
    % diag: diagonal correction
    % x: input vector
    num_factors = length(factors);
    sizes = zeros(1,num_factors);
    for fact_idx=1:num_factors
        sizes(fact_idx) = factors{fact_idx}.N;
    end
    X = reshape(x,sizes);
    for fact_idx=1:num_factors
        X = refold(rskelf_mv(factors{fact_idx}, unfold(X, fact_idx, sizes)), fact_idx, sizes);
    end
    y = reshape(X, [prod(sizes), 1]);
    y = y + diag.*x;
end
 
function y = kron_hodlr_mult(factors, diag, x)
    % subroutine for computing (Kron(factors)+diag)x
    % factors: cell array containing kronecker factors in HODLR form
    % diag: diagonal correction
    % x: input vector
    num_factors = length(factors);
    sizes = zeros(1,num_factors);
    for fact_idx=1:num_factors
        foo = factors{fact_idx}.sz;
        sizes(fact_idx) = foo(1);
    end
    X = reshape(x,sizes);
    for fact_idx=1:num_factors
        X = refold(mtimes(factors{fact_idx},unfold(X, fact_idx, sizes)), fact_idx, sizes);
    end
    y = reshape(X, [prod(sizes), 1]);
    y = y + diag.*x;
end
 
function correction = detlemma_correction_rskelf(factors,noise)
    % computes logdet(Kron(factors)+noise*I) - logdet(Kron(factors))
    % using the determinant lemma for rank-1 updates and kronecker
    % structure
    % factors: kronecker factors of primary matrix, rskelf form
    % noise: coefficient on identity correction
    num_factors = length(factors);
    full_diag = 1;
    for fact_idx=1:num_factors
        full_diag = kron(full_diag,rskelf_diag(factors{fact_idx},1));
    end
    correction = sum(1+noise*full_diag);
end

function correction = detlemma_correction_hodlr(factors,noise)
    % computes logdet(Kron(factors)+noise*I) - logdet(Kron(factors))
    % using the determinant lemma for rank-1 updates and kronecker
    % structure
    % factors: kronecker factors of primary matrix, hodlr form
    % noise: coefficient on identity correction
    num_factors = length(factors);
    full_diag = 1;
    for fact_idx=1:num_factors
        full_diag = kron(full_diag,diag(inv(factors{fact_idx})));
    end
    correction = sum(1+noise*full_diag);
end
 
function y = toeplitzmult(A, x)
    % uses DFT to quickly multiply columns of x by toeplitz A
     [n, num_vecs] = size(x);
     column = [A(:,1); 0; flipud(A(2:end,1))];
     p = ifft(fft(column).*fft([x; zeros(n,num_vecs)]));
     y = p(1:n,:);
end
 
function matrixout = refold(vectorized, axis, sizes)
    % reshapes vectorized according to sizes along index given by axis
    % part of the Kronecker mat-vec and inverse machinery
    if axis==1
        matrixout = reshape(vectorized, sizes);
    else
        num_dims = length(sizes);
        temp = sizes(axis);
        sizes(axis) = [];
        matrixout = reshape(vectorized, [temp sizes]);
        matrixout = permute(matrixout, [2:axis, 1, (axis+1):num_dims]);
    end
end
 
function vectorized = unfold(matrixin, axis, sizes)
    % reshapes matrixin according to sizes along index given by axis
    % part of the Kronecker mat-vec and inverse machinery
    if axis == 1
        vectorized = reshape(matrixin, sizes(1), []);
    else
        num_dims = length(sizes);
        matrixin = permute(matrixin, [axis, 1:(axis-1), (axis+1):num_dims]);
        vectorized = reshape(matrixin, sizes(axis), []);
    end
end
