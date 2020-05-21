%% Numerical Experiment: fixed dataset scale, nonuniform grid, dimensionality scaling
% This experiment probes the scaling of different SKI paradigms 
% (RS-factorization exploiting, HODLR-factorization exploiting, and default
% SKI) for a kernel with axis-multiplicative structure as the problem 
% dimensionality grows, for a fixed dataset scale.

% This experiment is complementary to the experiment "fixed dimensionality,
% nonuniform grid, dataset scaling", which simply scales problem
% dimensionality with a fixed number of data points per axis instead of
% scaling dataset size with fixed dimensionality. That experiment is
% located in "experiment_data_scaling_nonuniform.m"

% Due to the non-uniform grid, Toeplitz structure is not achievable even
% with stationary kernels, therefore we use a non-stationary kernel - the
% Matern 3/2 kernel multiplied by a polynomial kernel. Additionally, 
% because the grid need not be uniform we can probe scaling behavior in a 
% case well-suited to rank-structured factorizations, namely a domain 
% comprised of two well-separated subdomains, [-1,-0.7]^3 and [0.7,1]^3.
 
%% Initialization
 
% problem dimensionalities
% 1D case means trivial Kronecker product, which changes some helper
% functions, so we omit.
min_num_dims = 2; 
max_num_dims = 5;
 
% dataset size parameters:
% each dimension will have n grid locations, 
% where n is in [data_low, data_high]
half_data_low = 5; half_data_high = 7;
 
% kernel parameters: Matern
maternparam = [1.2, 2.3];
% kernel parameters: Polynomial
polyparam = [0.8, 0.1, 2.0];
% noise variance:
noise = 1e-2;
 
% conjugate gradients error tolerance
cg_tol = 1e-4;
 
% RS factorization tolerance
rs_tol = 1e-4;
 
% containers for output information:
% timing
% default: 
% (:,1)-total time needed for solve without linsolve
% (:,2)-total time needed for linsolve via eigendecomp
% (:,3)-total time needed for linsolve via CG
default_times = zeros(max_num_dims, 3);  
% the rest skip the eigendecomp, so no linear solve shortcuts (CG only)
% therefore: 
% (:,1)-total time needed for solve without linsolve
% (:,2)-total time neede for linsolve via CG
rskelf_times = zeros(max_num_dims, 2);
hodlr_times = zeros(max_num_dims, 2);
% problem sizes
problem_sizes = zeros(max_num_dims, 1);
% linear CG iterations, for reference
default_iters = zeros(max_num_dims, 1);
rskelf_iters = zeros(max_num_dims, 1);
hodlr_iters = zeros(max_num_dims, 1);
 
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
        foo = randsample(half_data_low:half_data_high,1);
        bar = randsample(half_data_low:half_data_high,1);
        size_counter(axis) = foo+bar;
        domain{axis} = [1-0.3*rand(foo,1); -1+0.3*rand(bar,1)];
    end
    problem_sizes(dimstep) = prod(size_counter);
    
    % Additional setup: generate a right-hand-side for mat-vecs and linear
    % solves. The values are arbitrary.
    rhs = -2.5 + 5*rand(problem_sizes(dimstep), 1);
    
    % generate component kernel matrices for Default, Toeplitz, HODLR, and
    % HODLR-Toeplitz methodologies
    K_cells = cell(1,num_dims);
    for axis=1:num_dims
        K_cells{axis} = Knonstat(domain{axis}, domain{axis}, maternparam, polyparam);
    end
    
    fprintf('Beginning dimension-%i experiments with %i points\n', [num_dims, problem_sizes(dimstep)])
    
    %%%%%%%%%%%
    % SECTION % Default SKI, Toeplitz SKI
    %%%%%%%%%%%
    
    fprintf('Default SKI\n')
    
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
        Afunwrapper = @(i,j)Afun_nonstat(i,j,domain{axis}',maternparam,polyparam);
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
    rskelf_times(dimstep,2) = ttemp;
    
    clear rs_factors logdets Afunwrapper logdet_via_rs linsolve_via_cg
    
    %%%%%%%%%%%
    % SECTION % HODLR SKI
    %%%%%%%%%%%
    
    fprintf('HODLR SKI\n')
    
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
end
 
% save data for plotting later
save('data/dim_scaling_nonuniform.mat','default_times','rskelf_times',...
    'hodlr_times','problem_sizes','default_iters','rskelf_iters',...
    'hodlr_iters');
 
%% Helper functions

function K = Kmatern(x,y,maternparam)
    % x: column vector of 1D left arguments
    % y: column vector of 1D right arguments
    % maternparam: [magnitude scale, length-scale]
    % returns the 3/2 matern kernel applied to |x-y'|
    dr = abs(x-y');
    K = (maternparam(1)^2)*(1+sqrt(3)*dr/maternparam(2)).*exp(-sqrt(5)*dr/maternparam(2));
end

function K = Knonstat(x,y,maternparam,polyparam)
    % x: column vector of 1D left arguments
    % y: column vector of 1D right arguments
    % maternparam: [magnitude scale, length-scale]
    % polyparam: [coefficient, shift, exponent]
    % returns the elementwise product of Kmatern and the polynomial kernel
    K = ((polyparam(1)*x*y'+polyparam(2)).^polyparam(3)).*Kmatern(x,y,maternparam);
end

function K = Kmatern_rs(x,y,maternparam)
    % returns Kmatern but with inputs appropriate for the rskelf factorizer
    % in FLAM - specifically rskelf wants row vectors as domain inputs.
    dr = abs(x'-y);
    K = (maternparam(1)^2)*(1+sqrt(3)*dr/maternparam(2)).*exp(-sqrt(5)*dr/maternparam(2));
end

function K = Knonstat_rs(x,y,maternparam,polyparam)
    % returns Knonstat but with inputs appropriate for the rskelf
    % factorizer in FLAM - specifically rskelf wants row vectors as domain
    % inputs.
    K = ((polyparam(1)*x'*y+polyparam(2)).^polyparam(3)).*Kmatern_rs(x,y,maternparam);
end

function A = Afun_nonstat(i,j,domain,maternparam,polyparam)
    % i: indices of left arguments required
    % j: indices of right arguments required
    % domain: data locations
    % noise: model noise variance
    % maternparam: [magnitude scale, length-scale]
    % returns polynomial*matern-covariance matrix samples for rskelf
    A = Knonstat_rs(domain(:,i),domain(:,j),maternparam,polyparam);
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
    correction = sum(reallog(1+noise*full_diag));
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
    correction = sum(reallog(1+noise*full_diag));
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