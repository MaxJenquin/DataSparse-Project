%% Auxiliary Numerical Experiment: fixed dataset scale, uniform grid, 
% dimensionality scaling, no CG iteration

% This sub-experiment mirrors experiment_dim_scaling_uniform.m, but with no
% CG iteration to slow it down and prevent it from reaching
% high-dimensional results. Its counterpart in the nonuniform case is
% likewise the mirror of experiment_dim_scaling_nonuniform.m, and can be
% found in subexperiment_dim_scaling_nonuniform.m
 
%% Initialization
 
% problem dimensionalities
min_num_dims = 2; % 1D case means trivial Kronecker product, which changes some helper functions
max_num_dims = 8;
 
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
default_times = zeros(max_num_dims, 1);  
toeplitz_times = zeros(max_num_dims, 1);
rskelf_times = zeros(max_num_dims, 1);
hodlr_times = zeros(max_num_dims, 1);
toeplitz_hodlr_times = zeros(max_num_dims, 1);
% problem sizes
problem_sizes = zeros(max_num_dims, 1);
 
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
    default_times(dimstep) = ttemp;
    toeplitz_times(dimstep) = ttemp;
    
    clear eig_decomp vecs_decomp vecs_decomp_trps full_eigs foo logdet_via_eigs
    
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
    rskelf_times(dimstep) = ttemp;
    
    clear rs_factors logdets Afunwrapper logdet_via_rs 
    
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
    hodlr_times(dimstep) = ttemp;
    
    clear hr_factors logdets logdet_via_hr
    
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
    toeplitz_hodlr_times(dimstep) = ttemp;
    
    clear t_hr_factors t_logdets logdet_via_t_hr
end
 
% save data for plotting later
save('data/dim_scaling_uniform_extra.mat','default_times','toeplitz_times',...
    'rskelf_times','hodlr_times','toeplitz_hodlr_times','problem_sizes'...
    );
 
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
