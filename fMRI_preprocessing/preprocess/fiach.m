function results = fiach(input, t, tr, varargin)
%FIACH MATLAB/SPM implementation of the FIACH R pipeline.
%
%   RESULTS = FIACH(INPUT, T, TR) applies the retrospective noise-control
%   workflow from tierneytim/FIACH to realigned fMRI images. INPUT is either
%   one 4-D NIfTI/Analyze image or a cell array of 3-D images. T is the
%   deterministic outlier threshold in percent signal change and TR is in
%   seconds.
%
%   Name-value options mirror the R function where possible:
%     'RP'           realignment-parameter text file (default: [])
%     'MaxGap'       longest run spline-corrected before median scrubbing (1)
%     'Freq'         DCT high-pass period in seconds (128)
%     'NMads'        stochastic outlier threshold in scaled MADs (1.96)
%     'Mask'          optional aligned brain-mask path or 3-D mask array ([])
%     'UseUserMask'   use Mask; automatic when Mask is supplied ([])
%     'DefaultMask'  use the FIACH k-means brain mask (true)
%     'QuantileMask' fallback-mask multiplier of the mean median signal (.7)
%     'GMMMethod'    'R' (default) or 'ModifiedICH'
%     'GMMPosteriorThreshold'  posterior threshold for 'ModifiedICH' (.5)
%     'Overwrite'    allow replacement of existing output files (false)
%     'MakePlots'    save a non-interactive GMM diagnostic PNG (false)
%
%   Output files follow the R package convention as closely as SPM permits:
%     filt_<input>.nii                 corrected functional data
%     <input-folder>_fiach_diagnostics/mask.nii
%     <input-folder>_fiach_diagnostics/rtsnr.nii
%     <input-folder>_fiach_diagnostics/median.nii and noise_mask.nii
%     <input-folder>_fiach_diagnostics/metrics.txt
%     noise_basis6.txt and (when RP is supplied) gs.txt and fd_noise.txt
%
%   Notes
%   -----
%   This is a direct MATLAB translation of the main R function and its
%   private helpers, not a wrapper around run_fiach_DC_modified_ICH.m.
%   It deliberately uses FIACH's DCT implementation rather than SPM's
%   spm_filter so that the high-pass step matches R/highPass.R. SPM12 must
%   be on the MATLAB path for image I/O. SPM writes uncompressed NIfTI
%   output even when the input is .nii.gz.
%
%   Upstream reference: https://github.com/tierneytim/FIACH

    fiach_require_spm();
    opts = fiach_parse_options(varargin{:});
    inputFiles = fiach_normalise_input(input);

    validateattributes(t,  {'numeric'}, {'scalar', 'real', 'finite', '>=', 0}, mfilename, 't', 2);
    validateattributes(tr, {'numeric'}, {'scalar', 'real', 'finite', '>', 0},  mfilename, 'tr', 3);
    validateattributes(opts.MaxGap, {'numeric'}, {'scalar', 'integer', '>=', 1});
    validateattributes(opts.Freq, {'numeric'}, {'scalar', 'real', 'finite', '>', 0});
    validateattributes(opts.NMads, {'numeric'}, {'scalar', 'real', 'finite', '>=', 0});

    for iFile = 1:numel(inputFiles)
        sourceFile = fiach_strip_volume_index(inputFiles{iFile});
        if exist(sourceFile, 'file') ~= 2
            error('fiach:InputNotFound', 'Functional image not found: %s', sourceFile);
        end
    end

    % R's multi-file mode is a sequence of 3-D images. Supporting a cell
    % array of 4-D files would make the time/order correspondence ambiguous.
    V = spm_vol(char(inputFiles));
    if isempty(V)
        error('fiach:NoVolumes', 'SPM did not find any readable image volumes.');
    end
    if numel(inputFiles) > 1 && numel(V) ~= numel(inputFiles)
        error('fiach:Multi4DInput', ['For multiple inputs, provide one 3-D image per file. ' ...
            'For a 4-D run, provide one file path.']);
    end

    data = double(spm_read_vols(V));
    if ndims(data) < 4
        data = reshape(data, [size(data, 1), size(data, 2), size(data, 3), 1]);
    end
    data(~isfinite(data)) = 0; % FIACH::zeroNa

    spatialSize = size(data);
    spatialSize = spatialSize(1:3);
    nTime = size(data, 4);
    nVoxels = prod(spatialSize);
    mat = reshape(data, nVoxels, nTime).'; % time x voxel, as in R::arrMat

    [inputDir, ~, outputExt] = fiach_file_parts(inputFiles{1});
    diagnosticDir = [inputDir '_fiach_diagnostics'];
    if exist(diagnosticDir, 'dir') ~= 7
        mkdir(diagnosticDir);
    end

    % ---- FIACH mask creation ------------------------------------------------
    voxelMedians = median(mat, 1);
    if opts.UseUserMask
        [maskVec, maskSource] = fiach_resolve_user_mask(opts.Mask, spatialSize, V(1));
    else
        maskVec = [];
        maskSource = '';
    end
    if isempty(maskVec)
        [maskVec, maskFit] = fiach_kmeans_mask(voxelMedians);
        if maskFit < 0.8 || ~opts.DefaultMask
            warning('fiach:FallbackMask', ['Using FIACH''s median-signal fallback mask ' ...
                '(k-means fit = %.3f). Check for receive-field bias.'], maskFit);
            maskVec = voxelMedians > mean(voxelMedians) * opts.QuantileMask;
            maskSource = 'FIACH fallback mask';
        else
            maskSource = 'FIACH k-means mask';
        end
    else
        maskFit = NaN;
    end
    maskVec = logical(maskVec(:)).';
    if ~any(maskVec)
        error('fiach:EmptyMask', 'The FIACH brain mask contains no voxels.');
    end
    smallBrain = mat(:, maskVec);

    % ---- High-pass filter: exact DCT construction in R/highPass.R ----------
    smallHp = fiach_high_pass(smallBrain, opts.Freq, tr);
    hpMat = mat;                     % R preserves out-of-mask data unchanged.
    hpMat(:, maskVec) = smallHp;

    % ---- Robust tSNR ---------------------------------------------------------
    brainMedians = voxelMedians(maskVec);          % median BEFORE high-pass
    brainMads = fiach_col_mad(smallHp);            % scaled MAD AFTER high-pass
    robustTsnr = brainMedians ./ brainMads;
    finiteTsnr = robustTsnr(isfinite(robustTsnr));
    if numel(finiteTsnr) < 10
        error('fiach:InsufficientTSNR', 'Fewer than 10 finite brain-voxel tSNR values were found.');
    end

    gmmMethod = fiach_normalise_gmm_method(opts.GMMMethod);
    switch gmmMethod
        case 'R'
            [mixture, noiseThreshold] = fiach_noise_segmentation(finiteTsnr);
            noisyInBrain = robustTsnr < noiseThreshold;
            noisyInBrain(~isfinite(robustTsnr)) = false;
            applySpatialNoiseExpansion = true;
        case 'ModifiedICH'
            [mixture, noisyInBrain] = fiach_modified_ich_noise_segmentation( ...
                robustTsnr, opts.GMMPosteriorThreshold);
            noiseThreshold = NaN; % This method classifies by posterior probability.
            % The selector changes GMM fitting/segmentation only. Retain the
            % R FIACH spatial expansion so downstream PCA remains comparable.
            applySpatialNoiseExpansion = true;
        otherwise
            error('fiach:InvalidGMMMethod', 'Unsupported GMM method: %s', gmmMethod);
    end
    noiseMaskVec = false(1, nVoxels);
    noiseMaskVec(maskVec) = noisyInBrain;
    noiseMask = reshape(noiseMaskVec, spatialSize);
    if applySpatialNoiseExpansion
        noiseMask = fiach_expand_noise_mask(noiseMask);
    end
    noiseMaskVec = logical(noiseMask(:)).';

    % ---- Noise regressors: scale noisy voxel series, then PCA --------------
    noisyWithinBrain = noiseMaskVec(maskVec);
    noise = smallHp(:, noisyWithinBrain);
    noiseRegs = fiach_noise_pca(noise, 6);

    % ---- FIACH bad-data detection and correction ----------------------------
    naBrain = fiach_bad_data(smallHp, brainMedians, brainMads, opts.NMads, t);
    nInitiallyFlagged = nnz(isnan(naBrain));

    % The R code replaces bad first/last observations before spline fitting.
    firstBad = isnan(naBrain(1, :));
    lastBad = isnan(naBrain(end, :));
    naBrain(1, firstBad) = brainMedians(firstBad);
    naBrain(end, lastBad) = brainMedians(lastBad);

    filteredBrain = fiach_na_spline(naBrain, opts.MaxGap);
    filteredBrain = fiach_replace_remaining_nans(filteredBrain, brainMedians);
    outputMat = hpMat;
    outputMat(:, maskVec) = filteredBrain;
    correctedData = reshape(outputMat.', [spatialSize, nTime]);

    % ---- Write FIACH-compatible derivative files ----------------------------
    outputFiles = fiach_write_functional(correctedData, V, inputFiles, outputExt, opts.Overwrite);
    diagFiles = fiach_write_diagnostics(diagnosticDir, V(1), maskVec, noiseMaskVec, ...
        robustTsnr, voxelMedians, spatialSize, opts.Overwrite);

    denom = nTime * (size(smallBrain, 2) - nnz(noisyWithinBrain));
    if denom > 0
        percentChanged = 100 * nInitiallyFlagged / denom;
    else
        percentChanged = NaN;
    end
    peakTsnr = mixture.peakTsnr;
    metricFile = fullfile(diagnosticDir, 'metrics.txt');
    fiach_write_metrics(metricFile, percentChanged, peakTsnr, opts.Overwrite);

    % The R package writes its first six PC scores under this conventional name.
    noiseFile = fullfile(inputDir, 'noise_basis6.txt');
    fiach_write_matrix(noiseFile, noiseRegs, opts.Overwrite);

    gsFile = '';
    fdNoiseFile = '';
    if ~isempty(opts.RP)
        rp = fiach_read_motion_parameters(opts.RP, nTime);
        fd = fiach_framewise_displacement(rp);
        globalSignal = mean(fiach_zero_nan((filteredBrain - brainMedians) ./ brainMads), 2);
        gsFile = fullfile(inputDir, 'gs.txt');
        fiach_write_matrix(gsFile, globalSignal, opts.Overwrite);
        fdNoiseFile = fullfile(inputDir, 'fd_noise.txt');
        fiach_write_matrix(fdNoiseFile, [fd, noiseRegs], opts.Overwrite);

        % In the R implementation noise_basis6 contains motion plus PCs when
        % an RP file is supplied. Replace the PCs-only version accordingly.
        fiach_write_matrix(noiseFile, [rp, noiseRegs], true);
    end

    if opts.MakePlots
        fiach_write_mixture_plot(diagnosticDir, finiteTsnr, mixture, noiseThreshold);
    end

    results = struct();
    results.input = inputFiles;
    results.functionalFiles = outputFiles;
    results.diagnosticDirectory = diagnosticDir;
    results.diagnosticFiles = diagFiles;
    results.noiseBasisFile = noiseFile;
    results.globalSignalFile = gsFile;
    results.fdNoiseFile = fdNoiseFile;
    results.metricFile = metricFile;
    results.maskFit = maskFit;
    results.maskSource = maskSource;
    results.useUserMask = opts.UseUserMask;
    results.gmmMethod = gmmMethod;
    results.noiseMaskSpatialExpansion = applySpatialNoiseExpansion;
    results.noiseThreshold = noiseThreshold;
    if strcmp(gmmMethod, 'ModifiedICH')
        results.noisePosteriorThreshold = opts.GMMPosteriorThreshold;
    else
        results.noisePosteriorThreshold = NaN;
    end
    results.noiseRegressors = noiseRegs;
    results.metrics = struct('percentDataChanged', percentChanged, 'peakTSNR', peakTsnr, ...
        'nInitiallyFlagged', nInitiallyFlagged);
    results.mixture = mixture;

    fprintf(['FIACH complete: %.3f%% of eligible samples initially flagged; ' ...
        '%d noise regressors written.\n'], percentChanged, size(noiseRegs, 2));
end


function opts = fiach_parse_options(varargin)
    parser = inputParser;
    parser.FunctionName = 'fiach';
    addParameter(parser, 'RP', [], @(x) isempty(x) || ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'MaxGap', 1, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'Freq', 128, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'NMads', 1.96, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'Mask', [], @(x) isempty(x) || ischar(x) || ...
        (isstring(x) && isscalar(x)) || isnumeric(x) || islogical(x));
    addParameter(parser, 'UseUserMask', [], @(x) isempty(x) || (islogical(x) && isscalar(x)));
    addParameter(parser, 'DefaultMask', true, @(x) islogical(x) && isscalar(x));
    addParameter(parser, 'QuantileMask', .7, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    addParameter(parser, 'GMMMethod', 'R', @(x) ischar(x) || (isstring(x) && isscalar(x)));
    addParameter(parser, 'GMMPosteriorThreshold', .5, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
    addParameter(parser, 'Overwrite', false, @(x) islogical(x) && isscalar(x));
    addParameter(parser, 'MakePlots', false, @(x) islogical(x) && isscalar(x));
    parse(parser, varargin{:});
    opts = parser.Results;
    if isstring(opts.RP)
        opts.RP = char(opts.RP);
    end
    if isstring(opts.Mask)
        opts.Mask = char(opts.Mask);
    end
    if isstring(opts.GMMMethod)
        opts.GMMMethod = char(opts.GMMMethod);
    end
    % Supplying a mask is normally all a caller needs to do. The explicit
    % flag also lets a batch script retain a mask path while choosing the
    % built-in FIACH mask for a particular run.
    if isempty(opts.UseUserMask)
        opts.UseUserMask = ~isempty(opts.Mask);
    end
    if opts.UseUserMask && isempty(opts.Mask)
        error('fiach:MissingUserMask', ...
            'UseUserMask is true, but no Mask was supplied.');
    end
end


function fiach_require_spm()
    required = {'spm_vol', 'spm_read_vols', 'spm_create_vol', 'spm_write_vol', 'spm_type'};
    for i = 1:numel(required)
        if exist(required{i}, 'file') ~= 2
            error('fiach:SPMNotFound', ['SPM12 is required for NIfTI I/O. Add the SPM12 ' ...
                'directory to the MATLAB path before calling fiach. Missing: %s'], required{i});
        end
    end
end


function files = fiach_normalise_input(input)
    if isstring(input)
        files = cellstr(input(:));
    elseif ischar(input)
        files = cellstr(input);
    elseif iscell(input) && all(cellfun(@(x) ischar(x) || (isstring(x) && isscalar(x)), input))
        files = cellfun(@char, input(:), 'UniformOutput', false);
    else
        error('fiach:InvalidInput', 'INPUT must be a path, a character matrix, or a cell array of paths.');
    end
    files = cellfun(@strtrim, files, 'UniformOutput', false);
    files = files(~cellfun(@isempty, files));
    if isempty(files)
        error('fiach:EmptyInput', 'At least one functional image must be supplied.');
    end
end


function pathWithoutIndex = fiach_strip_volume_index(pathWithIndex)
    token = regexp(pathWithIndex, '^(.*),[0-9]+$', 'tokens', 'once');
    if isempty(token)
        pathWithoutIndex = pathWithIndex;
    else
        pathWithoutIndex = token{1};
    end
end


function [folder, base, extension] = fiach_file_parts(inputPath)
    inputPath = fiach_strip_volume_index(inputPath);
    [folder, name, extension] = fileparts(inputPath);
    if isempty(folder)
        folder = pwd;
    end
    if strcmpi(extension, '.gz')
        [~, base, innerExtension] = fileparts(name);
        if strcmpi(innerExtension, '.nii')
            extension = '.nii';
        else
            extension = innerExtension;
        end
    else
        base = name;
    end
    if strcmpi(extension, '.hdr')
        extension = '.img';
    end
    if isempty(extension)
        extension = '.nii';
    end
end


function [maskVec, source] = fiach_resolve_user_mask(maskOption, spatialSize, Vref)
% Resolve a user-supplied mask without silently reslicing it. A supplied
% mask must already be in the functional-image grid and orientation.
    maskVec = [];
    source = '';
    if isempty(maskOption)
        return;
    end

    if ischar(maskOption) || (isstring(maskOption) && isscalar(maskOption))
        maskPath = fiach_strip_volume_index(char(maskOption));
        if exist(maskPath, 'file') ~= 2
            error('fiach:MaskNotFound', 'User mask not found: %s', maskPath);
        end
        Vmask = spm_vol(maskPath);
        if numel(Vmask) ~= 1
            error('fiach:MaskMustBe3D', 'User mask must contain exactly one 3-D volume: %s', maskPath);
        end
        if ~isequal(Vmask.dim(1:3), Vref.dim(1:3))
            error('fiach:MaskDimensionMismatch', ['User-mask dimensions [%s] do not match ' ...
                'functional dimensions [%s]. Reslice the mask to functional space first.'], ...
                num2str(Vmask.dim(1:3)), num2str(Vref.dim(1:3)));
        end
        if isfield(Vmask, 'mat') && isfield(Vref, 'mat') && ...
                max(abs(Vmask.mat(:) - Vref.mat(:))) > 1e-4
            error('fiach:MaskGeometryMismatch', ['User-mask affine does not match the functional image. ' ...
                'Coregister/reslice it to functional space before running FIACH.']);
        end
        maskData = spm_read_vols(Vmask);
        source = maskPath;
    else
        maskData = maskOption;
        source = 'user-supplied numeric mask';
    end

    if ~isequal(size(maskData), spatialSize)
        error('fiach:MaskDimensionMismatch', ['User-mask dimensions [%s] do not match ' ...
            'functional dimensions [%s].'], num2str(size(maskData)), num2str(spatialSize));
    end
    maskVec = logical(isfinite(maskData) & maskData ~= 0);
    maskVec = maskVec(:).';
    if ~any(maskVec)
        error('fiach:EmptyUserMask', 'The supplied user mask contains no non-zero finite voxels.');
    end
end


function filtered = fiach_high_pass(x, frequency, tr)
% Direct port of R/highPass.R. Rows are time points; columns are voxels.
    nTime = size(x, 1);
    % R rounds .5 ties to even; MATLAB rounds them away from zero.
    nBasis = fiach_round_like_r(2 * (nTime * tr / (frequency + 1)));
    if nBasis < 2
        filtered = x;
        return;
    end
    n = (0:nTime-1).';
    basis = zeros(nTime, nBasis - 1);
    for i = 2:nBasis
        basis(:, i - 1) = sqrt(2 / nTime) * cos(pi * (2 * n + 1) * (i - 1) / (2 * nTime));
    end
    filtered = x - basis * (basis.' * x);
end


function madValues = fiach_col_mad(x)
% FIACH C++ colMad: median(abs(x - median(x))) * 1.4826.
    centres = median(x, 1);
    madValues = median(abs(x - centres), 1) * 1.4826;
end


function value = fiach_round_like_r(value)
% R's round() uses round-to-even. FIACH's basis count is always non-negative.
    lower = floor(value);
    fraction = value - lower;
    tolerance = 4 * eps(max(1, abs(value)));
    if abs(fraction - .5) <= tolerance
        value = lower + mod(lower, 2);
    else
        value = floor(value + .5);
    end
end


function [mask, fit] = fiach_kmeans_mask(x)
% Direct deterministic 1-D equivalent of R/kmeansMask.R.
    samples = fiach_quantile(x, linspace(0, 1, 1000));
    [labels, centres] = fiach_kmeans_1d(samples);
    transition = find(abs(diff(labels)) > 0, 1, 'first');
    if isempty(transition)
        threshold = centres(1);
    else
        threshold = samples(transition);
    end
    mask = x > threshold;
    totalSS = sum((samples - mean(samples)).^2);
    withinSS = sum((samples(labels == 1) - centres(1)).^2) + ...
        sum((samples(labels == 2) - centres(2)).^2);
    if totalSS == 0
        fit = 0;
    else
        fit = max(0, (totalSS - withinSS) / totalSS);
    end
end


function [labels, centres] = fiach_kmeans_1d(x)
% Lloyd's method avoids making the port depend on the Statistics Toolbox.
    x = x(:);
    centres = [fiach_quantile(x, .25); fiach_quantile(x, .75)];
    if centres(1) == centres(2)
        labels = ones(size(x));
        return;
    end
    labels = zeros(size(x));
    for iter = 1:200
        distances = [abs(x - centres(1)), abs(x - centres(2))];
        [~, newLabels] = min(distances, [], 2);
        newCentres = centres;
        for k = 1:2
            if any(newLabels == k)
                newCentres(k) = mean(x(newLabels == k));
            end
        end
        if isequal(newLabels, labels) && max(abs(newCentres - centres)) < eps(max(abs(centres)) + 1)
            labels = newLabels;
            centres = newCentres;
            break;
        end
        labels = newLabels;
        centres = newCentres;
    end
    [centres, order] = sort(centres);
    relabel = zeros(2, 1);
    relabel(order) = 1:2;
    labels = relabel(labels);
end


function method = fiach_normalise_gmm_method(value)
    value = lower(strtrim(char(value)));
    switch value
        case {'r', 'rcode', 'r-compatible', 'original'}
            method = 'R';
        case {'modifiedich', 'modified_ich', 'ich', 'fiach_modified_ich'}
            method = 'ModifiedICH';
        otherwise
            error('fiach:InvalidGMMMethod', ['GMMMethod must be ''R'' or ''ModifiedICH''. ' ...
                'Received: %s'], char(value));
    end
end


function [model, threshold] = fiach_noise_segmentation(rtsnr)
% MATLAB implementation of .noiseSeg, .quantPriors and .tooClose in fiach.R.
    rtsnr = rtsnr(isfinite(rtsnr));
    resampled = fiach_quantile(rtsnr, linspace(0, .98, 1000));
    proportions = .01:.05:.95;
    models = repmat(struct('mu', [], 'sigma', [], 'lambda', [], 'LL', [], ...
        'niter', [], 'valid', false), 1, numel(proportions));

    for i = 1:numel(proportions)
        p = proportions(i);
        cut = fiach_quantile(resampled, p);
        lower = resampled(resampled <= cut);
        upper = resampled(resampled >= cut);
        if numel(lower) < 2 || numel(upper) < 2
            continue;
        end
        initial = struct('mu', [mean(lower), mean(upper)], ...
            'sigma', [std(lower, 0), std(upper, 0)], 'lambda', [p, 1-p]);
        if any(~isfinite([initial.mu, initial.sigma])) || any(initial.sigma <= 0)
            continue;
        end
        candidate = fiach_em_gmm_1d(resampled, initial);
        candidate.valid = ~fiach_gmm_too_close(candidate);
        models(i) = candidate;
    end

    usable = models([models.valid]);
    if isempty(usable)
        threshold = fiach_quantile(rtsnr, .05);
        model = struct('mu', NaN, 'sigma', NaN, 'lambda', NaN, 'LL', NaN, ...
            'niter', 0, 'fallback', true, 'peakTsnr', NaN);
        return;
    end

    score = zeros(1, numel(usable));
    for i = 1:numel(usable)
        score(i) = usable(i).LL * abs(diff(sort(usable(i).mu)));
    end
    [~, bestIndex] = min(score); % This reproduces which.min(LL * muDiff) in R.
    model = usable(bestIndex);
    [~, largest] = max(model.lambda);
    threshold = fiach_normal_quantile(.05, model.mu(largest), model.sigma(largest));
    model.fallback = false;
    model.peakTsnr = model.mu(largest);
end


function [model, noisy] = fiach_modified_ich_noise_segmentation(rtsnr, posteriorThreshold)
% Port of TSNR_gmm + fiach_segment_noisy from run_fiach_DC_modified_ICH.m.
% It uses the Statistics and Machine Learning Toolbox fitgmdist routine,
% histogram-derived starts, 98%% upper-tail trimming, and posterior-based
% allocation to the lower-mean component.
    if exist('fitgmdist', 'file') ~= 2 || exist('statset', 'file') ~= 2
        error('fiach:ICHGMMToolboxRequired', ['GMMMethod ''ModifiedICH'' requires ' ...
            'the Statistics and Machine Learning Toolbox (fitgmdist/statset).']);
    end

    fitData = rtsnr(isfinite(rtsnr) & rtsnr ~= 0);
    fitData = sort(double(fitData(:)));
    if numel(fitData) < 4
        error('fiach:ICHGMMInsufficientData', 'ModifiedICH GMM requires at least four finite, non-zero tSNR values.');
    end
    index98 = max(1, round(.98 * numel(fitData)));
    trimmed = fitData(1:index98);
    lowerBound = min(trimmed);
    upperBound = max(trimmed);
    if ~isfinite(lowerBound) || ~isfinite(upperBound) || lowerBound == upperBound
        error('fiach:ICHGMMBadRange', 'ModifiedICH GMM requires non-degenerate tSNR values.');
    end

    nBins = 200;
    edges = linspace(lowerBound, upperBound, nBins + 1);
    centres = (edges(1:end-1) + edges(2:end)) / 2;
    counts = histcounts(trimmed, edges);
    nPriors = 19;
    lowerMeans = nan(1, nPriors); upperMeans = nan(1, nPriors);
    lowerVars = nan(1, nPriors);  upperVars = nan(1, nPriors);
    for i = 1:nPriors
        cut = max(1, round(.05 * i * nBins));
        lowerWeights = counts(1:cut);
        upperWeights = counts(cut+1:end);
        if sum(lowerWeights) > 1
            lowerCentres = centres(1:cut);
            lowerMeans(i) = (lowerWeights * lowerCentres.') / sum(lowerWeights);
            lowerDelta = lowerCentres - lowerMeans(i);
            lowerVars(i) = max((lowerWeights * (lowerDelta.^2).') / sum(lowerWeights), eps);
        end
        if sum(upperWeights) > 1
            upperCentres = centres(cut+1:end);
            upperMeans(i) = (upperWeights * upperCentres.') / sum(upperWeights);
            upperDelta = upperCentres - upperMeans(i);
            upperVars(i) = max((upperWeights * (upperDelta.^2).') / sum(upperWeights), eps);
        end
    end

    bestNll = inf;
    bestFit = [];
    for i = 4:nPriors
        startValues = [lowerMeans(i), upperMeans(i), lowerVars(i), upperVars(i)];
        if any(~isfinite(startValues))
            continue;
        end
        initialMu = sort([lowerMeans(i); upperMeans(i)]);
        initialVars = sort([lowerVars(i), upperVars(i)]);
        initial = struct('mu', initialMu, 'Sigma', reshape(initialVars, [1 1 2]));
        try
            candidate = fitgmdist(trimmed, 2, 'Start', initial, ...
                'Options', statset('MaxIter', 500), 'CovarianceType', 'diagonal');
        catch
            continue;
        end
        if candidate.NegativeLogLikelihood < bestNll
            bestNll = candidate.NegativeLogLikelihood;
            bestFit = candidate;
        end
    end
    if isempty(bestFit)
        error('fiach:ICHGMMNoFit', 'ModifiedICH GMM could not find a valid two-component fit.');
    end

    originalMu = bestFit.mu(:).';
    originalSigma = sqrt(reshape(bestFit.Sigma, 1, []));
    originalLambda = bestFit.ComponentProportion(:).';
    [mu, order] = sort(originalMu);
    sigma = originalSigma(order);
    lambda = originalLambda(order);
    [~, dominant] = max(lambda);
    if isprop(bestFit, 'NumIterations')
        nIterations = bestFit.NumIterations;
    else
        nIterations = NaN;
    end
    model = struct('mu', mu, 'sigma', sigma, 'lambda', lambda, ...
        'LL', -bestFit.NegativeLogLikelihood / numel(trimmed), 'niter', nIterations, ...
        'valid', true, 'fallback', false, 'peakTsnr', mu(dominant), ...
        'method', 'ModifiedICH', 'posteriorThreshold', posteriorThreshold);

    values = rtsnr(:);
    finiteValues = isfinite(values);
    values(~finiteValues) = 0;
    probabilities = posterior(bestFit, values);
    lowerMeanComponent = order(1);
    % The legacy fiach_segment_noisy code only writes a posterior-selected
    % voxel to its mask when its tSNR is strictly positive. Keep that detail
    % here so the ModifiedICH selector is behaviourally compatible.
    noisy = probabilities(:, lowerMeanComponent) > posteriorThreshold & values > 0;
    noisy(~finiteValues) = false;
    noisy = noisy(:).';
end


function isTooClose = fiach_gmm_too_close(model)
    [~, biggestVariance] = max(model.sigma);
    [~, smallestVariance] = min(model.sigma);
    univariateLike = abs(model.mu(biggestVariance) - model.mu(smallestVariance)) < model.sigma(biggestVariance);
    tinyComponent = any(model.lambda < .01);
    isTooClose = univariateLike || tinyComponent;
end


function model = fiach_em_gmm_1d(x, initial)
% A toolbox-free two-component Gaussian EM analogue of FIACH's C++ gmm call.
    x = x(:);
    mu = initial.mu(:).';
    sigma = max(initial.sigma(:).', sqrt(eps));
    lambda = initial.lambda(:).';
    lambda = lambda / sum(lambda);
    lastLL = -inf;
    maxIter = 1000;
    for iter = 1:maxIter
        probability = zeros(numel(x), 2);
        for k = 1:2
            probability(:, k) = lambda(k) ./ (sigma(k) * sqrt(2*pi)) .* ...
                exp(-.5 * ((x - mu(k)) ./ sigma(k)).^2);
        end
        evidence = max(sum(probability, 2), realmin('double'));
        responsibilities = probability ./ evidence;
        ll = mean(log(evidence));
        counts = sum(responsibilities, 1);
        lambda = max(counts / numel(x), eps);
        lambda = lambda / sum(lambda);
        % Both components need their own weighted mean. Matrix right
        % division (/) collapses this 1-by-2 vector to a scalar; use ./.
        mu = (responsibilities.' * x).' ./ counts;
        variance = sum(responsibilities .* (x - mu).^2, 1) ./ counts;
        sigma = sqrt(max(variance, eps));
        if abs(ll - lastLL) < 1e-8
            break;
        end
        lastLL = ll;
    end
    model = struct('mu', mu, 'sigma', sigma, 'lambda', lambda, 'LL', ll, ...
        'niter', iter, 'valid', true);
end


function maskOut = fiach_expand_noise_mask(maskIn)
% Port of .nnInterp calls in fiach.R. It combines 3-D and each planar
% majority filter after a one-voxel dilation, with replicated edge values.
    maskIn = double(maskIn ~= 0);
    full = fiach_nn_majority(maskIn, true(3, 3, 3));
    xyPlane = false(3, 3, 3); xyPlane(:, :, 2) = true; % R method 2
    xzPlane = false(3, 3, 3); xzPlane(:, 2, :) = true; % R method 3
    yzPlane = false(3, 3, 3); yzPlane(2, :, :) = true; % R method 4
    maskOut = full | fiach_nn_majority(maskIn, xyPlane) | ...
        fiach_nn_majority(maskIn, xzPlane) | fiach_nn_majority(maskIn, yzPlane);
end


function output = fiach_nn_majority(mask, kernel)
    candidates = fiach_neighbour_sum(mask, true(3, 3, 3)) > 0;
    neighbourhood = fiach_neighbour_sum(mask, kernel ~= 0);
    cutoff = floor(nnz(kernel) / 2);
    output = candidates & neighbourhood > cutoff;
end


function total = fiach_neighbour_sum(mask, kernel)
% 3-D sums with replicated borders, matching the edge treatment in .nnInterp.
    total = zeros(size(mask));
    [xOffsets, yOffsets, zOffsets] = ind2sub(size(kernel), find(kernel));
    xOffsets = xOffsets - ceil(size(kernel, 1) / 2);
    yOffsets = yOffsets - ceil(size(kernel, 2) / 2);
    zOffsets = zOffsets - ceil(size(kernel, 3) / 2);
    nx = size(mask, 1); ny = size(mask, 2); nz = size(mask, 3);
    for i = 1:numel(xOffsets)
        xIndex = min(max((1:nx) + xOffsets(i), 1), nx);
        yIndex = min(max((1:ny) + yOffsets(i), 1), ny);
        zIndex = min(max((1:nz) + zOffsets(i), 1), nz);
        total = total + mask(xIndex, yIndex, zIndex);
    end
end


function scores = fiach_noise_pca(noise, nComponents)
% FIACH scales each noisy voxel before prcomp; scores = prcomp(...)$x.
    if isempty(noise)
        error('fiach:EmptyNoiseMask', 'The noise mask contains no brain voxels after spatial expansion.');
    end
    centres = mean(noise, 1);
    scales = std(noise, 0, 1); % R::scale uses the sample SD.
    keep = sum(noise, 1) > 0 & isfinite(scales) & scales > eps;
    if ~any(keep)
        error('fiach:NoPCAVoxels', 'No finite, non-constant noisy voxel series are available for PCA.');
    end
    scaled = (noise(:, keep) - centres(keep)) ./ scales(keep);
    [u, s, ~] = svd(scaled, 'econ');
    scores = u * s;
    nComponents = min([nComponents, size(scores, 1), size(scores, 2)]);
    scores = scores(:, 1:nComponents);
end


function marked = fiach_bad_data(x, medians, mads, nMads, deterministicThreshold)
% Direct translation of FIACH_src.cpp::badData.
    upper = medians + (medians / 100 * deterministicThreshold) + nMads * mads;
    lower = medians - (medians / 100 * deterministicThreshold) - nMads * mads;
    marked = x;
    outlier = x < lower | x > upper;
    marked(outlier) = NaN;
end


function corrected = fiach_na_spline(mat, maxGap)
% MATLAB analogue of R/naSpline.R. R uses FMM endpoint conditions, whereas
% interp1's spline method is not-a-knot; only short missing runs are filled.
    corrected = mat;
    time = (1:size(mat, 1)).';
    for voxel = 1:size(mat, 2)
        isMissing = isnan(corrected(:, voxel));
        if ~any(isMissing)
            continue;
        end
        transitions = diff([false; isMissing; false]);
        starts = find(transitions == 1);
        stops = find(transitions == -1) - 1;
        eligible = false(size(isMissing));
        for run = 1:numel(starts)
            if stops(run) - starts(run) + 1 <= maxGap
                eligible(starts(run):stops(run)) = true;
            end
        end
        good = ~isMissing;
        if any(eligible) && nnz(good) >= 2
            corrected(eligible, voxel) = interp1(time(good), corrected(good, voxel), ...
                time(eligible), 'spline');
        end
    end
end


function matrix = fiach_replace_remaining_nans(matrix, medians)
    for voxel = 1:size(matrix, 2)
        missing = isnan(matrix(:, voxel));
        matrix(missing, voxel) = medians(voxel);
    end
end


function values = fiach_quantile(x, probabilities)
% R's default quantile(type = 7), implemented without toolbox dependencies.
    x = sort(x(isfinite(x)));
    if isempty(x)
        values = NaN(size(probabilities));
        return;
    end
    probabilities = min(max(probabilities, 0), 1);
    positions = 1 + (numel(x) - 1) * probabilities;
    lower = floor(positions);
    upper = ceil(positions);
    weight = positions - lower;
    values = x(lower) .* (1 - weight) + x(upper) .* weight;
end


function q = fiach_normal_quantile(p, mu, sigma)
% norminv implemented using erfcinv, avoiding a Statistics Toolbox requirement.
    q = mu - sqrt(2) * sigma * erfcinv(2 * p);
end


function outputFiles = fiach_write_functional(data, V, inputFiles, extension, overwrite)
    if numel(inputFiles) == 1
        [folder, base, ~] = fiach_file_parts(inputFiles{1});
        outputFiles = {fullfile(folder, ['filt_' base extension])};
        fiach_write_4d(data, V(1), outputFiles{1}, overwrite);
    else
        outputFiles = cell(numel(inputFiles), 1);
        for timePoint = 1:numel(inputFiles)
            [folder, base, thisExtension] = fiach_file_parts(inputFiles{timePoint});
            outputFiles{timePoint} = fullfile(folder, ['filt_' base thisExtension]);
            fiach_write_3d(data(:, :, :, timePoint), V(timePoint), outputFiles{timePoint}, overwrite, 'native');
        end
    end
end


function files = fiach_write_diagnostics(folder, Vref, maskVec, noiseMaskVec, rtsnr, medians, spatialSize, overwrite)
    mask = reshape(uint8(maskVec), spatialSize);
    noiseMask = reshape(uint8(noiseMaskVec), spatialSize);
    rtsnrMap = zeros(spatialSize);
    medianMap = zeros(spatialSize);
    rtsnrMap(maskVec) = rtsnr;
    medianMap(maskVec) = medians(maskVec);
    rtsnrMap(~isfinite(rtsnrMap)) = 0;

    files = struct();
    files.mask = fullfile(folder, 'mask.nii');
    files.noiseMask = fullfile(folder, 'noise_mask.nii');
    files.rtsnr = fullfile(folder, 'rtsnr.nii');
    files.median = fullfile(folder, 'median.nii');
    fiach_write_3d(mask, Vref, files.mask, overwrite, 'uint8');
    fiach_write_3d(noiseMask, Vref, files.noiseMask, overwrite, 'uint8');
    fiach_write_3d(rtsnrMap, Vref, files.rtsnr, overwrite, 'float32');
    fiach_write_3d(medianMap, Vref, files.median, overwrite, 'float32');
end


function fiach_write_3d(data, Vref, outputFile, overwrite, datatype)
    fiach_assert_writable(outputFile, overwrite);
    Vout = Vref;
    Vout.fname = outputFile;
    Vout.n = [1 1];
    if strcmpi(datatype, 'native')
        Vout.dt = Vref.dt;
    else
        Vout.dt = [spm_type(datatype), Vref.dt(2)];
    end
    Vout.pinfo = [1; 0; 0];
    Vout = spm_create_vol(Vout);
    spm_write_vol(Vout, data);
end


function fiach_write_4d(data, Vref, outputFile, overwrite)
    fiach_assert_writable(outputFile, overwrite);
    nTime = size(data, 4);
    Vout = repmat(Vref, nTime, 1);
    for timePoint = 1:nTime
        Vout(timePoint).fname = outputFile;
        Vout(timePoint).n = [timePoint 1];
        Vout(timePoint).dt = Vref.dt; % FIACH R preserves the input datatype.
        Vout(timePoint).pinfo = [1; 0; 0];
        Vout(timePoint) = spm_create_vol(Vout(timePoint));
        spm_write_vol(Vout(timePoint), data(:, :, :, timePoint));
    end
end


function fiach_assert_writable(outputFile, overwrite)
    if exist(outputFile, 'file') == 2 && ~overwrite
        error('fiach:OutputExists', ['Output already exists: %s\n' ...
            'Use ''Overwrite'', true to replace it.'], outputFile);
    end
end


function fiach_write_metrics(file, percentChanged, peakTsnr, overwrite)
    fiach_assert_writable(file, overwrite);
    fid = fopen(file, 'w');
    if fid < 0
        error('fiach:WriteFailed', 'Could not write %s.', file);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%% Data Changed\tPeak TSNR\n%.12g\t%.12g\n', percentChanged, peakTsnr);
end


function fiach_write_matrix(file, matrix, overwrite)
    fiach_assert_writable(file, overwrite);
    try
        writematrix(matrix, file, 'Delimiter', 'tab');
    catch writeException
        error('fiach:WriteFailed', 'Could not write %s: %s', file, writeException.message);
    end
end


function rp = fiach_read_motion_parameters(file, nTime)
    if exist(file, 'file') ~= 2
        error('fiach:RPNotFound', 'Realignment parameter file not found: %s', file);
    end
    rp = readmatrix(file, 'FileType', 'text');
    rp = rp(:, any(isfinite(rp), 1));
    if size(rp, 1) ~= nTime
        error('fiach:RPTimeMismatch', ['The number of functional volumes (%d) does not match ' ...
            'the number of rows in RP (%d).'], nTime, size(rp, 1));
    end
    if size(rp, 2) < 6
        error('fiach:InvalidRP', 'RP must contain at least 6 motion-parameter columns.');
    end
    rp = rp(:, 1:6);
end


function fd = fiach_framewise_displacement(rp)
% Direct port of R/fd.R: rotations converted with sin(theta) * 50 mm.
    displacement = rp(:, 1:3);
    angularDisplacement = sin(rp(:, 4:6)) * 50;
    totalDisplacement = [displacement, angularDisplacement];
    fd = [0; sqrt(sum(diff(totalDisplacement, 1, 1).^2, 2))];
end


function x = fiach_zero_nan(x)
    x(~isfinite(x)) = 0;
end


function fiach_write_mixture_plot(folder, rtsnr, model, threshold)
    figureHandle = figure('Visible', 'off');
    cleanupFigure = onCleanup(@() close(figureHandle));
    histogram(rtsnr, 100, 'Normalization', 'pdf');
    hold on;
    if ~model.fallback
        x = linspace(min(rtsnr), max(rtsnr), 500);
        y = zeros(size(x));
        for k = 1:numel(model.mu)
            component = model.lambda(k) ./ (model.sigma(k) * sqrt(2*pi)) .* ...
                exp(-.5 * ((x - model.mu(k)) ./ model.sigma(k)).^2);
            plot(x, component, 'LineWidth', 1.5);
            y = y + component;
        end
        plot(x, y, 'k--', 'LineWidth', 1.5);
    end
    if isfinite(threshold)
        xline(threshold, 'r-', 'LineWidth', 1.5);
    end
    xlabel('Robust tSNR'); ylabel('Density');
    title('FIACH noise segmentation');
    exportgraphics(figureHandle, fullfile(folder, 'segmentation.png'), 'Resolution', 150);
end
