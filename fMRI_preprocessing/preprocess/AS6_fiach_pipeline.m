%% AS6_fiach_pipeline.m
% Headless FIACH step for wrapper_function_preprocess.sh (called by
% AS6_call_FIACH.sh, which defines `cfg` and then runs this script).
%
% Runs the packaged fiach.m unchanged, then:
%   - writes noisyPCs_<input>.tsv (the 6 FIACH PCs, no motion) for AS7
%   - writes QC maps/plots of which samples FIACH changed
%   - task mode: rclean_<input>.nii = FIACH output, nothing regressed
%                (motion + PCs go into the task GLM via AS7)
%   - rs mode:   rclean_<input>.nii = FIACH output with motion + PCs
%                projected out voxelwise (temporal mean retained)
%
% Required cfg fields:
%   fiachRoot spmRoot fiachDir input mask motionPar TR TEms B0 mode
%   gmmMethod posteriorThreshold contrastProfile nMads freq maxGap
%   rsIncludeMotion

if ~exist('cfg', 'var') || ~isstruct(cfg)
    error('AS6_fiach_pipeline:NoCfg', 'cfg must be defined before running this script.');
end
cfg.mode = lower(char(cfg.mode));
if ~ismember(cfg.mode, {'task', 'rs'})
    error('AS6_fiach_pipeline:BadMode', 'cfg.mode must be ''task'' or ''rs'' (got ''%s'').', cfg.mode);
end


addpath(cfg.fiachRoot);          % preprocess/ — holds fiach.m and fiach_bold_contrast.m
if ~isempty(cfg.spmRoot)
    addpath(cfg.spmRoot);        % only when SPM_DIR was set
end
if exist('spm', 'file') ~= 2
    error('AS6_fiach_pipeline:NoSPM', ...
        'SPM is not on the MATLAB path. Add it once in MATLAB (addpath + savepath) or run with SPM_DIR=/path/to/spm12.');
end
spm('defaults', 'fmri');
fprintf('SPM used    : %s\n', which('spm'));
fprintf('fiach.m used: %s\n', which('fiach'));


[~, inputName, ~] = fileparts(cfg.input);
filtFile    = fullfile(cfg.fiachDir, ['filt_' inputName '.nii']);
cleanFile   = fullfile(cfg.fiachDir, ['rclean_' inputName '.nii']);
pcsFile     = fullfile(cfg.fiachDir, ['noisyPCs_' inputName '.tsv']);
rpFile      = fullfile(cfg.fiachDir, ['rp_' inputName '.txt']);
summaryFile = fullfile(cfg.fiachDir, 'FIACH_summary.txt');

fprintf('\n=== FIACH (%s mode) ===\n', cfg.mode);
fprintf('Input : %s\n', cfg.input);
if isempty(cfg.mask)
    fprintf('Mask  : FIACH automatic\n');
else
    fprintf('Mask  : %s\n', cfg.mask);
end
fprintf('TR %.4g s | TE %d ms | B0 %g T | GMM %s\n', cfg.TR, cfg.TEms, cfg.B0, cfg.gmmMethod);


%% 1) Motion parameters: mcflirt .par -> SPM/FIACH rp order -----------------
% mcflirt .par columns: rx ry rz (rad) tx ty tz (mm)
% FIACH expects SPM order: tx ty tz (mm) rx ry rz (rad) - used for FD and noise_basis6
par = readmatrix(cfg.motionPar, 'FileType', 'text');
par = par(:, all(isfinite(par), 1));   % drop empty trailing column from mcflirt spacing
if size(par, 2) ~= 6
    error('AS6_fiach_pipeline:BadPar', 'Expected 6 columns in %s, found %d.', ...
        cfg.motionPar, size(par, 2));
end
rp = par(:, [4 5 6 1 2 3]);
writematrix(rp, rpFile, 'FileType', 'text', 'Delimiter', 'tab');


%% 2) FIACH (packaged implementation, unchanged) ----------------------------
t = fiach_bold_contrast(cfg.B0, cfg.TEms, 'Profile', cfg.contrastProfile);
fprintf('Deterministic threshold T = %.4f %% (%s profile)\n', t, cfg.contrastProfile);

args = {'RP', rpFile, ...
    'MaxGap', cfg.maxGap, ...
    'Freq', cfg.freq, ...
    'NMads', cfg.nMads, ...
    'GMMMethod', cfg.gmmMethod, ...
    'GMMPosteriorThreshold', cfg.posteriorThreshold, ...
    'Overwrite', true, ...          % reruns are controlled by the wrapper checkpoints
    'MakePlots', true};             % -> <diagnostics>/segmentation.png
if ~isempty(cfg.mask)
    args = [args, {'Mask', cfg.mask}];
end

result = fiach(cfg.input, t, cfg.TR, args{:});
diagDir = result.diagnosticDirectory;
fprintf('Mask source: %s | GMM: %s\n', result.maskSource, result.gmmMethod);

% PCs only (no motion) in the same file name/format the old run_fiach wrote
writematrix(result.noiseRegressors, pcsFile, 'FileType', 'text', 'Delimiter', 'tab');
fprintf('Wrote %s\n', pcsFile);


%% 3) Which samples did FIACH change? (QC) ----------------------------------
% fiach.m does not export its bad-data mask, so the flags are recomputed
% here with the same definitions (raw median, DCT high-pass from
% R/highPass.R, scaled MAD after high-pass, median +/- (T% + NMads*MAD)).
% The total is checked against FIACH's own nInitiallyFlagged below.
Vmask = spm_vol(fullfile(diagDir, 'mask.nii'));
mask = spm_read_vols(Vmask) > 0;
spatialSize = Vmask.dim(1:3);
maskIdx = find(mask);
nBrain = numel(maskIdx);

V = spm_vol(cfg.input);
nT = numel(V);
brain = zeros(nT, nBrain);
for frame = 1:nT
    vol = spm_read_vols(V(frame));
    vol(~isfinite(vol)) = 0;
    brain(frame, :) = vol(maskIdx);
end
brainMedians = median(brain, 1);

% DCT high-pass basis count, rounding .5 to even as R does
nBasisRaw = 2 * (nT * cfg.TR / (cfg.freq + 1));
if abs(nBasisRaw - floor(nBasisRaw) - .5) <= 4 * eps(max(1, nBasisRaw))
    nBasis = floor(nBasisRaw) + mod(floor(nBasisRaw), 2);
else
    nBasis = floor(nBasisRaw + .5);
end
hp = brain;
if nBasis >= 2
    n = (0:nT-1).';
    basis = zeros(nT, nBasis - 1);
    for b = 2:nBasis
        basis(:, b - 1) = sqrt(2 / nT) * cos(pi * (2 * n + 1) * (b - 1) / (2 * nT));
    end
    hp = brain - basis * (basis.' * brain);
end
clear brain;

brainMads = median(abs(hp - median(hp, 1)), 1) * 1.4826;
upperBound = brainMedians + brainMedians / 100 * t + cfg.nMads * brainMads;
lowerBound = brainMedians - brainMedians / 100 * t - cfg.nMads * brainMads;
flagged = hp < lowerBound | hp > upperBound;          % nT x nBrain
clear upperBound lowerBound;

nFlagged = nnz(flagged);
if nFlagged ~= result.metrics.nInitiallyFlagged
    warning('AS6_fiach_pipeline:FlagMismatch', ...
        'QC flag count (%d) differs from FIACH''s nInitiallyFlagged (%d).', ...
        nFlagged, result.metrics.nInitiallyFlagged);
end
flagCount = sum(flagged, 1);
[~, worst] = maxk(flagCount, 3);
hpWorst = hp(:, worst);
clear hp;

% Per-voxel maps (same grid/affine as the FIACH mask)
Vout = Vmask;
if isfield(Vout, 'private'), Vout = rmfield(Vout, 'private'); end
Vout.dt = [spm_type('float32') spm_platform('bigend')];
Vout.pinfo = [1; 0; 0];
Vout.n = [1 1];

countMap = zeros(spatialSize);
countMap(maskIdx) = flagCount;
Vout.fname = fullfile(cfg.fiachDir, 'FIACH_changed_count.nii');
Vout.descrip = 'FIACH: number of changed samples per voxel';
spm_write_vol(Vout, countMap);

pctMap = 100 * countMap / nT;
Vout.fname = fullfile(cfg.fiachDir, 'FIACH_changed_pct.nii');
Vout.descrip = 'FIACH: % of frames changed per voxel';
spm_write_vol(Vout, pctMap);


%% 4) FIACH output, in-mask -------------------------------------------------
Vfilt = spm_vol(filtFile);
brainFilt = zeros(nT, nBrain, 'single');
for frame = 1:nT
    vol = spm_read_vols(Vfilt(frame));
    brainFilt(frame, :) = single(vol(maskIdx));
end


%% 5) Mode-specific final output -> rclean_<input>.nii ----------------------
nuisanceRank = NaN;
residualDf = NaN;
cleanedWorst = [];
if strcmp(cfg.mode, 'task')
    % Task: keep regressors for the GLM; rclean is the FIACH output as-is
    copyfile(filtFile, cleanFile, 'f');
    fprintf('Task mode: %s = FIACH output (no regression)\n', cleanFile);
else
    % rs: project nuisance subspace out of each in-mask voxel
    designRaw = result.noiseRegressors;
    if cfg.rsIncludeMotion
        designRaw = [rp, designRaw];
    end
    centred = designRaw - mean(designRaw, 1);
    columnNorm = sqrt(sum(centred .^ 2, 1));
    usable = columnNorm > max(columnNorm) * eps * nT;
    scaled = centred(:, usable) ./ columnNorm(usable);
    [U, S, ~] = svd(scaled, 'econ');
    spectrum = diag(S);
    nuisanceRank = nnz(spectrum > max(size(scaled)) * eps(max(spectrum)));
    Q = U(:, 1:nuisanceRank);
    residualDf = nT - nuisanceRank - 1;
    writematrix(designRaw, fullfile(cfg.fiachDir, 'FIACH_rs_design.txt'), ...
        'FileType', 'text', 'Delimiter', 'tab');

    cleaned = zeros(size(brainFilt), 'single');
    chunkSize = 10000;
    for first = 1:chunkSize:nBrain
        last = min(first + chunkSize - 1, nBrain);
        values = double(brainFilt(:, first:last));
        cleaned(:, first:last) = single(values - Q * (Q.' * values));
    end
    cleanedWorst = double(cleaned(:, worst));

    % 4-D float32, out-of-mask voxels left as FIACH wrote them
    if exist(cleanFile, 'file'), delete(cleanFile); end
    Vtmpl = Vfilt(1);
    if isfield(Vtmpl, 'private'), Vtmpl = rmfield(Vtmpl, 'private'); end
    Vtmpl.fname = cleanFile;
    Vtmpl.dt = [spm_type('float32') spm_platform('bigend')];
    Vtmpl.pinfo = [1; 0; 0];
    Vtmpl.descrip = 'FIACH + motion/PC nuisance projection (rs)';
    for frame = 1:nT
        Vo = Vtmpl;
        Vo.n = [frame 1];
        Vo = spm_create_vol(Vo);
        vol = spm_read_vols(Vfilt(frame));
        vol(maskIdx) = double(cleaned(frame, :));
        spm_write_vol(Vo, vol);
    end
    clear cleaned;
    fprintf('rs mode: regressed %d nuisance dimensions (residual df %d) -> %s\n', ...
        nuisanceRank, residualDf, cleanFile);
end


%% 6) QC figures -------------------------------------------------------------
% % of in-mask voxels changed per frame
fig = figure('Color', 'w', 'Position', [100 100 900 350]);
plot(100 * mean(flagged, 2), 'LineWidth', 1.1, 'Color', [0.1 0.35 0.75]);
grid on;
xlabel('Frame');
ylabel('% in-mask voxels changed');
title('FIACH corrections per frame');
saveas(fig, fullfile(cfg.fiachDir, 'FIACH_QC_corrections_per_frame.png'));
close(fig);

% Timecourses of the three most-corrected voxels
for k = 1:numel(worst)
    v = worst(k);
    [vx, vy, vz] = ind2sub(spatialSize, maskIdx(v));
    bad = find(flagged(:, v));
    fig = figure('Color', 'w', 'Position', [100 100 900 350]);
    plot(hpWorst(:, k), '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
    hold on;
    plot(double(brainFilt(:, v)), '-', 'Color', [0.85 0.25 0.15], 'LineWidth', 1.1);
    legendText = {'before (high-passed input)', 'after FIACH'};
    if ~isempty(cleanedWorst)
        plot(cleanedWorst(:, k), '-', 'Color', [0.15 0.55 0.25], 'LineWidth', 1.1);
        legendText{end+1} = 'after nuisance regression (rs)'; %#ok<SAGROW>
    end
    plot(bad, hpWorst(bad, k), 'kx', 'MarkerSize', 6);
    legendText{end+1} = 'changed samples'; %#ok<SAGROW>
    legend(legendText, 'Location', 'best');
    grid on;
    xlabel('Frame');
    ylabel('Signal intensity');
    title(sprintf('Voxel [%d %d %d]: %d / %d samples changed', vx, vy, vz, flagCount(v), nT));
    saveas(fig, fullfile(cfg.fiachDir, sprintf('FIACH_QC_worstvoxel_%d.png', k)));
    close(fig);
end


%% 7) Summary for the QC report ---------------------------------------------
noiseMask = spm_read_vols(spm_vol(fullfile(diagDir, 'noise_mask.nii'))) > 0;
fid = fopen(summaryFile, 'w');
fprintf(fid, 'Mode\t%s\n', cfg.mode);
fprintf(fid, 'GMM method\t%s\n', result.gmmMethod);
fprintf(fid, 'Mask source\t%s\n', result.maskSource);
fprintf(fid, 'B0 (T)\t%g\n', cfg.B0);
fprintf(fid, 'TR (s)\t%g\n', cfg.TR);
fprintf(fid, 'TE (ms)\t%d\n', cfg.TEms);
fprintf(fid, 'Threshold T (%%)\t%.4f\n', t);
fprintf(fid, 'NMads\t%g\n', cfg.nMads);
fprintf(fid, 'High-pass period (s)\t%g\n', cfg.freq);
fprintf(fid, 'MaxGap (volumes)\t%d\n', cfg.maxGap);
fprintf(fid, 'Volumes\t%d\n', nT);
fprintf(fid, 'Brain-mask voxels\t%d\n', nBrain);
fprintf(fid, 'Noise-mask voxels\t%d\n', nnz(noiseMask));
fprintf(fid, 'Peak tSNR\t%.4f\n', result.metrics.peakTSNR);
fprintf(fid, 'FIACH %% data changed\t%.4f\n', result.metrics.percentDataChanged);
fprintf(fid, 'Samples changed\t%d\n', nFlagged);
fprintf(fid, 'Samples changed (%% of in-mask)\t%.4f\n', 100 * nFlagged / numel(flagged));
fprintf(fid, 'Voxels with >=1 change\t%d\n', nnz(flagCount));
fprintf(fid, 'Frames with >=1 change\t%d\n', nnz(any(flagged, 2)));
fprintf(fid, 'Max changes in one voxel\t%d\n', max(flagCount));
if strcmp(cfg.mode, 'rs')
    if cfg.rsIncludeMotion
        fprintf(fid, 'rs regressors\t6 motion + %d FIACH PCs\n', size(result.noiseRegressors, 2));
    else
        fprintf(fid, 'rs regressors\t%d FIACH PCs\n', size(result.noiseRegressors, 2));
    end
    fprintf(fid, 'Nuisance rank\t%d\n', nuisanceRank);
    fprintf(fid, 'Residual df\t%d\n', residualDf);
end
fclose(fid);
fprintf('Wrote %s\n=== FIACH done ===\n', summaryFile);