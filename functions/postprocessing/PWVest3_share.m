% function qDiffSqr = PWVest3_share(inParams,D,F,tres,Q)% Estimates an "underlying" velocity waveform and a pwvPWVest2(inParams,D,F,tres,Q)
% n = number of cross-sections
% m = number of timepoints in one velocity waveform
% Input:
% inparams: Initial guesses for the velocity waveform and PWV, waveform
%           datapoints first and PWV last. 
%           Size [1 x m+1]
% distance: Vector containing the distanses from the seed-point
%           cross-sections of the arterial tree to the subsequent
%           cross-sections. Usually sorted after distance. In meters. 
%           Size [n x 1]
% waveMat:  Matrix containing the velocity waveform for each
%           cross-section. Same sorting as in "distance".
%           Usually normalized to have zero mean and unit std.
%           Size [n x m]
% tRes:     Time between consecutive frames in seconds.
% scaling:  Weight factor for each cross-section.
%           "cross-sections-area / scalingFactor^2", where scalingFactor is
%           what was used to variance-normalize waveforms.
%           scaling can probably be set to ones(n,1) for a single vessel of
%           approximatly equal size along its path.
%           Size [n x 1]
% Example Usage:
% fun1=@(inParams)PWVest3_share(inParams,d,waveMat,tRes,w); 
% pwv0 = 10; %initial guess of pwv
% mean_flow = mean(waveMat); %initial guess of waveform
% initialGuess=[mean_flow, pwv0]; 
% options = optimset('Display','iter', 'TolCon', 1e-7, 'TolX', 1e-7, 'TolFun', 1e-7,'DiffMinChange', 1e-3);
% [params,exitflag,output] = fminunc(fun1,initialGuess, options);%,options1);
% pwv = params(end) % have a look at the PWV
% FOR REFERENCE SEE BJÖRNFOT et al JCBFM 2021
function qDiffSqr = PWVest3_share(inParams,D,F,tres,Q)
    m = size(inParams,2)-1; %number of samples
    tV = 0:tres:(tres*3*m-tres); %Triples the time series, (basically to avoid overlap when sampling the shift)
    velocity = inParams(1:m); %velocity guess
    pwv = inParams(m+1); % current pwv
    region = [m+1:m*2]; % region, is 21:40 (middle of a 3x time series
    deltaT = repmat(tV,length(D),1)-repmat(D,1,length(tV))/pwv; %First term, repeats time array for each row, secong term estimates a time shift for each D (row) dt=D/pwv
    vShift=interp1(tV,repmat(velocity,1,3),deltaT,'linear'); %triple time, triple flow (velocity??), 
    weightm = repmat(Q,1,m); %Weight matrix for each point
    vShift = vShift(:,region); %collect only the middle waveform
    qDiffSqr=sum(weightm(:).*(vShift(:)-F(:)).^2);
end
