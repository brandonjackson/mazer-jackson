function [] = preel(pf, lag, winsize, varargin)
% PREEL plots a summary of the reel response for the given PF file
% 
% The work is done by Touryan's load_reel and reel_summary files. Make sure
% that these are in available in your path.
%
%    [] = preel(pf, latency, winsize, [options])
%
%  INPUT
%    file      - p2m PF struct
%    latency   - temporal latency (ms)
%    winsize   - temporal integration window (ms)
%
%  OUTPUT
%    (figure)


Reel_Info = load_reel(pf,lag,winsize);
reel_summary(Reel_Info);

end

