function [] = preel(pf, lag, winsize, varargin)
% PREEL plots a summary of the reel response for the given PF file

Reel_Info = load_reel(pf,lag,winsize);
reel_summary(Reel_Info);

end

