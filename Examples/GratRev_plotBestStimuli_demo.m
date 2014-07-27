% Plot Best Stimuli
% Sorts stimuli based on their evoked response, then plots the kernel,
% details of the response distribution, and a grid of the best stimuli.
GR = GratRev(pffind('romeo0295*gratrev*003'),50,50);
GR.plotBestStimuli();