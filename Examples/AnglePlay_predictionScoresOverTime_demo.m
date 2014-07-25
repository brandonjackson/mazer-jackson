% Prediction Scores Over Time
% Compares how well the AnglePlay kernel explains variance versus the
% maximum amount of explainable variance for a variety of offsets and
% winsizes.

pf = pffind('romeo0295*curvplay');
AnglePlay.predictionScoresOverTime(pf);
