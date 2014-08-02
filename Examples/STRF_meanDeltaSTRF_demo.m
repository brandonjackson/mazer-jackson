% Find mean STRF and clusters from batch pstrf_delta Results

% Step 1: Specify the files to perform batch analysis on (as a cell
% array whose contents will be fed into dbfind to get the actual file name)
files = cell(10,1);
files{1} = 'bert0255.curvplay.002';
files{2} = 'bert0265.curvplay.005';
files{3} = 'bert0266.curvplay.004';
files{4} = 'bert0270.curvplay.007';
files{5} = 'romeo0284.curvplay.004';
files{6} = 'romeo0284.curvplay.005';
files{7} = 'romeo0291.curvplay.002';
files{8} = 'romeo0291.curvplay.003';
files{9} = 'romeo0295.curvplay.004';
files{10} = 'romeo0300.curvplay.005';

% Step 2: Perform Batch of pstrf_delta Analyses, Saving Results to Directory
batch_filedir = 'Examples/batch_pstrf_delta_results/';
mkdir(batch_filedir);
STRF.batch_pstrf_delta(files,'auto',batch_filedir);

% Step 3: Analyze Batch Results
STRF.meanDeltaSTRF(batch_filedir,'auto',2);
