function pf = pffind( query )
% PFFIND Loads p2m files and performs basic processing

% strip down full path if present
if strncmp(query,'/auto/',5)==1
    [~,exper_name,exper_no] = fileparts(query);
    query = [exper_name exper_no];
end

pf = dbfind(query);

if PFUtil.spikeCount(pf) < 100
    error('Insufficient Spikes Found in P2M File');
end

if strcmp(PFUtil.taskName(pf),'curvplay')
    pf = PFUtil.copyParam(pf,'IMAGE_INFO');
end
pf = PFUtil.removeBadTrials(pf);

end