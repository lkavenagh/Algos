%% Setup
[~, universe] = xlsread('S&PList.xlsx');
universe = strtrim(universe(:,1));

buffer = 1e3;
startingCash = 1e4;

%costPerTrade = 6.95;
costPerShare = 0.005;
costMin = 1;

%% Get returns
dates = today-100:today;
dates = dates(isbusday(dates));

pClose = nan(size(universe,1), size(dates,2));

fprintf('Downloading returns...\n')

conn = yahoo;
for i = 1:length(universe)
	disp(length(universe)-i+1)
	try
		tmp = fetch(conn, universe{i}, 'Close', datestr(dates(1)), datestr(dates(end)));
	catch
		tmp = [nan nan];
	end
	dateIdx = ismember(dates,tmp(:,1));
	pClose(i,dateIdx) = tmp(:,2);
end

tmp = fetch(conn, universe{i}, 'Close', datestr(dates(1)), datestr(dates(end)));
dateIdx = ismember(dates,tmp(:,1));
spxClose = nan(1,length(dates));
spxClose(1,dateIdx) = tmp(:,2);

tmp = [zeros(size(pClose,1), 1) pClose(:,1:end-1)];
ret = (pClose - tmp) ./ tmp;
ret(:,1) = 0;

tmp = [0 spxClose(1:end-1)];
spxRet = (spxClose - tmp) ./ tmp;
spxRet(1) = 0;

beta = ones(size(ret));
for i = 1:size(ret,1)
	disp(size(ret,1)-i+1)
	for j = 61:length(dates)
		tmp = corrcoef(spxRet(j-60:j-1)', ret(i,j-60:j-1)');
		beta(i,j) = tmp(2);
	end
end

alpha = ret - beta.*repmat(spxRet, size(universe,1), 1);

%% Run simulation

% Pick out investable stocks
% Run at beginning of quarter, or if no stocks were entered in beginning of
% quarter.

% 		picks = find(sum(alpha(:,i-freq+1:i-2),2) < 0 & sum(alpha(:,i-2:i-1),2) > 0, 10, 'first');
% 		picks = find(pClose(:,i-10) < pClose(:,i-(freq*2)+1) & pClose(i-1) > pClose(i-10));
% 		picks = find(all(alpha(:,i-7:i-1) > 0, 2));
picks = find(all(alpha(:,end-6:end-2) < 0,2) & alpha(:,end-1) > 0);
% 		picks = find(all(beta(:,i-7:i-1) > 0.6 & beta(:,i-7:i-1) < 0.8, 2));
% 		picks = find(all(alpha(:,i-5:i-1)>0,2));
% 		picks = find(all((pClose(:,i-5:i-1)-pClose(:,i-6:i-2))>0,2));

if isempty(picks)
	fprintf('No stocks were selected today.\n');
	return
end

% Print out tickers
fprintf('The following tickers should be entered today:\n')
for i = 1:length(picks)
	fprintf('\t%s\n', universe{i});
end

idx = qty(:,end-1) ~= 0;
fprintf('The following tickers should be exited today:\n');