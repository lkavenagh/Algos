%% Setup
[~, universe] = xlsread('PennyList.xlsx');
universe = strtrim(universe(:,1));

buffer = 1e3;
startingCash = 2e4;
maxExp = 1e3;
freq = 60;
maxTickers = 25;

% costPerTrade = 6.95;
costPerShare = 0.005;
costMin = 1;
costMaxPerc = 0.005;

%% Get returns
dates = getLastBusinessDates({'IBM'}, 1000, today-1, 0);

pClose = nan(size(universe,1), size(dates,2));

fprintf('Downloading returns...\n')

pClose = XFLDownload(universe, 'IQ_CLOSEPRICE_ADJ', dates);
spxClose = XFLDownload({'^SPX'}, 'IQ_CLOSEPRICE', dates);

tmp = [zeros(size(pClose,1), 1) pClose(:,1:end-1)];
ret = (pClose - tmp) ./ tmp;
ret(:,1) = 0;
ret(isnan(ret)) = 0;

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
cash = startingCash;
qty = zeros(size(ret));
exp = zeros(size(ret));
comm = zeros(size(ret));
numTrades = zeros(1, size(ret,2));
picks = [];

for i = freq:length(dates)
	% Pick out investable stocks
	if mod(i,freq) == 0
 		picks = find(pClose(:,i-1) < 0.05);
		r = round(randbetween(maxTickers, 1, length(picks)));
		picks = picks(r);
	else
		% Remove tickers that have risen enough
		qty(:,i) = qty(:,i-1);
		exp(:,i) = qty(:,i-1) .* pClose(:,i-1);
		
		toRemove = pClose(:,i-1) > 0.5;
		if sum(qty(toRemove,i))>0
			sum(qty(toRemove,i))
		end
		qty(toRemove,i) = 0;
		expRemove = sum(exp(toRemove,i));
		exp(toRemove,i) = 0;
		comm(:,i) = abs(qty(:,i)-qty(:,i-1)) .* costPerShare;
		comm(comm(:,i) > 0 & comm(:,i) < costMin, i) = costMin;

		cash = cash - sum(comm(:,i)) + expRemove;
		numTrades(i) = sum(qty(:,i) ~= qty(:,i-1));
		continue
	end
	
	% Add new tickers to portfolio
	
	yesterdaysExp = nansum(qty(:,i-1) .* pClose(:,i-1));
	cash = cash + yesterdaysExp;
	
	if cash < buffer
		fprintf(2, 'Ran out of money on %s\n', datestr(dates(i)));
		break
	end
	tmp = (cash - buffer) ./ length(picks);
	tmp = min(tmp, maxExp);
	exp(picks,i) = tmp;
	qty(:,i) = round(exp(:,i) ./ pClose(:,i-1));
	qty(isnan(qty)) = 0;
	
	exp(qty==0) = 0;
	
	if nansum(exp(:,i))==0 && ~isempty(picks)
		fprintf(2, 'Ran out of money on %s\n', datestr(dates(i)));
		break
	end
	
	todaysExp = nansum(qty(:,i) .* pClose(:,i-1));
	if isnan(todaysExp)
		todaysExp = 0;
	end
	cash = cash - todaysExp;
	
	numTrades(i) = sum(qty(:,i) ~= qty(:,i-1));
	
	tmp = qty(picks,i) .* costPerShare;
	tmp(tmp < costMin) = costMin;
	maxCost = costMaxPerc .* exp(picks,i);
	tmp(tmp > maxCost) = maxCost;
	comm(picks,i) = tmp;
	
	cash = cash - sum(comm(:,i));
	
% 	if todaysExp > buffer
% 		adjust = false;
% 		while cash < buffer
% 			% Rebalance down
% 			adjust = true;
% 			tmp = floor(qty(:,i) .* 0.95);
% 			exp1 = nansum(qty(:,i) .* pClose(:,i-1));
% 			exp2 = nansum(tmp .* pClose(:,i-1));
% 			cash = cash + (exp1-exp2);
% 			qty(:,i) = tmp;
% 			
% 			todaysExp = nansum(qty(:,i) .* pClose(:,i-1));
% 		end
% 		
% 		while cash > buffer && ~adjust
% 			% Rebalance up
% 			tmp = ceil(qty(:,i) .* 1.05);
% 			exp1 = nansum(qty(:,i) .* pClose(:,i-1));
% 			exp2 = nansum(tmp .* pClose(:,i-1));
% 			cash = cash - (exp2-exp1);
% 			qty(:,i) = tmp;
% 			
% 			todaysExp = sum(qty(:,i) .* pClose(:,i-1));
% 		end
% 	end
	
end
%% Calculate returns
fprintf('Backtest results: %s to %s\n', datestr(dates(1)), datestr(dates(end)));
fprintf('Starting cash: %s\n', util.Disp.AsDollars(startingCash));
fprintf('Max daily positions: %2.0f\n', maxTickers);

fprintf('\nAverage positive/negative universe return: %2.2f%%/%2.2f%%\n', 100*nanmean(ret(ret>0)), 100*nanmean(ret(ret<0)));
dollarRet = exp .* ret;

sharesTraded = diff(exp,1,2);

dailyRet = nansum(dollarRet,1);
dailyRet = dailyRet - sum(comm);

dailyExp = nansum(exp);
dailyExp(dailyExp==0) = nan;

y_qty = [nan(size(qty,1),1) qty(:,1:end-1)];
qtyDiff = qty-y_qty;
dailyTradeOut = sum(qtyDiff<0,1);
dailyTradeIn = sum(qtyDiff>0,1);
% fprintf('Return series:\n');
% for i = 1:length(dailyRet)
% 	if dailyExp(i) ~= 0;
% 		fprintf('%s: $%2.0f\n', datestr(dates(i)), dailyRet(i));
% 	end
% end

fprintf('\nTotal return: %s\n', util.Disp.AsDollars(sum(dailyRet)));

fprintf('Total exposure in longs: %s\n', util.Disp.AsDollars(nansum(exp(:,end-1))));
fprintf('Total LTD return: %2.2f%%\n', 100*nansum(dailyRet)/nanmean(dailyExp));

fprintf('\nTotal positive/negative days: %2.0f/%2.0f\n', sum(dailyRet>0), sum(dailyRet<0));
fprintf('Positive:negative ratio: %2.2f\n', sum(dailyRet>0)/sum(dailyRet<0));

fprintf('Average number of positions held: %2.1f\n', nanmean(nansum(qty>0,1)));
fprintf('Average number of daily trades: %2.1f\n', nanmean(dailyTradeOut+dailyTradeIn));

fprintf('Dollar return per day: %s\n', util.Disp.AsDollars(sum(dailyRet)/sum(dailyExp~=0)));

fprintf('\n')

% plot(dates,cumsum(dailyRet))
% datetick
% util.Plot.FormatPlot(gcf)
% 
% plot(dates,nansum(exp))
% datetick
% util.Plot.FormatPlot(gcf)