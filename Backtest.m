%% Setup
[~, universe] = xlsread('S&PList.xlsx');
universe = strtrim(universe(:,1));

buffer = 1e3;
startingCash = 1e4;
freq = 60; %days

%costPerTrade = 6.95;
costPerShare = 0.005;
costMin = 1;

%% Get returns

dates = getLastBusinessDates(universe(1), 400, today-1, 0);

pClose = nan(size(universe,1), size(dates,2));

fprintf('Downloading returns...\n')

pClose = XFLDownload(universe, 'IQ_CLOSEPRICE_ADJ', dates);
spxClose = XFLDownload({'^SPX'}, 'IQ_CLOSEPRICE', dates);

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
cash = startingCash;
qty = zeros(size(ret));
exp = zeros(size(ret));
comm = zeros(size(ret));
numTrades = zeros(1, size(ret,2));
picks = [];

for i = freq:length(dates)
	% Pick out investable stocks
	if mod(i,freq) == 0
% 		picks = find(sum(alpha(:,i-freq+1:i-2),2) < 0 & sum(alpha(:,i-2:i-1),2) > 0, 10, 'first');
% 		picks = find(pClose(:,i-10) < pClose(:,i-(freq*2)+1) & pClose(i-1) > pClose(i-10));
% 		picks = find(all(alpha(:,i-7:i-1) > 0, 2));
		picks = find(all(alpha(:,i-6:i-2) < 0,2) & alpha(:,i-1) > 0);
% 		picks = find(all(beta(:,i-7:i-1) > 0.6 & beta(:,i-7:i-1) < 0.8, 2));
% 		picks = find(all(alpha(:,i-5:i-1)>0,2));
% 		picks = find(all((pClose(:,i-5:i-1)-pClose(:,i-6:i-2))>0,2));
	else
		% Remove badly performing tickers
		qty(:,i) = qty(:,i-1);
		exp(:,i) = qty(:,i-1) .* pClose(:,i-1);
		
% 		toRemove = alpha(:,i) < -0.01 & alpha(:,i-1) < -0.01;
		toRemove = false(size(alpha(:,i)));
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
	exp(picks,i) = (cash - buffer) ./ length(picks);
	qty(:,i) = round(exp(:,i) ./ pClose(:,i));
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
	
	comm(:,i) = abs(qty(:,i)-qty(:,i-1)) .* costPerShare;
	comm(comm(:,i) > 0 & comm(:,i) < costMin, i) = costMin;
	
	cash = cash - sum(comm(:,i));
	
	if todaysExp > buffer
		adjust = false;
		while cash < buffer
			% Rebalance down
			adjust = true;
			tmp = floor(qty(:,i) .* 0.95);
			exp1 = nansum(qty(:,i) .* pClose(:,i-1));
			exp2 = nansum(tmp .* pClose(:,i-1));
			cash = cash + (exp1-exp2);
			qty(:,i) = tmp;
			
			todaysExp = nansum(qty(:,i) .* pClose(:,i-1));
		end
		
		while cash > buffer && ~adjust
			% Rebalance up
			tmp = ceil(qty(:,i) .* 1.05);
			exp1 = nansum(qty(:,i) .* pClose(:,i-1));
			exp2 = nansum(tmp .* pClose(:,i-1));
			cash = cash - (exp2-exp1);
			qty(:,i) = tmp;
			
			todaysExp = sum(qty(:,i) .* pClose(:,i-1));
		end
	end
	
end
%% Calculate returns

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

fprintf('\nTotal return: $%2.0f\n', sum(dailyRet));

fprintf('Total exposure in longs: %s\n', util.Disp.AsDollars(nansum(exp(:,end-1))));
fprintf('Total LTD return: %2.2f%%\n', 100*nansum(dailyRet)/nanmean(dailyExp));

fprintf('\nTotal positive/negative days: %2.0f/%2.0f\n', sum(dailyRet>0), sum(dailyRet<0));
fprintf('Positive:negative ratio: %2.2f\n', sum(dailyRet>0)/sum(dailyRet<0));

fprintf('Average number of positions held: %2.1f\n', nanmean(nansum(qty>0,1)));
fprintf('Average number of daily trades: %2.1f\n', nanmean(dailyTradeOut+dailyTradeIn));

fprintf('Dollar return per day: %s\n', util.Disp.AsDollars(sum(dailyRet)/sum(dailyExp~=0)));

fprintf('\n')

plot(dates,cumsum(dailyRet))
datetick
util.Plot.FormatPlot(gcf)