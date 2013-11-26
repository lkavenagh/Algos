for k = 1:100
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
% 				sum(qty(toRemove,i))
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
	
	dollarRet = exp .* ret;
	dailyRet = nansum(dollarRet,1);
	dailyRet = dailyRet - sum(comm);
	montecarlo(k,:) = dailyRet;
end

fprintf('Backtest results: %s to %s\n', datestr(dates(1)), datestr(dates(end)));
fprintf('Number of iterations: %2.0f\n', size(montecarlo,1))
fprintf('Starting cash: %s\n', util.Disp.AsDollars(startingCash));
fprintf('Max daily positions: %2.0f\n', maxTickers);
fprintf('Rebalance every %2.0f days (remove >50%% moves every day)\n', freq);

fprintf('\nAverage positive/negative universe return: %2.2f%%/%2.2f%%\n', 100*nanmean(ret(ret>0)), 100*nanmean(ret(ret<0)));

fprintf('\nAverage total return over whole period: %s\n', util.Disp.AsDollars(nanmean(nansum(montecarlo,2))));

fprintf('Average dollar return per day: %s\n', util.Disp.AsDollars(nanmean(nanmean(montecarlo,1))));

fprintf('Number of simulations with negative returns: %2.0f out of %2.0f\n', sum(sum(montecarlo,2)<0), size(montecarlo,1));

fprintf('\n')


plot(cumsum(montecarlo,2)')
util.Plot.FormatPlot(gcf)