% https://fr.mathworks.com/help/econ/simulate-responses-of-an-estimated-varx-model.html

load Data_USEconModel
varNaN = any(ismissing(DataTable),1); % Variables containing NaN values
seriesWithNaNs = series(varNaN)

idx = all(~ismissing(DataTable(:,{'UNRATE' 'M1SL'})),2);

rGDP = DataTable.GDP(idx)./(DataTable.GDPDEF(idx)/100);
rM1SL = DataTable.M1SL(idx)./(DataTable.GDPDEF(idx)/100);

dLRGDP = diff(log(rGDP));              % rGDP growth rate
dLRM1SL = diff(log(rM1SL));            % rM1SL growth rate
d3MTB = diff(DataTable.TB3MS(idx));    % Change in short-term interest rate (3MTB)
dUNRATE = diff(DataTable.UNRATE(idx)); % Change in unemployment rate

T = numel(d3MTB);    % Total sample size
oosT = 12;           % Out-of-sample size
estT = T - oosT;     % Estimation sample size
estIdx = 1:estT;     % Estimation sample indices
oosIdx = (T - 11):T; % Out-of-sample indices
dates = dates((end - T + 1):end);

EstY = [dLRGDP(estIdx) dLRM1SL(estIdx) d3MTB(estIdx)]; % In-sample responses
estX = dUNRATE(estIdx);                                % In-sample exogenous data
n = size(EstY,2);

OOSY = [dLRGDP(oosIdx) dLRM1SL(oosIdx) d3MTB(oosIdx)]; % Out-of-sample responses
oosX = dUNRATE(oosIdx);                                % Out-of-sample exogenous data


Mdl = varm(n,4);


EstMdl = estimate(Mdl,EstY,'X',estX);
summarize(EstMdl)

numPaths = 1000;
Y0 = EstY((end-3):end,:);
rng(1); % For reproducibility
YSim = simulate(EstMdl,oosT,'X',oosX,'Y0',Y0,'NumPaths',numPaths);


YSimBar = mean(YSim,3);
YSimQrtl = quantile(YSim,[0.05 0.25 0.5 0.75 0.95],3);
RepDates = repmat(dates(oosIdx),1,1000);
respNames = {'dLRGDP' 'dLRM1SL' 'd3MTB'};

figure;
for j = 1:n
    subplot(3,1,j);
    h1 = plot(dates(oosIdx),squeeze(YSim(:,j,:)),'Color',0.75*ones(3,1));
    hold on;   
    h2 = plot(dates(oosIdx),YSimBar(:,j),'.-k','LineWidth',2);
    h3 = plot(dates(oosIdx),squeeze(YSimQrtl(:,j,:)),':r','LineWidth',1.5);
    h4 = plot(dates((end - 30):end),[EstY((end - 18):end,j);OOSY(:,j)],...
        'b','LineWidth',2);
    title(sprintf('%s',respNames{j}));
    datetick;
    axis tight;
    hold off;
end
legend([h1(1) h2(1) h3(1) h4],{'Simulated Series','Simulation Mean',...
    'Simulation Quartile','Data'},'Location',[0.4 0.1 0.01 0.01],...
    'FontSize',8);

MdlUNRATE = arima('ARLags',1:4);
EstMdlUNRATE = estimate(MdlUNRATE,estX,'Display','off');

XSim = simulate(EstMdlUNRATE,oosT,'Y0',estX(end-3:end),...
    'NumPaths',numPaths);

YSimRX = zeros(oosT,n,numPaths); % Preallocate

for j = 1:numPaths
    YSimRX(:,:,j) = simulate(EstMdl,oosT,'X',XSim(:,j),'Y0',Y0);
end

YSimBarRX = mean(YSimRX,3);
YSimQrtlRX = quantile(YSimRX,[0.05 0.25 0.5 0.75 0.95],3);

figure;
for j = 1:n;
    subplot(3,1,j);
    h1 = plot(dates(oosIdx),squeeze(YSimRX(:,j,:)),'Color',0.75*ones(3,1));
    hold on;   
    h2 = plot(dates(oosIdx),YSimBarRX(:,j),'.-k','LineWidth',2);
    h3 = plot(dates(oosIdx),squeeze(YSimQrtlRX(:,j,:)),':r','LineWidth',1.5);
    h4 = plot(dates((end - 30):end),[EstY((end - 18):end,j);OOSY(:,j)],...
        'b','LineWidth',2);
    title(sprintf('%s with Simulated Unemployment Rate',respNames{j}));
    datetick;
    axis tight;
    hold off;
end
legend([h1(1) h2(1) h3(1) h4],{'Simulated Series','Simulation Mean',...
    'Simulation Quartile','Data'},'Location',[0.4 0.1 0.01 0.01],...
    'FontSize',8)