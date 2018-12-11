function [ MeanParams, paramarray ] = CalibAvg( CalibCell, TrueParams )
%CALIBAVG Gives Statistics on the Results

nCal = numel(CalibCell);
nState = numel(TrueParams.nu);

paramarray = zeros([nCal,3*nState]);
nuArray = zeros([nCal,nState]);
QArray = zeros([nCal,nState,nState]);

Tparamarray = [TrueParams.mu,TrueParams.kappa,TrueParams.ThetaValues];

for k=1:numel(CalibCell)
    tempParam = CalibCell{k};
    paramarray(k,:) = [tempParam.mu,tempParam.kappa,tempParam.ThetaValues];
    nuArray(k,:) = tempParam.nu;
    QArray(k,:,:) = tempParam.Q;
end

paramMDist = mean(paramarray,1);

MeanParams.mu = paramMDist(1:nState);
MeanParams.kappa = paramMDist((nState+1):(2*nState));
MeanParams.ThetaValues = paramMDist((2*nState+1):(3*nState));
MeanParams.Q = mean(QArray,1);
MeanParams.nu = mean(nuArray,1);

paramMDist = paramMDist - Tparamarray;
paramDiff = bsxfun(@minus,paramarray,Tparamarray);
paramVDist = var(paramDiff,1);

end

