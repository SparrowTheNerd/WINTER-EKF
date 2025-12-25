%  Processes CSV data input from InertialExport,
%  an OpenRocket simulation extension

function ORDat = ORDataImport(dT)
ORraw = readtable("OpenRocket Data/simulation-001.csv");

TT = table2timetable(ORraw, 'RowTimes', seconds(ORraw.t));
% retime resamples time data and resolves irregular timing, like from OR
TT_100Hz = retime(TT, seconds(0:dT:max(ORraw.t)), 'linear');

ORDat = timetable2table(TT_100Hz, 'ConvertRowTimes', true);
ORDat = removevars(ORDat,"Time");
clear TT TT_100Hz ORraw
end
