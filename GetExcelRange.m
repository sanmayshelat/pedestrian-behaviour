% From: https://nl.mathworks.com/matlabcentral/answers/179438-how-to-read-dynamic-range-with-xlsread-function

function xlsrange = GetExcelRange(startrow, endrow, startcol, endcol)
   xlsrange = sprintf('%s%d:%s%d', GetExcelColumn(startcol), startrow, GetExcelColumn(endcol), endrow);
end