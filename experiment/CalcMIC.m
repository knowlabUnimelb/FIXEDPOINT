itemNumber = data(:,7); % specifies values in the 7th column of variable x

HH = mean(data(itemNumber==1, 13)); % the 13th column for all rows of x which have a value of 5 in the 7th column (i.e., the RTs for the Ex item)
HL = mean(data(itemNumber==2, 13));
LH = mean(data(itemNumber==3, 13));
LL = mean(data(itemNumber==4, 13));



MIC = (LL- LH) - (HL - HH)