function designOptions = dg_designOptions()
% deduces from the information what design the user wants to do 

% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% design options
% design code is in binary:
% 1st column tells if there are repeated-measure factors
% 2nd column tells if there are independent groups

designStruct = struct();
designStruct(1).codeBin = num2str([0 0]);
designStruct(1).lbl = 'continuous association';

designStruct(2).codeBin = num2str([0 1]);
designStruct(2).lbl = 'group comparison(s)';

designStruct(3).codeBin = num2str([1 0]);
designStruct(3).lbl = 'condition comparison(s)';

designStruct(4).codeBin = num2str([1 1]);
designStruct(4).lbl = 'group comparison(s) and condition comparison(s)';

designOptions = struct2table(designStruct);

% transform code from binary to decimal
designOptions.codeDec = nan(size(designOptions,1),1);
for rIdx = 1:size(designOptions,1)
    designOptions{rIdx,"codeDec"} = bin2dec(char(designOptions.codeBin(rIdx,:)));
end
designOptions = sortrows(designOptions(:,[3 1 2]),"codeDec");
