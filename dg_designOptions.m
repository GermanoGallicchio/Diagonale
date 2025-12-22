function designOptions = dg_designOptions()
% Returns table of design codes and their interpretations
%
% Creates a lookup table mapping 2-digit binary design codes. 
% Used for reporting and validation.
%
% OUTPUT:
%   designOptions - table with columns:
%                   .codeBin: string representation of binary code ('0 0', '0 1', etc)
%                   .codeDec: decimal representation of binary code
%                   .lbl: description of the design
%
%                   Binary codes represent:
%                   1. [0 0] - Correlation/association analysis
%                   2. [0 1] - Independent observation comparison
%                   3. [1 0] - Repeated measures comparison
%                   4. [1 1] - Mixed design (independent + repeated measures)
%
%                   In other words, each binary digit tells about the design:
%   Digit 1: 0 = no repeated measures; 1 = has repeated measures
%   Digit 2: 0 = only one group of observations; 1 = has independent groups of observations, 
%
% Author: Germano Gallicchio (germano.gallicchio@gmail.com)

%% design options
% design code is in binary:
% 1st column tells if there are repeated-measure factors
% 2nd column tells if there are independent groups of observations

designStruct = struct();
designStruct(1).codeBin = num2str([0 0]);
designStruct(1).lbl = 'continuous association';

designStruct(2).codeBin = num2str([0 1]);
designStruct(2).lbl = 'independent observation comparison(s)';

designStruct(3).codeBin = num2str([1 0]);
designStruct(3).lbl = 'repeated measure comparison(s)';

designStruct(4).codeBin = num2str([1 1]);
designStruct(4).lbl = 'indepdendent and repeated measure comparison(s)';

designOptions = struct2table(designStruct);

% transform code from binary to decimal
designOptions.codeDec = nan(size(designOptions,1),1);
for rIdx = 1:size(designOptions,1)
    designOptions{rIdx,"codeDec"} = bin2dec(char(designOptions.codeBin(rIdx,:)));
end
designOptions = sortrows(designOptions(:,[3 1 2]),"codeDec");
