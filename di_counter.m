function di_counter(counter,total)
% SYNTAX:
%   di_counter(counter,total)
%
% DESCRIPTION:
%   Utility function displaying progress during iteratiions (e.g., Monte Carlo loops).
%   Prints iteration counts at logarithmically spaced intervals (to fit the screen most of the times), estimates time to
%   completion (ETA) at iteration 100 for relatively long runs (more than 1000 iterations), and signals
%   completion. Uses persistent variables to track state across calls.
%
% INPUT:
%   counter - [numeric] current iteration number (1, 2, ..., total)
%   total   - [numeric] total number of iterations expected
%
% OUTPUT:
%   (none; prints to command window)
%
% NOTE: function used in both Diagonale and Polar
%
% AUTHOR:
%   Germano Gallicchio (germano.gallicchio@gmail.com)

    persistent counting startTime

    if total > 1
        % set persistent variables

        % display message
        if counter==1
            fprintf(['\\ counting (of ' num2str(total) ' total): '])
        end

        % start counting
        if counter==1  &&  total>=1000
            counting = true;
            startTime = tic;
        end

        % display ETA
        if counter==100
            if counting
                elapsedTime = toc(startTime);
                remainingTime = (elapsedTime*total/counter) - elapsedTime;
                remainingMinutes = floor(remainingTime / 60);
                remainingSeconds = mod(remainingTime, 60);
                fprintf([ '     ' char(9203) ' ETA ' num2str(remainingMinutes) 'min, ' num2str(round(remainingSeconds)) 'sec ' char(9203) '     '])
                counting = false;
            end
        end

        % display the count
        nSteps = 20;
        if any(counter==round(logspace(log10(1),log10(total),nSteps)))
            fprintf([num2str(counter) ' '])
        end


        % display end of count
        if counter==total; fprintf(' ...completed \n'); end
    end
end
