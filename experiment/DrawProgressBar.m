function f = DrawProgressBar(screenparms, nMax,str)
%ProgressBar    Draw Ascii progress bar on Screen.
%   progBar = DrawProgressBar(nMax, str) creates a progress bar and returns a
%   pointer to a function handle which can then be called to update it.
%
%   To update, call progBar(currentStep)
%
%   Example:
%      n = 500;
%      progBar = ProgressBar(n,'computing...');
%      for tmp = 1:n
%        progBar(tmp);
%        pause(.01)
%      end

%   by David Szotten 2008
%   $Revision: 1.2 $  $Date: 2008/04/17 09:15:32 $
%   merged with utility by us / CSSM by: DN 2008
%  2008-09-16  DN Added elapsed time and estimated time left
%  2009-10-06  DL Adapted for Psychtoolbox

if nargin>1
    head = sprintf('%s\n', str);
else
    head = '';
end

lastPercentileWritten = 0;
pstrlen               = 0;
tstrlen               = 0;
pstr = '';

label = sprintf('| 0%s50%s100%% |\n',repmat(' ',1,21),repmat(' ',1,18));

Screen(screenparms.window, 'DrawText', sprintf('%s', head), 50, 50);
Screen(screenparms.window, 'DrawText', sprintf('%s', label), 50, 100);
% hllen = length([head label]);

t     = datenum(clock);

f = @updateBar;
    function updateBar(nCurrent)
        Screen(screenparms.window, 'DrawText', sprintf('%s', head), 50, 50);
        Screen(screenparms.window, 'DrawText', sprintf('%s', label), 50, 100);
        
        %what percentile are we up to
        currentPercentile = round(50 * nCurrent/nMax);
        
%         Screen(screenparms.window, 'DrawText', sprintf('%s',repmat(char(8),1,tstrlen)), 50, 150); % remove time string
        
        % compute time info
        ttn = datenum(clock)-t;
        tt  = datevec(ttn);
        dtt = ttn/nCurrent;
        ttleft = datevec(dtt*(nMax-nCurrent));
        tstr1 = sprintf('\nElapsed time:        %dh %dm %ds\n ', tt(4),tt(5),round(tt(6)));
        tstr2 = sprintf('Estimated time left: %dh %dm %ds',ttleft(4),ttleft(5),round(ttleft(6)));
%         tstrlen = length(tstr1) + length(tstr2);

        %have we passed another percentile?
        if (currentPercentile > lastPercentileWritten )

            %we may have increased by several percentiles,
            %so keep writing until we catch up
            percentileToWrite = lastPercentileWritten + 1;
            while(percentileToWrite <= currentPercentile)
                %for every 10th, use a '+' instead
                if( mod(percentileToWrite,5) == 0 )
                    pstr = [pstr, sprintf('%s','+')];
%                         Screen(screenparms.window, 'DrawText', sprintf('%s','+'), 50, 150);
                else
                    pstr = [pstr, sprintf('%s','.')];
%                     Screen(screenparms.window, 'DrawText', sprintf('%s','.'), 50, 150);
                end
                percentileToWrite = percentileToWrite + 1;
                pstrlen = pstrlen + 1;
            end
            Screen(screenparms.window, 'DrawText', pstr, 50, 150);
            
            %update status
            lastPercentileWritten = currentPercentile;

            % write time string
            Screen(screenparms.window, 'DrawText', sprintf('%s',tstr1), 50, 250);
            % write time string
            Screen(screenparms.window, 'DrawText', sprintf('%s',tstr2), 50, 300);
        else
            % write time string
            Screen(screenparms.window, 'DrawText', sprintf('%s',tstr1), 50, 250);
            % write time string
            Screen(screenparms.window, 'DrawText', sprintf('%s',tstr2), 50, 300);
        end
        
        
%         %are we done?
        if nCurrent==nMax
            %clear bar
%             pause(1)
    Screen('FillRect', screenparms.window, screenparms.color);
    Screen('Flip', screenparms.window); 
        end
    end
end

