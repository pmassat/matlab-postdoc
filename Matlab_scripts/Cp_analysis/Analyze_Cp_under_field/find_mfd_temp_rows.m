function [rows, wt, varargout] = find_mfd_temp_rows(t_calc, utmfd, Tc, Tmfd, h, varargin )

%% Parse function arguments
    defaultTref = [0 0];
    defaultPrintTref = false;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    validArrayPosNum = @(x) isnumeric(x) && all(x >= 0);
    
    p = inputParser;
    addRequired(p, 't_calc', validScalarPosNum);
    addRequired(p, 'utmfd', validArrayPosNum);
    addRequired(p, 'Tc', validScalarPosNum);
    addRequired(p, 'Tmfd', @istable);
    addRequired(p, 'h', validScalarPosNum);
    addOptional(p, 'tref', defaultTref, validArrayPosNum);
    addParameter(p, 'printTref', defaultPrintTref, @islogical);

    parse(p, t_calc, utmfd, Tc, Tmfd, h, varargin{:});
    
    tref = p.Results.tref;
    printTref = p.Results.printTref;

%% Function code
    % Find value of temperature in COMSOL mfd closest to that of actual data
%     utmfd = unique(Tmfd.T_K);
    Trcomp = utmfd/Tc - t_calc;
    [abst, mfdtidx] = sort(abs(Trcomp));

    % if the 2 closest temperatures are both above or both below the one of
    % interest, just use the single closest, otherwise use both
    if sign(Trcomp(mfdtidx(1)))==sign(Trcomp(mfdtidx(2)))
        t(1:2) = utmfd(mfdtidx(1));
        wt = [1,0];
    else
        t = utmfd(mfdtidx(1:2));% 2 closest temperatures
        wt = 1-abst(1:2)/sum(abst(1:2));% weight is 1 minus relative difference in temperature
    end

    % Check that the above code identifies the two closest temperatures correctly
    if printTref==true
        if ~all(t==tref)
            sprintf('T=%.2gK, Tref=[%.2g,%.2g]K',...
                t_calc*Tc, t)
        end
        varargout{1} = t;% output new value of tref as varargout{1}
    end

    % Find the row in Tmfd that matches both t and h 
    rows = find(ismember(Tmfd.T_K,t) & Tmfd.Hext_Oe==h);

    % if a row contains NaN, ignore this row, and only use the other 
    % closest temperature, with a weight of 1
    if any(isnan(Tmfd.hc(rows,:)),'all')
        rows = rows(~any(isnan(Tmfd.hc(rows,:)),2));
        wt = [1,0];
    end

end