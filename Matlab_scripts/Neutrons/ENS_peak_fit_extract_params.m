function [g,relErrGn,varargout] = ENS_peak_fit_extract_params(data_struct,dsfield,varargin)
% function [g,relErrGn,a,relErrAn,s,relErrSn] = ENS_peak_fit_extract_params(data_struct,dsfield,varargin)
    if isfield(data_struct,dsfield)
% dsfield should be of type 'cfit' or 'sfit'
    message1 = "Positive and negative error bars are ";
    message2 = "the same for parameter ";

%% Compute confidence intervals of fit
    % dfe = cell2mat( arrayfun(@(c) c.gof.dfe, data_struct.', 'Uniform', 0) );
    % levelCfd = 2*tcdf(-1,dfe) is the level of confidence of the fit
    % I dont see why one should use this quantity to compute error bars, as in
    % cft = confint(data_struct(i).fitresult,1-2*tcdf(-1,data_struct(i).gof.dfe))% confidence intervals of fit parameters for a single dataset
    % which was suggested on a Matlab forum webpage,
    % so I simply use the default 95% confidence bounds
    cft = cell2mat(arrayfun(@(c) confint(c.(dsfield)),data_struct.','Uniform',0));
    % 95% confidence bounds of all fitted datasets
    cftm = cft(1:2:end,:);% extract lower bounds of interval
    cftp = cft(2:2:end,:);% extract upper bounds of interval

%% Extract fit parameter gamma
    g = extract_structure_field(data_struct,dsfield,'gamma');

%% Extract fit parameter alpha if it exists and compute error on gamma
        try % if the fit object has a parameter called "alpha", extract this parameter
            a = extract_structure_field(data_struct,dsfield,'alpha');
            varargout{1} = a;
            errAn = abs(a-cftm(:,2));errAp = abs(a-cftp(:,2));% negative and positive
            relErrAn = errAn./a;relErrAp = errAp./a;% negative and positive
            varargout{2} = relErrAn;
% along with error on gamma, which will come in third position in the cfit
% object and hence the corresponding erros bars will be in 3rd position in
% the confidence interval
            % Use confidence intervals to determine the error on the fit parameters
            errGn = abs(g-cftm(:,3));errGp = abs(g-cftp(:,3));% standard errors
            relErrGn = errGn./g;relErrGp = errGp./g;% relative errors
            if nargin>2
% Check that negative and positive errors on alpha are equal when rounded to the 6th 
% digit in relative value (they can be different if not rounded due to
% numerical discretization)
                if isequal(round(relErrAn,6),round(relErrAp,6))%
                    disp(message1 + message2 + "alpha")
                else warning(message1 + "NOT " + message2 + "alpha")
                end
            end
        catch % otherwise extract only parameter gamma, which then comes in second position in the cfit object
            errGn = abs(g-cftm(:,2));errGp = abs(g-cftp(:,2));% standard errors
            relErrGn = errGn./g;relErrGp = errGp./g;% relative errors
        end

        if nargin>2
% Check that negative and positive errors on gamma are equal when rounded to the 6th 
% digit in relative value (they can be different if not rounded due to
% numerical discretization)
            if isequal(round(relErrGn,6),round(relErrGp,6))%
                disp(message1 + message2 + "gamma")
            else warning(message1 + "NOT " + message2 + "gamma")
            end
        end
        
%% Same for parameter sigma, if it exists
        try % if the fit object has a parameter called "sigma", extract this parameter
            s = extract_structure_field(data_struct,dsfield,'sigma');
            varargout{3} = s;
            errSn = abs(s-cftm(:,4));errSp = abs(s-cftp(:,4));
            relErrSn = errSn./s;relErrSp = errSp./s;
            varargout{4} = relErrSn;
            if nargin>2
% Check that negative and positive errors on sigma are equal when rounded to the 6th 
% digit in relative value (they can be different if not rounded due to
% numerical discretization)
                if isequal(round(relErrSn,6),round(relErrSp,6))%
                    disp(message1 + message2 + "sigma")
                else warning(message1 + "NOT " + message2 + "sigma")
                end
            end
        end % otherwise just ignore and continue (no 'catch' statement)

    else
        warning(['Structure does not have any field called "' dsfield '"'])
    end
end