function [g,relErrGn,varargout] = ENS_peak_fit_extract_params(data_struct,dsfield,varargin)
% function [g,relErrGn,a,relErrAn,s,relErrSn] = ENS_peak_fit_extract_params(data_struct,dsfield,varargin)
    if isfield(data_struct,dsfield)
% dsfield should be of type 'cfit' or 'sfit'

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
        % Use confidence intervals to determine the error on the fit parameters
        errGn = abs(g-cftm(:,3));errGp = abs(g-cftp(:,3));% standard errors
        relErrGn = errGn./g;relErrGp = errGp./g;% relative standard errors

%% Same for parameter alpha, if it exists
        try % if the fit object does not have any parameter called "alpha", just ignore and continue
            a = extract_structure_field(data_struct,dsfield,'alpha');
            varargout{1} = a;
            errAn = abs(a-cftm(:,2));errAp = abs(a-cftp(:,2));% negative and positive
            relErrAn = errAn./a;relErrAp = errAp./a;% negative and positive
            varargout{2} = relErrAn;
        end
        
%% Same for parameter sigma, if it exists
        try % if the fit object does not have any parameter called "sigma", just ignore and continue
            s = extract_structure_field(data_struct,dsfield,'sigma');
            varargout{3} = s;
            errSn = abs(s-cftm(:,4));errSp = abs(s-cftp(:,4));
            relErrSn = errSn./s;relErrSp = errSp./s;
            varargout{4} = relErrSn;
        end
%     if isfield(data_struct,dsfield)
% % When dsfield is of the type 'gof'
%         r = extract_structure_field(data_struct,dsfield,'rsquare');
%     end
%     gamma = cell2mat( arrayfun(@(c) c.fitresult.gamma, data_struct.', 'Uniform', 0) );
%     sigma = cell2mat( arrayfun(@(c) c.fitresult.sigma, data_struct.', 'Uniform', 0) );
%     alpha = cell2mat( arrayfun(@(c) c.fitresult.alpha, data_struct.', 'Uniform', 0) );
%     rsquare = cell2mat( arrayfun(@(c) c.gof.rsquare, data_struct.', 'Uniform', 0) );
% 
% 	if nargin>2
% Check that negative and positive errors are equal when rounded to the 6th 
% digit in relative value (they can be different if not rounded due to
% numerical discretization)
%         message = "Positive and negative error bars are the same for parameter ";
%         if isequal(round(relErrAn,6),round(relErrAp,6))%
%             disp(message + "alpha")
%         end
%         if isequal(round(relErrGn,6),round(relErrGp,6))%
%             disp(message + "gamma")
%         end
%         if isequal(round(relErrSn,6),round(relErrSp,6))%
%             disp(message + "sigma")
%         end
%     end
    else
        warning(['Structure does not have any field called "' dsfield '"'])
    end
end