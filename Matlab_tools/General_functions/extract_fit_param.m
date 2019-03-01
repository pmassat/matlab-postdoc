function [prm, prmStdErr] = extract_fit_param(data_struct,dsfield,paramNum)
    cfn = coeffnames(data_struct.dsfield);
    prm = extract_structure_field(data_struct,dsfield,cfn{paramNum});
    prmErr = abs(prm-cftm(:,paramNum));errAp = abs(prm-cftp(:,paramNum));% negative and positive
    relErrAn = prmErr./prm;relErrAp = errAp./prm;% negative and positive
end
