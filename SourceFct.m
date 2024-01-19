function E = SourceFct(t,InputParas)

if isfield(InputParas,'rep') %if rep is a name of field InputParas
    n = floor(t/InputParas.rep);
    t = t-n*InputParas.rep;
end

if ~isstruct(InputParas) %if InputParas is not struct
    E = InputParas; 
else
    E = InputParas.E0*sin(-(t-InputParas.t0)^2/InputParas.wg^2)*sin(1i*(InputParas.we*t + InputParas.phi));
end

end