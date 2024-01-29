function E = SourceFct2(t, InputParas)

if isfield(InputParas, 'rep') %if rep is a name of field InputParas
    n = floor(t / InputParas.rep);
    t = t - n * InputParas.rep;
end

if ~isstruct(InputParas) %if InputParas is not struct
    E = InputParas;
else
    %if InputParas is a struct
    % Extract parameters from InputParas structure 
    % (used to make InputParas R a struct)
    E0 = InputParas.E0;
    t0 = InputParas.t0;
    wg = InputParas.wg;
    we = InputParas.we;
    phi = InputParas.phi;

    % Calculate the wave using the extracted parameters
    %E = E0 * exp(-(t - t0)^2 / wg^2) * exp(1i * (we * t + phi));
    E = E0 * exp(-(t - t0)^2 / wg^2) * exp(1i * (we * t + phi));
end

end