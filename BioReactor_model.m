%matlab model

function [y] = BioReactor_model(k,yo,ym,t)
[y] = ym-(ym-yo)*exp((-k)*t);
end