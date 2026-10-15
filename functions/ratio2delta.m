function deltaD = ratio2delta(r)
%RATIO2DELTA Convert normalized isotope ratio to delta2H (per mil).

deltaD = 1000 .* (r - 1);

end