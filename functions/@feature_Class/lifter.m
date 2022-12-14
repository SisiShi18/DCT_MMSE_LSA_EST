function lcep = lifter(cep,L)
    [N,D] = size(cep);
    n = 0:D-1;
    lift = 1 + (L/2)*sin(pi*n/L);
    lcep = cep .* repmat(lift,N,1);
end