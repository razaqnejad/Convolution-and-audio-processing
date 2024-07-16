function y = trapezoidal_rule(func, ll, ul, n)
    dx = (ul-ll)/n;
    x = linspace(ll, ul, n+1);
    y = sum(func(x) .* ([1, 2*ones(1,n-1), 1] * dx/2));
end
