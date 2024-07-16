function y = upper_riemann(func, ll, ul, n)
    dx = (ul-ll)/n;
    x = linspace(ll, ul, n+1);
    y = sum(max(func(x(1:end-1)), func(x(2:end))) * dx);
end