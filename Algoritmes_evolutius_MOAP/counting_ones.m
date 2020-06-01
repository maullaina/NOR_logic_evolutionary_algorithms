function f = counting_ones( a )
%   counting_ones(a) a: vector{0,1}
%   Ex:  counting_ones(rand(1,10)> 0.5)
    f = sum(a > 0);
end

