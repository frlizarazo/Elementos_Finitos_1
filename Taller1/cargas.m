function b=cargas(x)
    if (2-0.05<x)&&(x<2+0.05)
        b=5000*x^2+50000;
    else
        b=5000*x^2;
    end
end