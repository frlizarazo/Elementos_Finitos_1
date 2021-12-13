function A=cono_Truncado(x)
    if  x<1
        r=-0.03*x+0.045;
        A=pi*r^2;
    else
        A=pi*0.015^2;
    end
end