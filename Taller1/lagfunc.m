function LaG=lagfunc(nef,nno_EF)
    nno  = (nno_EF-1)*nef + 1;
    LaG=zeros(nef,nno_EF);
    for i=1:nno_EF
        n=i-1;
        LaG(:,i)  = LaG(:,i)+(1+n:nno_EF-1:(nno-nno_EF+i))';
    end
end