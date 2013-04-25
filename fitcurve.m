function[estimates,model]=fitcurve(xdata,ydata,start_point)
model=@expfun;
estimates=fminsearch(model,start_point);
  function[sse,FittedCurve]=expfun(params) 
        A=params(1);
        tr=params(2);
        to=params(3);
        do=params(4);
        U=zeros(size(xdata-to)); 
        U((xdata-to)>=0)=1; 
        FittedCurve=(1- exp(-(xdata - to)./tr))*A; 
        FittedCurve=FittedCurve.*U+do; 
        ErorVector=FittedCurve-ydata; 
        sse=sum(ErorVector.^2);
   end
end