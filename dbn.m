function [y] = dbn(x,mindb) 
% [y] = db(x,mindb) 
% convert x to db normalized to 0dB max and optionally clipped to mindb - see also db(), undb()
   mx = max(abs(x));
   if nargin > 1,
      if (mindb > 0), mindb = - mindb; end;
   else
      mindb = -1000;
   end;
   if (mx >= eps), 
    x = x/mx;
   end;
   scl = (10/log(10)); 
   arg = x .* conj(x); 
   [m,n] = size(arg); 
   arg = max(arg,1E-20);
   if m==1|n==1, tmp = scl*log(arg); 
   else 
         tmp=ones(m,n); 
        if m>n, for i=1:n, tmp(:,i) = scl*log(arg(:,i)); end;
        else for i=1:m, tmp(i,:) = scl*log(arg(i,:)); end; 
        end;
   end; 
   minm = mindb * ones(m,n);
   y = max(tmp,minm);   