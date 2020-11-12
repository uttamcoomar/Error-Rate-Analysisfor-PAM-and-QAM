function [f]=transmitter_pmdc_PAM(M,Es)
A=sqrt((3*Es)/((M^2)-1));
symbol_tr=zeros(1,10000);                                                              
      for i=1:1:10000
        symbol_tr(i)=generate_PAM(M);
      end
symbol_tr=A.*symbol_tr;
f=symbol_tr;
end