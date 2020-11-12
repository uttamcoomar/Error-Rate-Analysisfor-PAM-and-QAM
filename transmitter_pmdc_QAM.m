function [f]=transmitter_pmdc_QAM(M,Es)
A=sqrt((3*Es)/((M)-1));
symbol_tr=zeros(10000,2);                                                              
      for i=1:1:10000
        cache_i=generate_QAM(M);
        symbol_tr(i,1)=cache_i(1);
        symbol_tr(i,2)=cache_i(2);
      end
symbol_tr=A.*symbol_tr;
f=symbol_tr;
end