function f=generate_QAM(M)
A=[-(sqrt(M)-1):2:(sqrt(M)-1)];
fi=randsrc(1,1,A);
fq=randsrc(1,1,A);
f=[fi fq];
end