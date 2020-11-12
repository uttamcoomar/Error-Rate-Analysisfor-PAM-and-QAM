prompt = {'Enter M:','PAM or QAM ?','SER comparison or BER comparison?','(For BER only) SNR per symbol or per bit'};
dlg_title = '  ';
num_lines = 1;
defaultans = {'2','PAM','SER','SNR per bit   '};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
M=str2num(vectorize(answer(1)));
choice=vectorize(answer(2));
choice1=vectorize(answer(3));
choice2=vectorize(answer(4));
Es=1;
Dbmax=40;
ser_sim=[];
count1=0;
if choice1=='SER'
switch choice
    case 'PAM'
        
for iEs=0:1:Dbmax
    count=0;
    symbol_tr=transmitter_pmdc_PAM(M,10^(iEs/10));
    rx_inner=zeros(1,length(symbol_tr));
    rx_outer_1=zeros(1,length(symbol_tr));
    rx_outer_2=zeros(1,length(symbol_tr));
    A=sqrt((3*(10^(iEs/10)))/((M^2)-1));
    symbol_tr_norm=round(symbol_tr/(abs(A)));
    symbol_rx=symbol_tr+(1/sqrt(2)*[randn(1,length(symbol_tr_norm))]);
    symbol_rx_norm=symbol_rx/(abs(A));
    for i=-(M-1):2:(M-1)
        if i==-(M-1)
            rx_outer_1(((symbol_rx_norm)<=-(M-2)))=-(M-1);
        end
        if i==(M-1)
            rx_outer_2(((symbol_rx_norm)>(M-2)))=(M-1);
        end
        if i~=-(M-1) && i~=(M-1)
            if i<0
                rx_inner((symbol_rx_norm)<=i+1 & (symbol_rx_norm)>i-1)=i;
            end
            if i>0
                rx_inner((symbol_rx_norm)<i+1 & (symbol_rx_norm)>=i-1)=i;
            end
        end
    end
    symbol_rx_norm=rx_inner+((rx_outer_1)+(rx_outer_2));
    count=nnz(symbol_rx_norm-symbol_tr_norm);
    count1=count1+1;
    ser_sim(count1)=count/(length(symbol_rx_norm));
end
p1=semilogy(0:1:Dbmax,ser_sim,'*','LineWidth',2,'MarkerSize',10,'DisplayName',sprintf('%d %s simulation SER',M,choice));
hold on
p2=semilogy(0:1:Dbmax,2*((M-1)/M)*(qfunc(sqrt((6/((M^2)-1))*(10.^((0:1:Dbmax)/10))))),'LineWidth',2,'DisplayName',sprintf('%d %s theoretical SER',M,choice));
hold on
xlabel('SNR in dB -->');
ylabel('Symbol Error Rate');
hold on
axis([0 Dbmax 10^-5 1]);
hold on
grid on
legend('-DynamicLegend');
hold all;

case 'QAM'
   for iEs=0:1:Dbmax
    count=0;
    f=transmitter_pmdc_QAM(M,(10^(iEs/10))/2);
    symbol_tr_i=f(:,1)';
    rx_inner_i=zeros(1,length(symbol_tr_i));
    rx_outer_1_i=zeros(1,length(symbol_tr_i));
    rx_outer_2_i=zeros(1,length(symbol_tr_i));
    A=sqrt((3*(10^(iEs/10))/2)/((sqrt(M)^2)-1));
    symbol_tr_norm_i=round(symbol_tr_i/(abs(A)));
    symbol_rx_i=symbol_tr_i+(1/sqrt(2)*randn(1,length(symbol_tr_norm_i)));
    symbol_rx_norm_i=symbol_rx_i/(abs(A));
    for i=-(sqrt(M)-1):2:(sqrt(M)-1)
        if i==-(sqrt(M)-1)
            rx_outer_1_i(((symbol_rx_norm_i)<=-(sqrt(M)-2)))=-(sqrt(M)-1);
        end
        if i==(sqrt(M)-1)
            rx_outer_2_i(((symbol_rx_norm_i)>(sqrt(M)-2)))=(sqrt(M)-1);
        end
        if i~=-(sqrt(M)-1) && i~=(sqrt(M)-1)
            if i<0
                rx_inner_i((symbol_rx_norm_i)<=i+1 & (symbol_rx_norm_i)>i-1)=i;
            end
            if i>0
                rx_inner_i((symbol_rx_norm_i)<i+1 & (symbol_rx_norm_i)>=i-1)=i;
            end
        end
    end
    symbol_rx_norm_i=rx_inner_i+((rx_outer_1_i)+(rx_outer_2_i));
    count_i=(symbol_rx_norm_i-symbol_tr_norm_i);
    count_i=abs((count_i)./2);
    count1=count1+1;
    symbol_tr_q=f(:,2)';
    rx_inner_q=zeros(1,length(symbol_tr_q));
    rx_outer_1_q=zeros(1,length(symbol_tr_q));
    rx_outer_2_q=zeros(1,length(symbol_tr_q));
    A=sqrt((3*((10^(iEs/10))/2))/((sqrt(M)^2)-1));
    symbol_tr_norm_q=round(symbol_tr_q/(abs(A)));
    symbol_rx_q=symbol_tr_q+(1/sqrt(2)*[randn(1,length(symbol_tr_norm_q))]);
    symbol_rx_norm_q=symbol_rx_q/(abs(A));
    for i=-(sqrt(M)-1):2:(sqrt(M)-1)
        if i==-(sqrt(M)-1)
            rx_outer_1_q(((symbol_rx_norm_q)<=-(sqrt(M)-2)))=-(sqrt(M)-1);
        end
        if i==(sqrt(M)-1)
            rx_outer_2_q(((symbol_rx_norm_q)>(sqrt(M)-2)))=(sqrt(M)-1);
        end
        if i~=-(sqrt(M)-1) && i~=(sqrt(M)-1)
            if i<0
                rx_inner_q((symbol_rx_norm_q)<=i+1 & (symbol_rx_norm_q)>i-1)=i;
            end
            if i>0
                rx_inner_q((symbol_rx_norm_q)<i+1 & (symbol_rx_norm_q)>=i-1)=i;
            end
        end
    end
    symbol_rx_norm_q=rx_inner_q+((rx_outer_1_q)+(rx_outer_2_q));
    count_q=symbol_rx_norm_q-symbol_tr_norm_q;
    count_q=abs(count_q./2);
    ser_sim(count1)=nnz(count_i+count_q)/(length(symbol_rx_norm_q));
end
p1=semilogy(0:1:Dbmax,ser_sim,'*','LineWidth',2,'MarkerSize',5,'DisplayName',sprintf('%d %s simulation SER',M,choice));
hold on
p2=semilogy(0:1:Dbmax,4*((sqrt(M)-1)/sqrt(M))*(qfunc(sqrt((3/((M)-1))*(10.^((0:1:Dbmax)/10))))),'LineWidth',2,'DisplayName',sprintf('%d %s theoretical SER',M,choice));
hold on
axis([0 Dbmax 10^-5 1]);
hold on
grid on
legend('-DynamicLegend');
hold all
end
end

if choice1=='BER'
switch choice
    case 'PAM'
        
for iEs=0:1:Dbmax
    count=0;
    symbol_tr=transmitter_pmdc_PAM(M,10^(iEs/10));
    rx_inner=zeros(1,length(symbol_tr));
    rx_outer_1=zeros(1,length(symbol_tr));
    rx_outer_2=zeros(1,length(symbol_tr));
    A=sqrt((3*(10^(iEs/10)))/((M^2)-1));
    symbol_tr_norm=round(symbol_tr/(abs(A)));
    symbol_rx=symbol_tr+(1/sqrt(2)*[randn(1,length(symbol_tr_norm))]);
    symbol_rx_norm=symbol_rx/(abs(A));
    for i=-(M-1):2:(M-1)
        if i==-(M-1)
            rx_outer_1(((symbol_rx_norm)<=-(M-2)))=-(M-1);
        end
        if i==(M-1)
            rx_outer_2(((symbol_rx_norm)>(M-2)))=(M-1);
        end
        if i~=-(M-1) && i~=(M-1)
            if i<0
                rx_inner((symbol_rx_norm)<=i+1 & (symbol_rx_norm)>i-1)=i;
            end
            if i>0
                rx_inner((symbol_rx_norm)<i+1 & (symbol_rx_norm)>=i-1)=i;
            end
        end
    end
    symbol_rx_norm=rx_inner+((rx_outer_1)+(rx_outer_2));
    count=nnz(symbol_rx_norm-symbol_tr_norm);
    count1=count1+1;
    ser_sim(count1)=count/(length(symbol_rx_norm));
end
if choice2=='SNR per symbol'
p1=semilogy(0:1:Dbmax,ser_sim./(log2(M)),'-*','LineWidth',2,'MarkerSize',5,'DisplayName',sprintf('%d %s simulation BER',M,choice));
xlabel('SNR in dB (per symbol)-->');
end
if choice2=='SNR per bit   '
p1=semilogy((0:1:Dbmax)./(log2(M)),ser_sim./(log2(M)),'-*','LineWidth',2,'MarkerSize',5,'DisplayName',sprintf('%d %s simulation BER',M,choice));
xlabel('SNR in dB (per bit)-->');
end
hold on
axis([0 Dbmax 10^-5 1]);
hold on
ylabel('Bit Error Rate -->');
hold on
grid on
legend('-DynamicLegend');
hold all

case 'QAM'
   for iEs=0:1:Dbmax
    count=0;
    f=transmitter_pmdc_QAM(M,(10^(iEs/10))/2);
    symbol_tr_i=f(:,1)';
    rx_inner_i=zeros(1,length(symbol_tr_i));
    rx_outer_1_i=zeros(1,length(symbol_tr_i));
    rx_outer_2_i=zeros(1,length(symbol_tr_i));
    A=sqrt((3*(10^(iEs/10))/2)/((sqrt(M)^2)-1));
    symbol_tr_norm_i=round(round(symbol_tr_i/(abs(A))));
    symbol_rx_i=symbol_tr_i+(1/sqrt(2)*randn(1,length(symbol_tr_norm_i)));
    symbol_rx_norm_i=symbol_rx_i/(abs(A));
    for i=-(sqrt(M)-1):2:(sqrt(M)-1)
        if i==-(sqrt(M)-1)
            rx_outer_1_i(((symbol_rx_norm_i)<=-(sqrt(M)-2)))=-(sqrt(M)-1);
        end
        if i==(sqrt(M)-1)
            rx_outer_2_i(((symbol_rx_norm_i)>(sqrt(M)-2)))=(sqrt(M)-1);
        end
        if i~=-(sqrt(M)-1) && i~=(sqrt(M)-1)
            if i<0
                rx_inner_i((symbol_rx_norm_i)<=i+1 & (symbol_rx_norm_i)>i-1)=i;
            end
            if i>0
                rx_inner_i((symbol_rx_norm_i)<i+1 & (symbol_rx_norm_i)>=i-1)=i;
            end
        end
    end
    symbol_rx_norm_i=rx_inner_i+((rx_outer_1_i)+(rx_outer_2_i));
    count_i=(symbol_rx_norm_i-symbol_tr_norm_i);
    count_i=abs((count_i)./2);
    count1=count1+1;
    symbol_tr_q=f(:,2)';
    rx_inner_q=zeros(1,length(symbol_tr_q));
    rx_outer_1_q=zeros(1,length(symbol_tr_q));
    rx_outer_2_q=zeros(1,length(symbol_tr_q));
    A=sqrt((3*((10^(iEs/10))/2))/((sqrt(M)^2)-1));
    symbol_tr_norm_q=round(symbol_tr_q/(abs(A)));
    symbol_rx_q=symbol_tr_q+(1/sqrt(2)*[randn(1,length(symbol_tr_norm_q))]);
    symbol_rx_norm_q=symbol_rx_q/(abs(A));
    for i=-(sqrt(M)-1):2:(sqrt(M)-1)
        if i==-(sqrt(M)-1)
            rx_outer_1_q(((symbol_rx_norm_q)<=-(sqrt(M)-2)))=-(sqrt(M)-1);
        end
        if i==(sqrt(M)-1)
            rx_outer_2_q(((symbol_rx_norm_q)>(sqrt(M)-2)))=(sqrt(M)-1);
        end
        if i~=-(sqrt(M)-1) && i~=(sqrt(M)-1)
            if i<0
                rx_inner_q((symbol_rx_norm_q)<=i+1 & (symbol_rx_norm_q)>i-1)=i;
            end
            if i>0
                rx_inner_q((symbol_rx_norm_q)<i+1 & (symbol_rx_norm_q)>=i-1)=i;
            end
        end
    end
    symbol_rx_norm_q=rx_inner_q+((rx_outer_1_q)+(rx_outer_2_q));
    count_q=symbol_rx_norm_q-symbol_tr_norm_q;
    count_q=abs(count_q./2);
    ser_sim(count1)=nnz(count_i+count_q)/(length(symbol_rx_norm_q));
   end
if choice2=='SNR per symbol'
p1=semilogy(0:1:Dbmax,ser_sim./(log2(M)),'-*','LineWidth',2,'MarkerSize',5,'DisplayName',sprintf('%d %s simulation BER',M,choice));
xlabel('SNR in dB (per symbol)-->');
end
if choice2=='SNR per bit   '
p1=semilogy((0:1:Dbmax)./(log2(sqrt(M))),ser_sim./(log2(M)),'-*','LineWidth',2,'MarkerSize',5,'DisplayName',sprintf('%d %s simulation BER',M,choice));
xlabel('SNR in dB (per bit)-->');
end
hold on
axis([0 Dbmax 10^-5 1]);
hold on
ylabel('Bit Error Rate -->');
hold on
grid on
legend('-DynamicLegend');
hold all
end
end