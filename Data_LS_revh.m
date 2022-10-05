function [A_est,D_est] = Data_LS_revh(B, H_est,Y_data_rx, cwsbl,NUM_PILOT,acUsrIndx,acData,par)
NUM_ALL_USERS = par.NUM_ALL_USERS;
NUM_RX = par.NUM_RX;
nd = par.nd;
nEvryUsrDataBits = par.nEvryUsrDataBits;
nOfTfcb = par.nOfTfcb;
Mod_Type=par.Mod_Type;
trellis = par.trellis;
mRS = par.mRS;
nRS = par.nRS ;
kRS = par.kRS;
nEvryUsrIDBits = par.nEvryUsrIDBits; 


if isempty(B)
    A_est=[];
    D_est=[];
    return;
end 

rxData=zeros(NUM_PILOT,nEvryUsrDataBits);
HX_pilot=zeros(NUM_PILOT,nd*NUM_RX);
D_mulrx = zeros(NUM_PILOT,nd*NUM_RX);
D= zeros(NUM_PILOT,nd);
eat_Data = zeros(NUM_ALL_USERS,nEvryUsrDataBits);
H_rev=zeros(NUM_PILOT,NUM_RX);

HX_pilot(B,:)= cwsbl(:,B) \ Y_data_rx;

H2 = mat2cell(H_est,[size(H_est,1)],NUM_RX);
H2 = repelem(H2,1,nd);
H2 = cell2mat(H2);

D_mulrx(B,:) = HX_pilot(B,:)./ H2(B,:);

D(B,:) = reshape(sum(reshape(D_mulrx(B,:).',NUM_RX,[]))/NUM_RX ,[],length(B)).';

%====================Demodulation====================     
if Mod_Type==1 % QPSK
    all_distance(1,:,:)=abs(D(B,:)-A(1));
    all_distance(2,:,:)=abs(D(B,:)-A(2));   
    all_distance(3,:,:)=abs(D(B,:)-A(3));      
    all_distance(4,:,:)=abs(D(B,:)-A(4));      
    [~,index]=min(all_distance);
    mapping=[0 1 2 3];
    bits_Aftr_DeMod=reshape(mapping(index),[length(B),nd]); 
    dec = de2bi(reshape(bits_Aftr_DeMod',[],1),'left-msb');
    bin=reshape(dec',6,[]);                 
elseif Mod_Type==2 % 16QAM 
    bits_Aftr_DeMod = reshape(D(B,:)',[],1);
    bits_Aftr_DeMod = qamdemod(bits_Aftr_DeMod,16);      
    bin = de2bi(bits_Aftr_DeMod',4,'left-msb');
    bin = reshape(bin',4,[]);                  
elseif Mod_Type==3 % 64QAM
    bits_Aftr_DeMod = reshape(D(B,:)',[],1);
    bits_Aftr_DeMod = qamdemod(bits_Aftr_DeMod,64);      
    bin = de2bi(bits_Aftr_DeMod',6,'left-msb');
    bin = bin';
end    

%%%%%%%%%%%%%%%%�ŵ�У��%%%%%%%%%%%%%%%%%%%%
if Mod_Type==1 % QPSK
    %����xΪ0 1 2 3  ���yΪ1+1i 1-1i -1+1i -1-1i
    dec = bi2de(bin','left-msb');
    sbl_Aftr_Mod = zeros(NUM_ALL_USERS,nd);
    sbl_Aftr_Mod(B,:) =reshape(A(dec(:,:)+1),nd,[]).';                
elseif Mod_Type==2 % 16QAM 
    dec = bi2de(bin','left-msb'); % Binary to decimal conversion
    sbl_Aftr_Mod = zeros(NUM_ALL_USERS,nd);
    sbl_Aftr_Mod(B,:) = reshape(qammod(dec,16),nd,[])';                   
elseif Mod_Type==3 % 64QAM
    dec = bi2de(bin','left-msb'); % Binary to decimal conversion
    sbl_Aftr_Mod = zeros(NUM_ALL_USERS,nd);
    sbl_Aftr_Mod(B,:) = reshape(qammod(dec,64),nd,[])';                
end   

Xsbl = repelem(sbl_Aftr_Mod,1,NUM_RX);

H_rev(B,:)=reshape(sum(reshape( HX_pilot(B,:)./ Xsbl(B,:),[],nd),2)/nd,[],NUM_RX);

H2 = mat2cell(H_rev,[size(H_est,1)],NUM_RX);
H2 = repelem(H2,1,nd);
H2 = cell2mat(H2);

D_mulrx(B,:) = HX_pilot(B,:)./ H2(B,:);

D(B,:) = reshape(sum(reshape(D_mulrx(B,:).',NUM_RX,[]))/NUM_RX ,[],length(B)).';

%====================Demodulation====================     
if Mod_Type==1 % QPSK
    all_distance(1,:,:)=abs(D(B,:)-A(1));
    all_distance(2,:,:)=abs(D(B,:)-A(2));   
    all_distance(3,:,:)=abs(D(B,:)-A(3));      
    all_distance(4,:,:)=abs(D(B,:)-A(4));      
    [~,index]=min(all_distance);
    mapping=[0 1 2 3];
    bits_Aftr_DeMod=reshape(mapping(index),[length(B),nd]); 
    dec = de2bi(reshape(bits_Aftr_DeMod',[],1),'left-msb');
    bin=reshape(dec',6,[]);                 
elseif Mod_Type==2 % 16QAM 
    bits_Aftr_DeMod = reshape(D(B,:)',[],1);
    bits_Aftr_DeMod = qamdemod(bits_Aftr_DeMod,16);      
    bin = de2bi(bits_Aftr_DeMod',4,'left-msb');
    bin = reshape(bin',4,[]);                  
elseif Mod_Type==3 % 64QAM
    bits_Aftr_DeMod = reshape(D(B,:)',[],1);
    bits_Aftr_DeMod = qamdemod(bits_Aftr_DeMod,64);      
    bin = de2bi(bits_Aftr_DeMod',6,'left-msb');
    bin = bin';
end   

% Deinterleaving
bits_Aftr_DeItrlv2 = matdeintrlv(bin,2,2);
bits_Aftr_DeItrlv(B,:) = reshape(bits_Aftr_DeItrlv2,[],length(B))';

% Convolutionally Decoding
for q = 1:length(B)
    decodeData(B(q),:) = vitdec(bits_Aftr_DeItrlv(B(q),:),trellis,5,'trunc','hard');
end

% RS decoding
decodeData2 = reshape(decodeData(B,:)',[],1);
rxRSmatrix = reshape(decodeData2,mRS,length(decodeData2)/mRS)';
rxRStrSbl = bi2de(rxRSmatrix,'left-msb');
hDec = comm.RSDecoder(nRS,kRS);
rxRS_decoded_Sbl = step(hDec, rxRStrSbl);
rxRS_decoded_Bits = de2bi(rxRS_decoded_Sbl',3,'left-msb');
rxData(B,:) = reshape(rxRS_decoded_Bits',nEvryUsrDataBits,length(B))';
usrID_mtx(B,:) = bi2de(rxData(B,nEvryUsrDataBits-nEvryUsrIDBits+1:nEvryUsrDataBits),'left-msb');
A_est = usrID_mtx(B,:);
B(find(A_est==0))=[];
A_est(find(A_est==0))=[];
eat_Data(A_est,:)= rxData(B,:);
D_est = acUsrIndx(find(all(bsxfun(@eq,eat_Data(acUsrIndx,:),acData(acUsrIndx,:)),2)));








