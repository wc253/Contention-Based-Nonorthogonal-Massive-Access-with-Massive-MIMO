function [A_est,B_est, D_est,H_est]= LTE_FOUR_STEP_MIMO(Y_pilot_rx,cwplt,HX,cwsbl,acData,acUsrIndx,pilot_acUsr,acUsr_pilot,par,pilotIndx,pilot_choose)
nd = par.nd;
mRS = par.mRS;
nRS = par.nRS;
kRS = par.kRS;
trellis = par.trellis;
nEvryUsrDataBits = par.nEvryUsrDataBits;
nEvryUsrIDBits = par.nEvryUsrIDBits; 
SNR = par.SNR;
NUM_RX = par.NUM_RX ;
nOfTfcb = par.nOfTfcb;
Mod_Type=par.Mod_Type;
A=par.A;

[~,n]=size(cwplt);
H_est = zeros(n,size(Y_pilot_rx,2));
% H_est=cwplt\Y_pilot_rx; 
% [H_sort,pos]=sort(sum(abs(H_est).^2,2),'descend');
% 
% err=10^(-SNR/10);
% B_est=pos(H_sort>err);
% 
% H_re_est=cwplt(:,B_est)\Y_pilot_rx; 
% H_est(B_est,:)=H_re_est;
B_est=pilot_choose;
H_est(B_est,:)=cwplt(:,B_est)\Y_pilot_rx;


% pilot_choose =sort(unique(pilotIndx));
% fprintf('num_est = %d,nun_real=%d,\n',length(B_est),length(pilot_choose));
% fprintf('acc = %d, \n',length(intersect(B_est,pilot_choose)) / length(pilot_choose));
A_est=[];
D_est=[];

for pilot_index = 1:length(B_est)
    H_user= H_est(B_est(pilot_index),:);
    H_user = mat2cell(H_user,size(H_user,1),NUM_RX);
    H_user = repelem(H_user,1,nd);
    H_user = cell2mat(H_user);   
    
    
    cwsbl_user=cwsbl(:,B_est(pilot_index));
    HX_multiuser=zeros(1,nd*NUM_RX);
    user_thispilot=cell2mat(pilot_acUsr(B_est(pilot_index)));
    for user_index = length(user_thispilot)
        HX_multiuser=HX_multiuser+HX(user_thispilot(user_index),:);
    end
    
    UserData = cwsbl_user*HX_multiuser;
    UserData = awgn(UserData,SNR,'measured');
    HX_est = cwsbl_user \ UserData;
    
    D_user= HX_est./H_user;
    D_user=sum(reshape(D_user,NUM_RX,[]))/NUM_RX;
    %====================Demodulation====================     
    if Mod_Type==1 % QPSK
        all_distance(1,:,:)=abs(D_user-A(1));
        all_distance(2,:,:)=abs(D_user-A(2));   
        all_distance(3,:,:)=abs(D_user-A(3));      
        all_distance(4,:,:)=abs(D_user-A(4));      
        [~,index]=min(all_distance);
        mapping=[0 1 2 3];
        bits_Aftr_DeMod=reshape(mapping(index),[1,nd]); 
        dec = de2bi(reshape(bits_Aftr_DeMod',[],1),'left-msb');
        bin=reshape(dec',6,[]);                 
    elseif Mod_Type==2 % 16QAM 
        bits_Aftr_DeMod = reshape(D_user',[],1);
        bits_Aftr_DeMod = qamdemod(bits_Aftr_DeMod,16);      
        bin = de2bi(bits_Aftr_DeMod',4,'left-msb');
        bin = reshape(bin',4,[]);                  
    elseif Mod_Type==3 % 64QAM
        bits_Aftr_DeMod = reshape(D_user',[],1);
        bits_Aftr_DeMod = qamdemod(bits_Aftr_DeMod,64);      
        bin = de2bi(bits_Aftr_DeMod',6,'left-msb');
        bin = bin';
    end    

    % Deinterleaving
    bits_Aftr_DeItrlv2 = matdeintrlv(bin,2,2);
    bits_Aftr_DeItrlv = reshape(bits_Aftr_DeItrlv2,[],1)';

    decodeData = vitdec(bits_Aftr_DeItrlv,trellis,5,'trunc','hard');

    % RS decoding
    decodeData2 = reshape(decodeData',[],1);
    rxRSmatrix = reshape(decodeData2,mRS,length(decodeData2)/mRS)';
    rxRStrSbl = bi2de(rxRSmatrix,'left-msb');
    hDec = comm.RSDecoder(nRS,kRS);
    rxRS_decoded_Sbl = step(hDec, rxRStrSbl);
    rxRS_decoded_Bits = de2bi(rxRS_decoded_Sbl',3,'left-msb');
    rxData = reshape(rxRS_decoded_Bits',nEvryUsrDataBits,1)';
    usrID_mtx = bi2de(rxData(nEvryUsrDataBits-nEvryUsrIDBits+1:nEvryUsrDataBits),'left-msb');
    A_est = [A_est usrID_mtx];
    
    D_est_user = acUsrIndx(find(all(bsxfun(@eq,rxData,acData(user_thispilot,:)),2))*user_index);
    D_est = [D_est D_est_user];
end





         












