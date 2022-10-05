% Beijing Jiaotong University
% Author: Bai
% Contention Based Massive Access Protocol for B5G： A Compressive Sensing Method
% Frequency-Selective Channel
% Channel bandwidth: 1.4MHz
% Subcarrier bandwidth: 15kHz
% Number of subcarriers: 6*12=72
% 128 FFT
% Protection of Bandwidth: 128-72=56

close all
clear all
clc
warning('off')

%====================Optional Parameters=====================
% Modulation:1-QPSK,2-16QAM,3-64QAM
Mod_Type = 2 ; 
% coherence time(ms) 
cohe_time=2;
%Signal to noise ratio
SNR = 20;

%======================Drawing Figure=========================
% Matrix storing point for drawing figure
NUM_ROW  = 10 ; % the first row for num of act users
drawPointMtx1 = zeros(NUM_ROW,1);%For pilot correct ratio
drawPointMtx2 = zeros(NUM_ROW,1);%For channel mse
drawPointMtx3 = zeros(NUM_ROW,1);%For Data correct ratio

%====================System Parameter=========================
% No. of slot
NUM_slot = cohe_time/0.5;
% No. of total symbols
m = 84*6*NUM_slot;
% Points of FFT
PO_FFT = 128;
% Cyclic Length
CYCLIC_LENGTH = 32;
% No. of TFCB
nOfTfcb = 6;
% No. of subcarriers
nOfSubCrr = 72;
% No. of protective bandwidth
nOfProtctBandwidth = PO_FFT-nOfSubCrr;

%====================User Parameter====================
% No.of all users
NUM_ALL_USERS = 1000000;
% Length of Pilot in each RB*NUM_slot
mp = 90;
% No. of data symbol
nd = 18; 
% Data code-word length
md = floor((m-mp)/nd); 
% No. of empty RB elements
eptymd = rem(m-mp,nd);
% transmit power(dbm) 
rho=23;
%total transmit energy
xi=sqrt((mp*rho)); 

%Number of bits per symbol in RS coding
mRS = 3;
%Number of symbols per codeword in RS coding
nRS = 6;
%Number of symbols per message in RS coding
kRS = 4;
%RS code rate
RS_CODE_RATE = kRS / nRS;

%(2,1,6)Convolutional code rate 1/2
CONV_CODE_RATE = 1/2;
% Convolutionally encoding parameter
CONSTR_LEN = 7;% Constraint length
CODE_GEN = [171 133];% Polynomial
trellis = poly2trellis(CONSTR_LEN, CODE_GEN);

if Mod_Type==1 % QPSK symbols set
    A=[1+1j,1-1j,-1+1j,-1-1j];
    BIT_PER_SBL = 2;
elseif Mod_Type==2 % 16QAM symbols set
    A=[-3+3j,-1+3j,1+3j,3+3j,-3+1j,-1+1j,1+1j,3+1j,-3-1j,-1-1j,1-1j,3-1j,-3-3j,-1-3j,1-3j,3-3j];
    BIT_PER_SBL = 4;
elseif Mod_Type==3 % 64QAM symbols set
    A = [1+1j,1+3j,1+5j,1+7j,3+1j,3+3j,3+5j,3+7j,5+1j,5+3j,5+5j,5+7j,7+1j,7+3j,7+5j,7+7j,1-1j,1-3j,1-5j,1-7j,3-1j,3-3j,3-5j,3-7j,5-1j,5-3j,5-5j,5-7j,7-1j,7-3j,7-5j,7-7j,-1+1j,-1+3j,-1+5j,-1+7j,-3+1j,-3+3j,-3+5j,-3+7j,-5+1j,-5+3j,-5+5j,-5+7j,-7+1j,-7+3j,-7+5j,-7+7j,-1-1j,-1-3j,-1-5j,-1-7j,-3-1j,-3-3j,-3-5j,-3-7j,-5-1j,-5-3j,-5-5j,-5-7j,-7-1j,-7-3j,-7-5j,-7-7j];
    BIT_PER_SBL = 6;
end

% every users' symbols
evryUsrBts_Aftr_RS = nd*BIT_PER_SBL*CONV_CODE_RATE;

% No. of data bits = Massage bits + Group ID bits + User ID bits 
nEvryUsrDataBits = evryUsrBts_Aftr_RS * RS_CODE_RATE;
nEvryUsrIDBits= length(de2bi(NUM_ALL_USERS));  % length(de2bi(NUM_ALL_USERS))

% No. of symbols each user transmit in all time slots
evryUsrSbl = evryUsrBts_Aftr_RS/(BIT_PER_SBL*CONV_CODE_RATE);

% paramaters that are known in BS (for activity detectioZn and data recovery)
par.nd =nd ;
par.mRS = mRS;
par.nRS = nRS;
par.kRS = kRS;
par.trellis = trellis;
par.nEvryUsrDataBits = nEvryUsrDataBits;
par.nEvryUsrIDBits = nEvryUsrIDBits;
par.SNR = SNR;
par.nOfTfcb = nOfTfcb;
par.NUM_ALL_USERS = NUM_ALL_USERS;
par.Mod_Type=Mod_Type;
par.A=A;


NUM_ACT_USERS = 60;

drawPointMtxROW = 1; 
for NUM_PILOT =  [NUM_ALL_USERS]   % No. of pilot [64 200 500 5000 50000 NUM_ALL_USERS]  [64 200 5000 NUM_ALL_USERS]  
    drawPointMtxColumn = 2;
    drawPointMtxROW = drawPointMtxROW + 2;
    drawPointMtx1(drawPointMtxROW,1) = NUM_PILOT;
    drawPointMtx2(drawPointMtxROW,1) = NUM_PILOT;
    drawPointMtx3(drawPointMtxROW,1) = NUM_PILOT;     
    for NUM_RX = 4:4:64 % 5:5:100
        par.NUM_RX = NUM_RX ;
        drawPointMtx1(1,drawPointMtxColumn) = NUM_RX;
        drawPointMtx2(1,drawPointMtxColumn) = NUM_RX;
        drawPointMtx3(1,drawPointMtxColumn) = NUM_RX;        
        NUM_ITR = 100;
        for itr = 1:NUM_ITR
            % Pilot code-word pool (non-orthometric)
            cwplt = xi*(1/mp)*(randn(mp,NUM_PILOT)+1j*randn(mp,NUM_PILOT));
            % Symbol code-word matrix (Corresponds to pilot)
            cwsbl = 2*randi([0,1],md,NUM_PILOT)-1;  
            if NUM_PILOT==NUM_ALL_USERS
                acUsrIndx = sort(randperm(NUM_ALL_USERS,NUM_ACT_USERS));
                acUsr_pilot = zeros(NUM_ALL_USERS, 1);
                acUsr_pilot(acUsrIndx,:)=acUsrIndx;                
            else 
                %=================Step 1: random send pilot & ID &data  =================
                % Indices of selected pilot (may be repeated)
                pilotIndx = unidrnd(NUM_PILOT,[1,NUM_ACT_USERS]);
                % Indices of selected pilot (not repeated)        
                pilot_choose =sort(unique(pilotIndx));
                % No of selected pilots             
                NUM_SELECT_PILOT=length(pilot_choose);
                % Indices of active users
                acUsrIndx = sort(randperm(NUM_ALL_USERS,NUM_ACT_USERS));
                %the active user index for each pilot
                pilot_acUsr = cell(NUM_PILOT,1);
                for i = 1:NUM_ACT_USERS
                     pilot_acUsr{pilotIndx(i)}= [pilot_acUsr{pilotIndx(i)};acUsrIndx(i)];
                end    
                %the pilot index for each active user
                acUsr_pilot = zeros(NUM_ALL_USERS, 1);
                acUsr_pilot(acUsrIndx,:)=pilotIndx;
            end    
            
            % User ID 
            usrID = zeros(NUM_ALL_USERS,nEvryUsrIDBits);
            usrID(acUsrIndx,:) = de2bi(acUsrIndx,nEvryUsrIDBits,'left-msb');
            % Generating data
            acData = zeros(NUM_ALL_USERS, nEvryUsrDataBits);
            acData(acUsrIndx,:) = randi([0,1],NUM_ACT_USERS,nEvryUsrDataBits);
            % Intigrating User ID into data
            acData(acUsrIndx,nEvryUsrDataBits-nEvryUsrIDBits+1:nEvryUsrDataBits) = usrID(acUsrIndx,:);            

            clear usrID
            
            % RS Coding
            RSmatrix = reshape(acData(acUsrIndx,:)',mRS,nEvryUsrDataBits*NUM_ACT_USERS/mRS)';
            RStrSbl = bi2de(RSmatrix,'left-msb');
            hEnc = comm.RSEncoder(nRS,kRS);% RS parameter
            RS_encoded_Sbl = step(hEnc, RStrSbl);
            RS_encoded_Bits = de2bi(RS_encoded_Sbl',mRS,'left-msb');
            acdata_Aftr_RS = zeros(NUM_ALL_USERS,evryUsrBts_Aftr_RS);
            acdata_Aftr_RS(acUsrIndx,:) = reshape(RS_encoded_Bits',evryUsrBts_Aftr_RS,NUM_ACT_USERS)';          
            
           
            clear RSmatrix
            clear RStrSbl
            clear hEnc
            clear RS_encoded_Sbl
            clear RS_encoded_Bits
            
            % Convolutionally encoding
            for q = 1:NUM_ACT_USERS
                codeData(acUsrIndx(q),:) = convenc(acdata_Aftr_RS(acUsrIndx(q),:)', trellis);% Coding data
            end
            matrix = reshape(codeData(acUsrIndx,:)',4,NUM_ACT_USERS*evryUsrBts_Aftr_RS/(BIT_PER_SBL*CONV_CODE_RATE));
            % Interleaving
            bits_Aftr_Itrlv = matintrlv(matrix,2,2)';
     
            clear acdata_Aftr_RS  
            clear codeData  
            clear matrix 
            
            %====================Modulation====================           
            if Mod_Type==1 % QPSK
                %输入x为0 1 2 3  输出y为1+1i 1-1i -1+1i -1-1i
                bits_Aftr_Itrlv = reshape(bits_Aftr_Itrlv',2,[])';
                dec = bi2de(bits_Aftr_Itrlv,'left-msb');
                sbl_Aftr_Mod = zeros(NUM_ALL_USERS,evryUsrSbl);
                sbl_Aftr_Mod(acUsrIndx,:) =reshape(A(dec(:,:)+1),[],NUM_ACT_USERS).';                
            elseif Mod_Type==2 % 16QAM 
                bits_Aftr_Itrlv = reshape(bits_Aftr_Itrlv',4,[])';
                dec = bi2de(bits_Aftr_Itrlv,'left-msb'); % Binary to decimal conversion
                sbl_Aftr_Mod = zeros(NUM_ALL_USERS,evryUsrSbl);
                sbl_Aftr_Mod(acUsrIndx,:) = reshape(qammod(dec,16),evryUsrSbl,NUM_ACT_USERS)';                   
            elseif Mod_Type==3 % 64QAM
                dec = bi2de(bits_Aftr_Itrlv,'left-msb'); % Binary to decimal conversion
                sbl_Aftr_Mod = zeros(NUM_ALL_USERS,evryUsrSbl);
                sbl_Aftr_Mod(acUsrIndx,:) = reshape(qammod(dec,64),evryUsrSbl,NUM_ACT_USERS)';                
            end            
                    
            clear bits_Aftr_Itrlv
            clear dec
            
            %=================Mapping to TFCB=================
            % possible active users' free space loss coefficient
            Lc = repmat(rand(NUM_ALL_USERS,1),1,NUM_RX);      
            
            % All users' M RXs channels in all TFCBs
            H = zeros(NUM_ALL_USERS,NUM_RX);
            %Rayleigh Channel
            H(acUsrIndx,:) = (1/sqrt(2))*(randn(NUM_ACT_USERS,NUM_RX)+1j*randn(NUM_ACT_USERS,NUM_RX));
            H = H .* Lc;
            
            Xsbl = repelem(sbl_Aftr_Mod,1,NUM_RX);
            Hsbl_est = mat2cell(H,[size(H,1)],NUM_RX);
            Hsbl_est = repelem(Hsbl_est,1,nd);
            Hsbl_est = cell2mat(Hsbl_est);
            HX = Hsbl_est .* Xsbl;
            
            clear sbl_Aftr_Mod
            clear Lc
%             clear Hsbl_est
            clear Xsbl
            
            
            if NUM_PILOT==NUM_ALL_USERS
                Y_pilot = cwplt*H;% Y of pilot
                Y_data = cwsbl*HX;% Y of data
            else
                % the channel for each pilot(superposed when collision occurs)
                H_pilot = zeros(NUM_PILOT,NUM_RX);
                HX_pilot = zeros(NUM_PILOT,nd*NUM_RX);
                for i = 1:NUM_ACT_USERS
                    H_pilot(pilotIndx(i),:) = H_pilot(pilotIndx(i),:) + H(acUsrIndx(i),:); 
                    HX_pilot(pilotIndx(i),:) = HX_pilot(pilotIndx(i),:) + HX(acUsrIndx(i),:); 
                end

                Y_pilot = cwplt*H_pilot;% Y of pilot
                Y_data = cwsbl*HX_pilot;% Y of data
            end
            
            Y_pilot_rx = awgn(Y_pilot,SNR,'measured');
            Y_data_rx = awgn(Y_data,SNR,'measured');                                
            

            
            %=================Step 2: solve the selected pilot and data =================
            if NUM_PILOT == 64
                [A_est,B_est, D_est,H_est]= LTE_FOUR_STEP_MIMO(Y_pilot_rx,cwplt,HX,cwsbl,acData,acUsrIndx,pilot_acUsr,acUsr_pilot,par,pilotIndx,pilot_choose);  
                drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) + length(intersect(B_est,pilot_choose)) / length(pilot_choose);                 
                drawPointMtx3(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx3(drawPointMtxROW,drawPointMtxColumn) + length(D_est)/ NUM_ACT_USERS;
            elseif NUM_PILOT==NUM_ALL_USERS
%                 tic
                [B_est, H_est] = OMP_blind(Y_pilot_rx,cwplt,SNR);  
                [A_est,D_est ] = Data_LS(B_est,H_est,Y_data_rx, cwsbl,NUM_PILOT,acUsrIndx,acData,par);
                drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) + length(intersect(B_est,acUsrIndx)) / NUM_ACT_USERS; 
                drawPointMtx3(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx3(drawPointMtxROW,drawPointMtxColumn) + length(D_est)/ NUM_ACT_USERS;   

                [B_est, H_est] = Joint_OMP_blind(Y_pilot_rx,cwplt,Y_data_rx,cwsbl,SNR);
                [A_est,D_est ] = CS_Data_LS_innercycle(B_est,H_est,Y_data_rx, cwsbl,NUM_PILOT,acUsrIndx,acData(acUsrIndx,:),cwplt,Y_pilot_rx,par);
%                 toc
                drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx1(drawPointMtxROW+1,drawPointMtxColumn) + length(intersect(B_est,acUsrIndx)) / NUM_ACT_USERS; 
                drawPointMtx3(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx3(drawPointMtxROW+1,drawPointMtxColumn) + length(D_est)/ NUM_ACT_USERS;   
%                 drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) + toc;
            else       
                [B_est, H_est] = OMP_blind(Y_pilot_rx,cwplt,SNR);    
                drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx1(drawPointMtxROW,drawPointMtxColumn) + length(intersect(B_est,pilot_choose)) / length(pilot_choose);
                [A_est,D_est ] = Data_LS(B_est,H_est,Y_data_rx, cwsbl,NUM_PILOT,acUsrIndx,acData,par);
                drawPointMtx3(drawPointMtxROW,drawPointMtxColumn) = drawPointMtx3(drawPointMtxROW,drawPointMtxColumn) + length(D_est)/ NUM_ACT_USERS;                 
                [B_est, H_est] = Joint_OMP_blind(Y_pilot_rx,cwplt,Y_data_rx,cwsbl,SNR);
                drawPointMtx1(drawPointMtxROW+1,drawPointMtxColumn) = drawPointMtx1(drawPointMtxROW+1,drawPointMtxColumn) + length(intersect(B_est,pilot_choose)) / length(pilot_choose);                 
                [A_est,D_est ] = Data_LS_innercycle(B_est,H_est,Y_data_rx, cwsbl,NUM_PILOT,acUsrIndx,acData,acUsr_pilot,cwplt,Y_pilot_rx,par);
                drawPointMtx3(drawPointMtxROW+1,drawPointMtxColumn) = drawPointMtx3(drawPointMtxROW+1,drawPointMtxColumn) + length(D_est)/ NUM_ACT_USERS;                             
            end           
            fprintf('NUM_PILOT = %d, ACT = %d, itr = %d \n',NUM_PILOT,NUM_ACT_USERS,itr);
                
        end
        drawPointMtxColumn = drawPointMtxColumn + 1;
    end
end

drawPointMtx1(2:NUM_ROW,2:end)  = drawPointMtx1(2:NUM_ROW,2:end) ./ NUM_ITR;
drawPointMtx3(2:NUM_ROW,2:end)  = drawPointMtx3(2:NUM_ROW,2:end) ./ NUM_ITR;

% save('figure_diff_acc_1000000user_l72','drawPointMtx1','drawPointMtx2','drawPointMtx3');
