function [B, H] = Joint_OMP_blind(Y_p, P,Y_d, S,SNR)
a=1.2;%阈值调整参数
max_itr = max(size(P));
t1 = 1;
B = [];
Y_r_p = Y_p;
[~,N_p]=size(P);
H = zeros(N_p,size(Y_p,2));
Y_r_d = Y_d;
[~,N_d]=size(S);
HX = zeros(N_d,size(Y_d,2));
norm_save(t1) = norm(Y_r_p,'fro')+norm(Y_r_d,'fro');
while 1
    t1 = t1 + 1;
    B_last = B;    
    [~,k] = max(sum(abs(P'*Y_r_p).^2,2)+sum(abs(S'*Y_r_d).^2,2));
    B = union(k,B);
    H(B,:) = P(:,B) \ Y_p;
    Y_r_p = Y_p - P(:,B)*H(B,:);
    HX(B,:) = S(:,B) \ Y_d;
    Y_r_d = Y_d - S(:,B)*HX(B,:);      
    norm_save(t1) = norm(Y_r_p,'fro')+norm(Y_r_d,'fro');
    
    if 20*log10( norm(P(:,B)*H(B,:),'fro')/norm(Y_r_p,'fro'))+20*log10( norm(S(:,B)*HX(B,:),'fro')/norm(Y_r_d,'fro')) > 2*a*SNR
        break;
    end    
       
    if norm_save(t1)/norm_save(t1-1) >= 1 
        B = B_last';
        H(B,:) = pinv(P(:,B))*Y_p;
        break;
    end
    
    if t1 >max_itr
        break;
    end
    
end
