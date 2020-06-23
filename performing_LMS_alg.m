function [recieved_seq,error,filter_coefs,Zn] = performing_LMS_alg(H , Xn , noise_var , known_symbols , step)



Zn = cell(size(Xn));

%calculating the Zn vectors------------------------------------------------
for i = 1 : length(Zn)
    
    Wn = sqrt(noise_var)*randn(size(H,1),1) + j*(sqrt(noise_var)*randn(size(H,1),1));
    Zn{i} = H*Xn{i}' + Wn;

end
%--------------------------------------------------------------------------

%LMS ALGORITHM-----------------------------------------------------------

Xn_hat = [];
err = [];

%filter coefs for each iteration of the algorithm
Fn = cell(size(known_symbols));

%initialization
Fn{1} = [ zeros(( size(H,1)-1)/2 ,1) ; 1 ; zeros(( size(H,1)-1)/2 ,1) ];

%Run the algorithm for T iterations(number of learning symbols) 
%TRAINING PHASE------------------------------------------------------
for i = 1 : length(known_symbols)
    
    Xn_hat = [ Xn_hat Fn{i}'*Zn{i} ];
    
    
    e = known_symbols(i) - Xn_hat(i);
    
    err = [ err e ];
    
    Fn{i+1} = Fn{i} + step*(Zn{i}.*conj(e));    
    
end
    
%This block of code is for the 4th question (obtaining the mse of the training)
recieved_seq = Xn_hat;
error = err;

%returning the filter coefs of the last iteration of the algorithm
filter_coefs = Fn{length(Fn)};

end
    
    
    