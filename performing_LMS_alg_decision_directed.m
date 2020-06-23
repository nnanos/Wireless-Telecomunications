function [recieved_seq,error,Zn] = performing_LMS_alg_decision_directed(H , Xn , noise_var , symbol_seq , T , symbols , step)

known_symbols = symbol_seq(1:T);

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
Fn = cell((size(symbol_seq,1)+1),1);

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


%decision directed mode 
for i = T + 1 : length(symbol_seq)
    
    Xn_hat = [ Xn_hat Fn{i}'*Zn{i} ];
    
    %we pass the Xn_hat by the decision device and the output is going to
    %be the desiered filter output
    %deciding which symbol (or signal) with respect to the l2 norm of the distance
    %MAXIMUM LIKELIHOOD---------------------------------------------------
    dist = [];
        for k = 1 : 4
            dist = [ dist  norm(Xn_hat(i) - symbols(k)) ];
        end 
        [m,ind] = min(dist);
        desiered_output = symbols(ind);

    %-------------------------------------------------------------------
    
    e = desiered_output - Xn_hat(i);
    
    %err = [ err e ];
    
    Fn{i+1} = Fn{i} + step*(Zn{i}.*conj(e));

end


%This block of code is for the 4th question (obtaining the mse of the training)
recieved_seq = Xn_hat;

err1 = [];

for i =  T + 1 : length(symbol_seq)
     err1 = [ err1 symbol_seq(i) - recieved_seq(i)];
end     

error = [ err err1 ];


end