function recieved_seq = get_eq_output_wiener_solution_known_chanel(H,Xn,noise_var)

Zn = cell(size(Xn));

%calculating the Zn vectors------------------------------------------------
for i = 1 : length(Zn)
    
    Wn = sqrt(noise_var)*randn(size(H,1),1) + j*(sqrt(noise_var)*randn(size(H,1),1));
    Zn{i} = H*Xn{i}' + Wn;

end
%--------------------------------------------------------------------------

%calculating the optimal(in the MMSE sense) coeficients of the linear
%equalizer -----------------------------------------------------------
d = zeros(size(H,2),1);
d((size(H,2)+1)/2) = 1;

f_star = (inv(H*H' + noise_var.*eye(size(H,1))))*(H*d);
%---------------------------------------------------------------------

%Calculating the estimated symbols at the recieving end (output of the EQ)
%This output(recieving_seq) is gonna be the input of the decision device 
for i = 1 : length(Zn)
    
    recieved_seq(i) = f_star'*Zn{i};

end   

end

