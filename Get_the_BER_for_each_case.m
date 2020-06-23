clear;

BER = zeros(length(1:35),1);
SNRdb = [];

x = input(['To obtain the BER curve for the ideal chanel press 1',...
       '\nTo obtain the BER curve for the PROAKIS A chanel press 2',...
       '\nTo obtain the BER curve for the PROAKIS B chanel press 3\n']);
   
if (x == 2) || (x == 3)        
    y = input(['To see the BER curve that is achieved with the wiener solution(known chanel) press 1',...
           '\nTo see the BER curve that is achieved with the LMS algorithm(unknown chanel) press 2\n']);
    if y == 2
        z = input(['To see the BER curve with the FREEZED coeficients press 1',...
            '\nTo see the BER curve in DECISION DIRECTED MODE press 2\n']);
    end    
end   
   
%number of known symbols (learning symbols)
T = 800;
%step of the LMS algorithm
step = 0.01;

for SNR = 1 : 35

SNRdb = [ SNRdb SNR ];    
k = 2;

bit_stream = randsrc(1,100000,[0 1; (1/2) (1/2)]);

%converting bit_stream to the corresponding decimals
if mod(length(bit_stream),k) == 0
    xsym = reshape(bit_stream,k,length(bit_stream)/k).';
else
    %pad a 0 at the msb
    bit_stream = [ 0 bit_stream ];
    xsym = reshape(bit_stream,k,length(bit_stream)/k).';
end

%creating the symbols that correspnds to the 4-PSK constellation diagram
symbols = [];
symbols = [ symbols exp(j*0) ];
symbols = [ symbols exp(j*(pi/2)) ];
symbols = [ symbols exp(j*(pi)) ];
symbols = [ symbols exp(j*((3*pi)/2)) ];

%creating the mapping: Bit_stream -> Symbols (gray encode) 
symbol_seq = zeros([size(xsym,1),1]);
for i = 1 : size(xsym,1)
    
        if xsym(i,:) == [0 0];
            symbol_seq(i) = symbols(1);
            
        elseif xsym(i,:) == [0 1];
            symbol_seq(i) = symbols(2);
            
        elseif xsym(i,:) == [1 1];
            symbol_seq(i) = symbols(3);    
            
        elseif xsym(i,:) == [1 0];
            symbol_seq(i) = symbols(4);
            
            
        end    
    
end    

%noise---------------------
noise_var = 1/(4*(10^(SNR/10)));
noise = sqrt(noise_var)*randn(size(symbol_seq)) + j*(sqrt(noise_var)*randn(size(symbol_seq))); 
%----------------------------------------------

%ideal chanel---------
%not using an equalizer
if x == 1
    l = -1:0.1:1;
    h1 = dirac(l);
    h1(11) = 1;
    recieved_seq = conv(symbol_seq,h1,'same') + noise;
    title_str = 'ideal chanel';
end    
%---------------------
%}



%proakis A chanell --------------------------------------------------------
if x == 2
    
    h2 = [ 0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0 0.21 0.03 0.07 ];
    %recieved_seq = conv(symbol_seq,h2,'same') + noise;

    [H,Xn] = get_chanel_matrix_and_Xn_vectors(31,h2,symbol_seq);
    
    %WIENER solution known chanel
    if y == 1
        recieved_seq = get_eq_output_wiener_solution_known_chanel(H,Xn,noise_var);
        title_str = 'WIENER solution known chanel (PROAKIS A)';
    end    
        
    if y == 2
        
        %ERWTHMA 5 A)
        %FREEZED coefs------------------------------------------------------------
        if z == 1    
            known_symbols = symbol_seq(1:T);
            [recieved_seq,error,filter_coefs,Zn] = performing_LMS_alg(H , Xn , noise_var , known_symbols , step);
            %obtaining the mse for the case of the freezed filter coefs

            for i = T + 1 : length(Zn)
    
                Xn_hat = filter_coefs'*Zn{i};
    
                error = [ error ( symbol_seq(i) - Xn_hat ) ];
  
            end 
        
            title_str = 'FREEZED coefs LMS algorithm (PROAKIS A)';
        end   
       
     %-------------------------------------------------------------------------

        if z == 2 
            %ERWTHMA 5 B)
            %DECISION DIRECTED MODE----------------------------------------------------
            [recieved_seq,error,Zn] = performing_LMS_alg_decision_directed(H , Xn , noise_var , symbol_seq , T , symbols , step);
            title_str = 'DECISION DIRECTED MODE LMS algorithm (PROAKIS A)';
            %-------------------------------------------------------------------------
        end
    end    
    %------------------------------------end of chanel A-----------------------
end


%PROAKIS B chanel----------------------------------------------------------
if x == 3
    
    h3 = [0.407 0.815 0.407];
    %recieved_seq = conv(symbol_seq,h2,'same') + noise;

    [H,Xn] = get_chanel_matrix_and_Xn_vectors(31,h3,symbol_seq);
    
    %WIENER solution known chanel
    if y == 1
        recieved_seq = get_eq_output_wiener_solution_known_chanel(H,Xn,noise_var);
        title_str = 'WIENER solution known chanel (PROAKIS B)';
    end
    
    if y == 2 
    
        %ERWTHMA 5 A)
        %FREEZED coefs------------------------------------------------------------
        if z == 1    
            known_symbols = symbol_seq(1:T);
            [recieved_seq,error,filter_coefs,Zn] = performing_LMS_alg(H , Xn , noise_var , known_symbols , step);
            %obtaining the mse for the case of the freezed filter coefs

            for i = T + 1 : length(Zn)
    
                Xn_hat = filter_coefs'*Zn{i};
    
                error = [ error ( symbol_seq(i) - Xn_hat ) ];
  
            end 
            
            title_str = 'FREEZED coefs LMS algorithm (PROAKIS B)';
        end   
        %-------------------------------------------------------------------------

        if z == 2 
            %ERWTHMA 5 B)
            %DECISION DIRECTED MODE----------------------------------------------------
            [recieved_seq,error,Zn] = performing_LMS_alg_decision_directed(H , Xn , noise_var , symbol_seq , T , symbols , step);
            title_str = 'DECISION DIRECTED MODE LMS algorithm (PROAKIS B)';
            %-------------------------------------------------------------------------
        end
    end    
    %------------------------------------end of chanel B-----------------------
end


%deciding which symbol (or signal) with respect to the l2 norm of the distance
%MAXIMUM LIKELIHOOD
dist = cell(size(symbol_seq));
symbol_seq_hat = zeros(size(symbol_seq));

for i = 1 : length(recieved_seq)
    for k = 1 : 4
        dist{i} = [ dist{i}  norm(recieved_seq(i) - symbols(k))^2 ];
    end 
    [m,ind] = min(dist{i});
    symbol_seq_hat(i) = symbols(ind);
end



%creating the inverse mapping : symbols -> bit stream

bit_stream_hat = [];

for i = 1 : length(symbol_seq_hat)
    
    if angle(symbol_seq_hat(i)) == 0
        temp = [0 0];
        bit_stream_hat = [ bit_stream_hat temp ];
        
    elseif  angle(symbol_seq_hat(i)) == (pi/2)
        temp = [0 1];
        bit_stream_hat = [ bit_stream_hat temp ];
        
    elseif  angle(symbol_seq_hat(i)) == pi
        temp = [1 1];
        bit_stream_hat = [ bit_stream_hat temp ]; 
        
    elseif  angle(symbol_seq_hat(i)) ~= (-(pi/2))
        temp = [1 0];
        bit_stream_hat = [ bit_stream_hat temp ]; 
       
    end
    
end    


BER(SNR) = length(find(bitxor(bit_stream_hat , bit_stream))) / length(bit_stream);

end

figure,semilogy(SNRdb,BER);
title(title_str);
xlabel('SNR(db)');
ylabel('BER');