clear;

T = input('Give the number of learning symbols : ');

step = input('Give the step of the LMS algorithm(learning rate) : ');

SNR = input('Give the desiered SNR : ');

x = input(['To obtain the MSE curve for the PROAKIS A chanel press 1',...
       '\nTo obtain the MSE curve for the PROAKIS B chanel press 2\n']);

z = input(['For the case of the FREEZED coefs press 1',...
       '\nFor the case of the DECISION DIRECTED MODE press 2\n']);   
   
k = 2;
%SNR = 20;

%number of known symbols (learning symbols)
%T = 800;
%step of the LMS algorithm
%step = 0.01;

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


%proakis A chanell --------------------------------------------------------
if x == 1
    
    h2 = [ 0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0 0.21 0.03 0.07 ];
    %recieved_seq = conv(symbol_seq,h2,'same') + noise;

    [H,Xn] = get_chanel_matrix_and_Xn_vectors(31,h2,symbol_seq);  
        
    
        
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
       
    %------------------------------------end of chanel A-----------------------
end

%PROAKIS B chanel----------------------------------------------------------
if x == 2
    
    h3 = [0.407 0.815 0.407];
    %recieved_seq = conv(symbol_seq,h2,'same') + noise;

    [H,Xn] = get_chanel_matrix_and_Xn_vectors(31,h3,symbol_seq);
   
    
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
    %------------------------------------end of chanel B-----------------------
end


%ploting MSE vs number of iteration (symbol at time n)
figure, plot(error.*conj(error));
title(title_str);
xlabel('Symbol');
ylabel('MSE');