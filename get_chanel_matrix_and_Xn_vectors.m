function [H,Xn] = get_chanel_matrix_and_Xn_vectors(filter_order,chanel_impulse_response,symbol_seq)

%GET CHANEL CONVOLUTION MATRIX---------------------------------------------
K_L = (filter_order - 1)/2;
col_size = 2*K_L + length(chanel_impulse_response);

H = zeros(filter_order , col_size);
H(1,1:length(chanel_impulse_response)) = chanel_impulse_response;

for i = 1 : filter_order-1
    H((i+1),(i+1):length(chanel_impulse_response)+ i) = chanel_impulse_response;
end
%--------------------------------------------------------------------------

%GET THE Xn VECTORS----------------------------------------------
tmp = K_L + (length(chanel_impulse_response)-1)/2;
padded_array = [ zeros(1,tmp) symbol_seq' zeros(1,tmp) ];

Xn = cell(size(symbol_seq));
for n = 0 : length(symbol_seq)-1
    Xn{n+1} = padded_array((1:size(H,2))+ n);
end
%----------------------------------------------------------------

end

