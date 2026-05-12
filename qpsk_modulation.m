%qpsk_modulation.m
classdef qpsk_modulation
    methods(Static)

        function [a, N_mod_symbols] = modulate(bits, params)

            bits_per_symbol = params.bits_per_symbol;

            bit_pairs = reshape(bits, bits_per_symbol, []).';
            N_mod_symbols = length(bits)/bits_per_symbol;
            a = zeros(1, N_mod_symbols);

            idx = 1;
            for i = 1:N_mod_symbols
                b1 = bit_pairs(i,1);
                b2 = bit_pairs(i,2);
                
                %Constellation pairs - Gray coding
                if (b1 == 0 && b2 == 0)
                    a(idx) = 1 + 1j;
                elseif (b1 == 1 && b2 == 0)
                    a(idx) = -1 + 1j;
                elseif (b1 == 1 && b2 == 1)
                    a(idx) = -1 - 1j;
                elseif (b1 == 0 && b2 == 1)
                    a(idx) = 1 - 1j;
                end

                idx = idx + 1;
            end
        end


        function bits_decoded = demodulate(data_decoded_symbols)

            bits_decoded = zeros(1, 2*length(data_decoded_symbols));

            idx = 1;
            for i = 1:length(data_decoded_symbols)
                sym = data_decoded_symbols(i);

                I = real(sym);
                Q = imag(sym);

                if (I >= 0 && Q >= 0)
                    bits_decoded(idx:idx+1) = [0 0];

                elseif (I < 0 && Q >= 0)
                    bits_decoded(idx:idx+1) = [1 0];

                elseif (I < 0 && Q < 0)
                    bits_decoded(idx:idx+1) = [1 1];

                else
                    bits_decoded(idx:idx+1) = [0 1];
                end

                idx = idx + 2;
            end
        end

    end
end