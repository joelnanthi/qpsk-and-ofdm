%channel_coding.m
%convolutional coding
classdef channel_coding
    methods(Static)

        function [coded_bits, N_coded_bits] = encode(bits, decoding_active)
            if (decoding_active == true)
                trellis = poly2trellis(6, [77 45]); %Encoder structure
                memory_length = 5; %5 zero bits

                bits_with_tail = [bits zeros(1, memory_length)]; %5 zeros at end
                coded_bits = convenc(bits_with_tail, trellis); %Convolutional encoder

                N_coded_bits = length(coded_bits); 

            elseif (decoding_active == false)

                coded_bits = bits;
                N_coded_bits = length(coded_bits);

            end
        end

        function data_decoded_bits = decode(bits_decoded, Nbits, decoding_active)
            if decoding_active == true
                trellis       = poly2trellis(6, [77 45]);
                memory_length = 5;
                tblen         = 30;

                length_data_bits          = Nbits;
                length_encoder_input_bits = length_data_bits + memory_length;
                length_coded_bits_needed  = 2 * length_encoder_input_bits;

                %% ── ADD THIS SAFETY PAD ──────────────────────────────────
                if length(bits_decoded) < length_coded_bits_needed
                    bits_decoded = [bits_decoded ...
                        zeros(1, length_coded_bits_needed - length(bits_decoded))];
                end
                %% ──────────────────────────────────────────────────────────

                bits_decoded      = bits_decoded(1:length_coded_bits_needed);
                data_decoded_bits = vitdec(bits_decoded,trellis,tblen,'term','hard');
                data_decoded_bits = data_decoded_bits(1:length_data_bits);
            else
                data_decoded_bits = bits_decoded(1:Nbits);
            end
        end

    end
end
