function crc_bits = compute_crc16(data_bits)
% Compute CRC-16-CCITT over a bit vector
% data_bits: row vector of bits
% Returns 16 CRC bits

gen = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1]; % x^16+x^12+x^5+1
msg = [double(data_bits(:).')  zeros(1,16)];

for i = 1:length(data_bits)
    if msg(i) == 1
        msg(i:i+16) = xor(msg(i:i+16), gen);
    end
end

crc_bits = msg(end-15:end);
end