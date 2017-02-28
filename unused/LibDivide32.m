classdef LibDivide32
% Fast integer division for unsigned integers up to (2^32 - 1)
    properties
        D;          % unsigned divisor
        multiplier;
        pre_shift;
        post_shift;
        full_shift;
        increment;
    end
    methods
        function result = LibDivide32(D, varargin)
            if length(varargin) < 1
                num_bits = 32; % our divident and divisor is 32 bits
            else
                num_bits = varargin{1};
            end
            UINT_BITS = 32;

            % The numerator must fit in a uint
            assert(num_bits > 0 && num_bits <= UINT_BITS);
    
            % Compute ceil(log_2 D) when D is not a power of 2
            % or return (k + 1) when D = 2^k
            ceil_log_2_D = uint64(0);
            tmp = D;
            while (tmp > 0)
                ceil_log_2_D = ceil_log_2_D + 1;                
                tmp = bitshift(tmp, -1); % >> 1
            end

            D = uint64(D);
            result.D = D;

            % D must be larger than zero and not a power of 2
            is_not_power_of_two = bitand(D, (D - 1));
            if ~is_not_power_of_two
                result.multiplier = uint64(1);
                result.increment = uint64(0);
                result.pre_shift = int32(0);
                result.post_shift = int32(0);
                result.full_shift = int32(ceil_log_2_D - 1);
                return
            end

            % The extra shift implicit in the difference between UINT_BITS and num_bits
            extra_shift = UINT_BITS - num_bits;
    
            % The initial power of 2 is one less than the first one that can possibly work
            initial_power_of_2 = bitshift(uint64(1), UINT_BITS - 1); % << (UINT_BITS - 1)

            % The remainder and quotient of our power of 2 divided by d
            quotient = idivide(initial_power_of_2, D, 'floor');
            remainder = mod(initial_power_of_2, D);
            
            % The magic info for the variant "round down" algorithm
            down_multiplier = uint64(0);
            down_exponent = 0; % small integer
            has_magic_down = 0; % boolean
    
            % Begin a loop that increments the exponent, until we find a power of 2 that works.

            % for (exponent = 0; ; exponent++)
            exponent = uint64(0);
            while 1
                % Quotient and remainder is from previous exponent; compute it for this exponent.
                if remainder >= D - remainder
                    % Doubling remainder will wrap around D
                    quotient = quotient * 2 + 1;
                    remainder = remainder * 2 - D;
                else
                    % Remainder will not wrap
                    quotient = quotient * 2;
                    remainder = remainder * 2;
                end
                
                % We're done if this exponent works for the round_up algorithm.
                % Note that exponent may be larger than the maximum shift supported,
                % so the check for >= ceil_log_2_D is critical.
                if (exponent + extra_shift) >= ceil_log_2_D
                    break
                end
                if (D - remainder) <= bitshift(uint64(1), exponent + extra_shift)
                    break
                end
                               
                % Set magic_down if we have not set it yet and this exponent works for the round_down algorithm
                if ~has_magic_down && remainder <= bitshift(uint64(1), (exponent + extra_shift))
                    has_magic_down = 1;
                    down_multiplier = quotient;
                    down_exponent = exponent;
                end
                exponent = exponent + 1;
            end % while 1 // for (exponent = 0 ; ; exponent ++)

            if (exponent < ceil_log_2_D)
                % magic_up is efficient
                result.multiplier = quotient + 1;
                result.pre_shift = 0;
                result.post_shift = exponent;
                result.increment = 0;
            elseif bitand(D, 1) == 1
                % Odd divisor, so use magic_down, which must have been set
                assert(has_magic_down == 1);
                result.multiplier = down_multiplier;
                result.pre_shift = 0;
                result.post_shift = down_exponent;
                result.increment = 1;
            else
                % Even divisor, so use a prefix-shifted dividend
                pre_shift = 0;
                shifted_D = D;
                while bitand(shifted_D, 1) == 0
                    shifted_D = bitshift(shifted_D, -1); % >>= 1;
                    pre_shift = pre_shift + 1;
                end
                result1 = LibDivide32(shifted_D, num_bits - pre_shift);
                result.multiplier = result1.multiplier;
                result.pre_shift = result1.pre_shift;
                result.post_shift = result1.post_shift;
                result.increment = result1.increment;               
                assert(result.increment == 0 && result.pre_shift == 0); % expect no increment or pre_shift in this path
                result.pre_shift = pre_shift;
            end
            result.pre_shift = int32(result.pre_shift);
            result.post_shift = int32(result.post_shift);
            result.full_shift = int32(result.post_shift + UINT_BITS);
        end
        function res = div(ld, a)
            res = uint64(a);
            if ld.pre_shift > 0
                res = bitshift(res, -ld.pre_shift);
            end
            res = res * ld.multiplier;
            if ld.increment
                res = res + ld.multiplier;
            end
            res = bitshift(res, -ld.full_shift);
            res = cast(res, class(a));
        end
        function res = mod(ld, a)
            res = uint64(a);
            cast(res - ld.div(res), class(a));
        end
    end
end
