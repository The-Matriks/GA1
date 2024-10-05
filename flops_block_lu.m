function flops = flops_block_lu(n, b)

    num_blocks = (n / b)^2;
    flops_per_block = (2/3) * b^3;
    flops = num_blocks * flops_per_block;
end