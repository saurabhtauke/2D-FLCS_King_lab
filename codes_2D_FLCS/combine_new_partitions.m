% To combine 2 paritions into one.
% add the partitions and name them new_part1 and new_part2


np = cell(size(new_part1));

for i = 1:32
    np{i} = struct();
    np{i}.sh1D_dec1 = new_part1{i}.sh1D_dec1 + new_part2{i}.sh1D_dec1;
    np{i}.sh1D_dec2 = new_part1{i}.sh1D_dec2 + new_part2{i}.sh1D_dec2;
    np{i}.sh2D = new_part1{i}.sh2D + new_part2{i}.sh2D;
end