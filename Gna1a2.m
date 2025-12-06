function [a1, a2] = Gna1a2(Np1,r0)

a1 = randperm(Np1,Np1);

for i = 1 : 99999999
    pos = (a1 == r0);
    if sum(pos) == 0
        break;
    else 
        a1(pos) = floor(rand(1, sum(pos)) * Np1) + 1;
    end
    if i > 1000 % this has never happened so far
        error('Can not genrate a1 in 1000 iterations');
    end
end

a2 = floor(rand(1, Np1) * Np1) + 1;

for i = 1 : 99999999
    pos = ((a2 == a1) | (a2 == r0));
    if sum(pos)==0
        break;
    else 
        a2(pos) = floor(rand(1, sum(pos)) * Np1) + 1;
    end
    if i > 1000 
        error('Can not genrate a2 in 1000 iterations');
    end
end
