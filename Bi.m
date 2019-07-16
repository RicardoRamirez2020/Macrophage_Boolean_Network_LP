function [B,mean] = Bi(itter,D)



for m = 1:D
    K = 2^(D-(m-1));
    [Bin,C] = kmeans(itter,K,'Distance','cityblock','EmptyAction', 'singleton');
    %if (m ~=D)
        
    for n = 1:length(Bin)
        itter(n) = C(Bin(n));
    end
    %end
    
end

mean = C; 

bool = false;
if (C(1)>C(2))
bool = true;     
end

for m = 1:length(Bin) 
    if(Bin(m) == 1)
    b(m) = bool; 
    
    else
    b(m) = not(bool);    
    end
end

for m = 1:length(Bin)
   if(b(m)==true)
       B(m) = 1;
   else
       B(m) = 0;
   end    
       
end

end