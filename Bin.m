Data = I4_LOG_9_abs;
Data  = permute(Data, [3 1 2]);

[ A, B, C] = size(Data);

MEAN  = size(B,C);

for x = 1:C

    for y = 1:B
        
        [Z, z ] = Bi(Data(:,y,x),2);
        MEAN(y,x) = sum(z)/2;
        
        
    end


end
MEAN(MEAN < -100000) = 0;

M_50_Total = zeros(50,27);
I4_50_9 = zeros(27,1,50);
for z = 1:50
    
    M_50_Total(z,:) = Mn((5*(z-1))+1,:);
    %I4_50_9(:,:,z) = I4(:,:,(2*(z-1))+1);

end

for z = 1:10
    
    M_50_9(z+40,:) = Mn((16*(z-1))+80,:);
   % I4_50_9(:,:,z+40) = I4(:,:,(16*(z-1))+80);

end








for x = 1:4
    for y = 1:27
        
        if (isnan(MEAN(x,y)) == 1)
        
            MEAN1(x,y) = .00001;
            
        
        else
            MEAN1(x,y) = MEAN(x,y);
        end
    end
end

    
for x = 1:50
   
    Msam(x,:,:) = M(1+((x-1)*5),:,:);
    
    
end