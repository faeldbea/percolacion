function perco
clear all;
rep = 1
sigmalist = [0.55];

for j = 1:size(sigmalist,2)
    
sigma = sigmalist(1,j);

for i = 1:rep
u(i,j) = fporo(sigma);
end
end

save all2.txt u -ascii

end

function fporo = fporo(sigma) 

dporo = 100 ; % pore lengh in nm
rporo = 5 ; % pore radius in nm

Aporo = pi*2*rporo*dporo ; % internal area in nm^2

NN = round(sigma*Aporo)

rr = 0.5; % molecule radius size in nm
rr2 = rr^2 ; 

[xc,yc,zc] = cylinder(rporo, 100);
zc(2, :) = zc(2, :)*dporo;
M = mesh(xc,yc,zc);
set(M,'facealpha',0)

tdist = 1.5; % maximum tunneling ditance
tdist2 = tdist^2;


%% Create system

for j = 1:NN
check = 0;

while (check == 0) 
    check = 1;
    % random position
    theta = 2.0*pi*rand;
    z = dporo*rand;
    x = cos(theta)*rporo;
    y = sin(theta)*rporo;
   
    if(z - rr < 0) 
        check = 0;
        break;
    end
    
    for i = 1:j-1
        dist = (x-pos(i,1))^2 + (y-pos(i,2))^2 + (z-pos(i,3))^2;
            if(dist < rr2)
                check = 0;
                break;
            end
    end

   
end
   pos(j,1) = x;
   pos(j,2) = y;
   pos(j,3) = z;
end

%% Generate neighbor list

neighbors = zeros(NN,NN);

for j = 1:NN
    Nv(j) = 0;
    for i = 1:NN
         if(i ~= j)
             dist = (pos(j,1)-pos(i,1))^2 + (pos(j,2)-pos(i,2))^2 + (pos(j,3)-pos(i,3))^2;
             if(dist < tdist2) 
                Nv(j) = Nv(j) + 1; % number of neighbors
                neighbors(j,Nv(j)) = i ; % i is neighbor of j
             end
         end
    end
end

Ns = 0;
Sneighbors = zeros(NN); 
for i = 1:NN
     if(pos(i,3) < tdist) 
                Ns = Ns + 1; % number of neighbors
                Sneighbors(Ns) = i ; % i is neighbor of the surface
     end
end



ox = zeros(NN,1); % ox has the oxidated molecules
%% percolation loop

if(Ns == 0) 
    fporo = 0
    return
end

% surface
for i = 1:Ns
ox(Sneighbors(i)) = 1;
end

% particle-particle
change = 1; 
while(change == 1)
    change = 0;
    for j = 1:NN % loop over all particles
        if(ox(j) == 1) % particle j is oxidized
           for i = 1:Nv(j) % loop over neighbors of particle j
               if(ox(neighbors(j,i)) == 0) % molecule is a neighbor but is not oxidized
                   ox(neighbors(j,i)) = 1; 
                   change = 1;
               end
           end
        end
    end
end

counter = 0;
for j = 1:NN
    if(ox(j) == 1)
      hold all;
      scatter3(pos(j,1), pos(j,2), pos(j,3), 'r*') 
      counter = counter + 1;
    end
    if(ox(j) == 0)
      hold all;
      scatter3(pos(j,1), pos(j,2), pos(j,3), 'bo') 
    end
end
    
fporo = counter/NN
return
end  
