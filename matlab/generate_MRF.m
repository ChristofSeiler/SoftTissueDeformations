% This functions generates Markov Randon Fields
% Input:
%		image size as an array (ex [32 32]);
%		nb_levels: Number of levels of the random field
% Examples:
%		generate_MRF([32,32],2);
%		generate_MRF([32,32],5);

function a = generate_MRF(size,nb_levels)

% Start with a random image
a = fix(nb_levels*rand(size(1),size(2)));

disp('This starts with a random image');
imagesc(a); axis image;
drawnow;
% Start with an initial uniform energy field
energ = 100*ones(size(1),size(2));

disp('enter a key to start');
pause;

% Create neigboorhood weights (influence of the nieghbors on the value of the central pixel
for i = 1 : 3
    for j = 1 : 3
        mask(i,j) = 1/3;
    end
end
mask(1,1) = 0;
%mask = [
%  -1/3 1/3 -1/3;
%  1/3 0 1/3;
%  -1/3 1/3 -1/3];

% Initialise temperature to 1
T = 1;
cnt = 0;

% Checks every pixel in the image in turn and calculated its
% associated energy for a new grey level value based on its
% neighboors. If the new energy is smaller accept the change. If it
% is higher, accept the change with a probability dependent on the
% Temparature (annealing optimisation).

for  i = 1 : 100
    for j = 2 : size(1)-1
        for k = 2 : size(2)-1
            % Propose a new grey level value
            new = fix(nb_levels*rand(1,1));

            % Calculates the energy associated to the new configuration as:
            % E=
            % SUM(mask(i,j)*delta(pixel-neighbour))-SUM(mask(i,j)*(1-delta(pixel-neighbour)))
            % where delta is the delta function.
            diff = ((a(j-1:j+1,k-1:k+1)-new)~=0).*mask;
            eq = ((a(j-1:j+1,k-1:k+1)-new)==0).*mask;
            %newenerg = sum(sum(diff-eq));
            newenerg = sum(sum(eq-diff));
            % calulates the ratio between the old and the new local
            % energy for this pixel
            r = exp(-newenerg)/exp(-energ(j,k));
            % Accept the change if the new enegy is lower
            if  (r > 1)
                a(j,k) = new;
                energ(j,k) = newenerg;
                % Accept the chnage if r <1 drawing a random number
            elseif (rand(1,1) > 1-T)
                a(j,k) = new;
                energ(j,k) = newenerg;
            end
            % else reject the change
        end
    end
    disp('iteration='); disp(i);
    % Reduce the temperature for convergence
    T = T /log((100+i))*log(100)
    % Show new image every 10 iterations.
    if (mod(i,10) == 0)
        imagesc(a(2:size(1)-1,2:size(2)-1)); axis image;
        drawnow;
        disp('Press a new key to start');
        pause;
    end
end
