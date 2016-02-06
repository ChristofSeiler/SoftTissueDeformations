% This functions generates Markov Randon Fields without simulated annealing or Gibbs sampling
% This is a cruder and faster version of the other prgram but wit no garantee of convergence.
% Input:
%		image size as an array (ex [32 32]);
%		nb_levels: Number of levels of the random field
% Examples:
%		generate_MRF([32,32],2);
%		generate_MRF([32,32],5);

function a = generate_MRF_quick(size,nb_levels)

% Start with a random image

a = fix(nb_levels*rand(size(1),size(2)));

disp('This starts with a random image');
imagesc(a); axis image;
drawnow;

% Start with an initial uniform energy field
energ = 100*ones(size(1),size(2));
disp('enter a key to start');
pause;

% This mask specifies the type of neighboorhood the Markov Random
% Field will use (influence of the nieghbors on the value of the central pixel)
% Here it is an anisotropic neighboorhood

mask = [
    1/3 1/3 1/3;
    1/3 0 1/3;
    1/3 1/3 1/3];

cnt = 0;

% Checks every pixel in the image in turn and calculated its
% associated energy for a new grey level value based on its
% neighboors. If the new energy is smaller accept the change.Else
% reject it.
for  i = 1 : 100
    for j = 2 : size(1)-1
        for k = 2 : size(2)-1
            newenerg = 0;
            % Generate a new random pixel value
            new = fix(nb_levels*rand(1,1));
            % Calculates the energy associated to the new configuration as:
            % E=
            % SUM(mask(i,j)*delta(pixel-neighbour))-SUM(mask(i,j)*(1-delta(pixel-neighbour)))
            % where delta is the delta function.
            diff = ((a(j-1:j+1,k-1:k+1)-new)~=0).*mask;
            eq = ((a(j-1:j+1,k-1:k+1)-new)==0).*mask;
            newenerg = sum(sum(diff-eq));
            % Directly compares energy and keep the best local energy.
            % This will converge to a local minimum but it is very fast!
            if (newenerg < energ(j,k))
                a(j,k) = new;
                energ(j,k) = newenerg;
            end
        end
    end
    disp('iteration='); disp(i);
    % Show evolution for a few iteration
    if (mod(i,10) == 0)
        imagesc(a(2:size(1)-1,2:size(2)-1)); axis image;
        drawnow;
        disp('Press a new key to start');
        pause;
    end
end
