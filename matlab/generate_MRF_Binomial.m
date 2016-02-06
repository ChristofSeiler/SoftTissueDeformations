% This functions generates a Binomial Markov Randon Fields
% Input:
%		image size as an array (ex [32 32]);
%		Betas: The 3 parameters of the binomial field
% Examples:
%		generate_MRF_Binomial([32,32],[-7.8 3.9 3.9])
%
% This implementation can only generate binary images
% Try implementing a grey level version as an exercise

function a = generate_MRF_Binomial(size, Betas)

% Start with a random image. The initialisation should not change
% the result
disp('This starts with a random image');

a = fix(rand(size(1),size(2))+0.5);
imagesc(a); axis image;
drawnow;
disp('enter a key to start');
pause;

% Enter the Beta parameters of the model:
B0 = Betas(1)
B1 = Betas(2)
B2 = Betas(3)

cnt = 0;

% Checks every pixel in the image in turn and calculated its
% associated energy Nt for a new grey level value based on its
% neighboors and the associated probability. Accept the change based on this 
% probability drawaing a random number

for  i = 1 : 100
    for j = 2 : size(1)-1
        for k = 2 : size(2)-1
            % Calculates the energy associated to the current configuration as:
            % Nt=
            % B0+B1*delta(pixel-neighbour on x)+B2*delta(pixel-neighbour on y)
            % where delta is the delta function.

            Nt = B0+B1*(a(j-1,k)+a(j+1,k))+B2*(a(j,k+1)+a(j,k-1));

            PNt = exp(Nt)/(1+exp(Nt));
            
            % This is the version implemented in the notes. Please
            % uncomment to try it out.
            
            %if (PNt > 0.5)
            %    a(j,k) = 1;
            %elseif (PNt < 0.5)
            %   a(j,k) = 0;
            %else
            %    a(j,k) = fix(rand(1,1)+0.5);
            %end
            
            % This is a more advanced and random version with a lower
            % convergence rate but better statistical properties
            % Draw a random number
            r = rand(1,1);
            % Make pixel 1 if random number < PNt , 0 otherwise
            if  (r  < PNt)
                a(j,k) = 1;
            else
                a(j,k) = 0;
            end
        end
    end
    disp('iteration='); disp(i);
    % Show new image every 10 iterations.
    if (mod(i,10) == 0)
        imagesc(a(2:size(1)-1,2:size(2)-1)); axis image;
        drawnow;
        disp('Press a new key to start');
        pause;
    end
end
