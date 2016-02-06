function [mrf,map,iman] = denoise_MRF(a)

figure(1);
nb_levels = 3;

imagesc(a);
figure(2);
mean = 0;
sigma = 0.5;  
b = a + mean + (sigma*(randn(size(a,1),size(a,2))));
imagesc(b);
var = sigma^2;
c = fix(b+0.5);
iman = b;
%c = fix(nb_levels * rand(size(a,1),size(a,2)));
c = ( c < nb_levels ).*c + ( c > nb_levels) * (nb_levels -1);
figure(1);
imagesc(c);
newenerg_noise = 0;
pause;
energ = -(b-c-mean).^2/(2*var);
T = 1;
i = 0;
change = 1;
while (i < 50 & (change ~= 0 | i == 1))
   change = 0;
   for j = 2 : size(a,1)-1
      for k = 2 : size(a,2)-1
         newenerg = 0;
         newenerg_noise = 0;
         min_energ = 10000000000;
         best = 0;
         new = fix(b(j,k)+0.5);
         while (new == fix(b(j,k)+0.5))
            new = fix(nb_levels*rand(1,1)-0.5);
         end
         for l = -1 : 1
            for m = -1 : 1
               if ( l == 0 & m == 0)
                  newenerg_noise = (b(j+l,k+m)-new-mean)^2/(2*var);
               else
                  if (new == fix(c(j+l,k+m)))
                     newenerg = newenerg - 1/3;
                  else	
                     newenerg = newenerg + 1/3;
                  end 
               end
            end	
         end
         %newenerg = newenerg + newenerg_noise;    
         if (newenerg_noise + newenerg/T < min_energ)
            min_energ = newenerg_noise+newenerg/T;
            best = new;
         end
         if (min_energ <= energ(j,k))
            c(j,k) = best;
            energ(j,k) = min_energ;
            change = change + 1;
         end
      end	     
    end
    T = T / 1.5
   %T = T / (5*log((100+i))*log(100)); 
   figure(2);
   if (i == 0)
     map = c;
   end
   imagesc(c(2:size(a,1)-1,2:size(a,2)-1));
   drawnow;
   i = i+1
   change
 end	
 mrf = c;
