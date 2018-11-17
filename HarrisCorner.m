function [I_out] = HarrisCorner(I_in,G_sgm,Gs_sgm,k_value,threshold)

%Create a guassian mask and divide it to 2 directions.
Gua = fspecial('gaussian',6 * G_sgm, G_sgm);
Gua_x = conv2( Gua, [-1 1], 'same');
Gua_y = conv2( Gua, [-1;1], 'same');

%Compute Ix,Iy.
I_in_x = conv2(I_in, Gua_x, 'same');
I_in_y = conv2(I_in, Gua_y, 'same');
Ix_2 = I_in_x.^2;
Iy_2 = I_in_y.^2;
IxIy = I_in_x .* I_in_y;

%Compute Sx,Sy.
Gus = fspecial('gaussian',6 * Gs_sgm, Gs_sgm);
Sx = conv2(Ix_2, Gus, 'same');
Sy = conv2(Iy_2, Gus, 'same');
Sxy = conv2(IxIy, Gus, 'same');

[row,col] = size(I_in);

%Compute R for every pixel;
R_matrix = zeros(row,col);

for i = 1:1:row
    for j = 1:1:col
        
        %Create the H matrix;
        H = [Sx(i,j) Sxy(i,j) ; Sxy(i,j) Sy(i,j)];
        
        %Find the R value;
        R_ini = det(H) - k_value * (trace(H) ^ 2);
        
        %Check R to Threshold;
        if R_ini >= threshold
            R_fin  = R_ini;
        else
            R_fin = 0;
        end
        
        %Input to R_matrix;
        R_matrix(i,j) = R_fin;
    end
end

%Compute nonmax suppression.
I_out = NonmaxSuppression(R_matrix);

end

function [ma_new] = NonmaxSuppression(ma_ori)

%Zero padding the origin matrix;
[row,col] = size(ma_ori);
ma_0pad = [0 zeros(1,col) 0; zeros(row,1) ma_ori zeros(row,1);0 zeros(1,col) 0];

for i=2:1:row+1
    for j = 2:1:col+1
      neigh_max = max([ma_0pad(i-1,j-1) ma_0pad(i,j-1) ma_0pad(i+1,j-1) ma_0pad(i-1,j) ma_0pad(i+1,j) ma_0pad(i-1,j+1) ma_0pad(i,j+1) ma_0pad(i+1,j+1)]); 
      
      %If has a larger neighbor, set to 0.
      if ma_0pad(i,j) < neigh_max
          ma_0pad(i,j) = 0;
      end
      
      ma_new = ma_0pad(2:row+1, 2:col+1);
    end
end

end

