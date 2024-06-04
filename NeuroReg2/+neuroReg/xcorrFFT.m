function Im_out = xcorrFFT(Im1,Im2)
% xcorrFFT provides a fast cross correlation calculation.
% Im_out = xcorrFFT(Im1,Im2)
% INPUT:
% Im1 should be Lx1-by-Lx2 2D matrix
% Im2 should be either a Lx2-by-Ly2 2D matrix or a [Lx2,Ly2,Lz2] 3D matrix.
% In the 2D case, it calculates the XCC between the two images (the same output as filter2(Im1,Im2)).
% In the 3D case, it calculates the XCC between Im1 and each layer of Im2.
% OUTPUT:
% Im_out: Size: [Lx1,Ly1,Lz2]
% [Lx1,Ly1]=size(Im1)
% [Lx2,Ly2,Lz2] = size(Im2)


[Lx1,Ly1]=size(Im1);
[Lx2,Ly2,Lz2] = size(Im2);
NewLx1 =Lx1; NewLy1 =Ly1; 
PadX = 1:Lx1;
PadY= 1:Ly1;
if Lx2>Lx1 | Ly2>Ly1  
    NewIm1 = zeros(max(Lx2,Lx1),max(Ly2,Ly1));
   [NewLx1,NewLy1]=size(NewIm1);
    dx = Lx1 - Lx2;
    dy = Ly1 - Ly2;
    
    if dy<0
        shiftY = -round(dy/2);
        PadY = [1+ shiftY :NewLy1-(shiftY-mod(dy,2))];
    end
    if dx<0
        shiftX = -round(dx/2);
         PadX = [1+shiftX:NewLx1-(shiftX-mod(dx,2))];
    end
    NewIm1(PadX,PadY) = Im1;
    warning('slice size exceeds stack size for this rotation')
    Im1 = NewIm1;
end


% middle points. for padding
Mx1 = round(NewLx1/2);
My1 = round(NewLy1/2);
Mx2 = round(Lx2/2);
My2 = round(Ly2/2);
x_start = Mx1-Mx2+1;
x_end = Mx1-Mx2+Lx2;
y_start = My1-My2+1;
y_end = My1-My2+Ly2;
f1 = fft2(Im1);
Im_out = zeros(Lx1,Ly1,Lz2);
for i = 1:Lz2
    Im2_pad = zeros(NewLx1,NewLy1);
    Im2_temp = Im2(:,:,i);
    Im2_pad(x_start:x_end,y_start:y_end)=Im2_temp;
    f2 = fft2(Im2_pad);
    value_temp_slice = fftshift(ifft2(f1 .* conj(f2)));
    value_temp_slice = value_temp_slice(PadX,PadY);
    Im_out(:,:,i) = value_temp_slice;
end



end

