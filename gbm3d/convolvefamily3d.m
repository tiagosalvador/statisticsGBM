function cval = convolvefamily3d(ind,val,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cval = convolvefamily3d(ind,val,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convolves the characteristic function of union of grains whose
% level set representation is stored in "ind" and "val" input vars
% with the kernel "KERNELalpha" and "KERNELbeta" (global variable,
% assumed present).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global KERNELalpha KERNELbeta Z; % These global variables *must* have been already
                 % allocated before calling this function! Both are
                 % of dimension "dims".

%% Convert level set data to volume of fluid representation:
[x,y,z] = ind2sub(dims,ind);
vf = ls2vf3d(int32(x),int32(y),int32(z),val,Z,dims(1),dims(2),dims(3));

%% Carry out the convolution:
Kualpha = zeros(dims);
Kualpha(ind) = vf; % Characteristic function of the union of grains.
Kualpha = real(ifftn(fftn(Kualpha).* KERNELalpha));
Kubeta = zeros(dims);
Kubeta(ind) = vf; % Characteristic function of the union of grains.
Kubeta = real(ifftn(fftn(Kubeta).* KERNELbeta));
cval{1} = Kualpha;
cval{2} = Kubeta;

