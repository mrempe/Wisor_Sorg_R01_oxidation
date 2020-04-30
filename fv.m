function fvout = fv(v,Cm,k,vr,vt,u,Isyn,Iapp,Igap) 

fvout = (1./Cm).*(k.*(v-vr).*(v-vt)-u+Iapp-Isyn+Igap);
