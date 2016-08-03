function TRI_2D_X(ibeg,iend,jbeg,jend,aa,bb,cc,dd)
for j=jbeg:jend
for i=ibeg:iend
if ibeg==1
cc(i,j)=cc(i,j)/bb(i,j);
dd(i,j)=dd(i,j)/bb(i,j);
else
cc(i,j)=cc(i,j)/(bb(i,j)-(aa(i,j)*cc(i-1,j)))
dd(i,j)=(dd(i,j)-aa(i,j)*dd(i-1,j))/(bb(i,j)-aa(i,j)*cc(i-1,j));
end
end
for i=iend-1:-1:ibeg
dd(i,j)=(dd(i,j)-cc(i,j)*dd(i+1,j));
end
end
return
end
function TRI_2D_Y(ibeg,iend,jbeg,jend,aa,bb,cc,dd)
for i=ibeg:iend
for j=jbeg:jend
if jbeg==1
cc(i,j)=cc(i,j)/bb(i,j);
dd(i,j)=dd(i,j)/bb(i,j);
else
cc(i,j)=cc(i,j)/(bb(i,j)-(aa(i,j)*cc(i,j-1)));
dd(i,j)=(dd(i,j)-aa(i,j)*dd(i,j-1))/(bb(i,j)-aa(i,j)*cc(i,j-1));
end
end
for j=jend-1:-1:jbeg
dd(i,j)=(dd(i,j)-cc(i,j)*dd(i,j+1));
end
end
return
