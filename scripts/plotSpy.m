clear all, close all
x = load("density.oct");
for i=1:size(x,1)
  for j=1:size(x,2)
    if x(i,j)<1.1
      x(i,j)=0;
    endif
  endfor
endfor
spy(x)

