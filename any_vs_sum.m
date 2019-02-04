% test speed of any.m vs sum.m

s=zeros(1,1e6);

I = rand(1e6,500);
I(find(I>0.5)) = 1;
I(find(I<=0.5)) = 0;




s=zeros(1,1e6);
tic
for i=1:1e6
	if any(I(i,:))
	s(i) = 1;
	end
end
T=toc;
disp(['time to run using if any: ', num2str(T)])



tic
for i=1:1e6
	if sum(I(i,:))>0
	s(i) = 1;
	end
end
T1=toc;
disp(['time to run using sum: ', num2str(T1)])
