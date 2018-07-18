% script to test speed of randn vs normrnd

N=5000

disp('randn: ')
tic
for i=1:N
	r = randn(2,1);
end
toc

disp('normrnd')
tic
for i=1:N
	s(1) = normrnd(0,1);
	s(2) = normrnd(0,1);
end
toc



