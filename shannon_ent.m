function hx = shannon_ent(x,r)
%hx = shannon_ent(x,r)
%
%This function will calculate Shannon's entropy as the sum of p(xi)*logb
%INPUTS: signals X and tolerance r [default is 0.001].
%      
%OUTPUT: hx entropy
%
%
%Define function name
% FuncName = 'shannon_ent';
%
% Written by Dr. Samhita Rhodes
% Adapted by Natalie Tipton

%Check input parameters
if(nargin<2)
	fprintf(1,'Must give all input parameters\nType help Entropy\n\n');
	return
end
if(~exist('r') | isempty(r))    % default value of r = 0.001
	r=.001;
end;

%Obtaining Entropy
x=(x-mean(x))./std(x);  % normalize data 
N = length(x);          % find length of data
x1 = unique(x);         % returns x with no repetititive values
hx = 0;                 % initialize entropy variable
for i=1:length(x1)
    nx(i) = sum(abs(x1(i)-x)<=r);   % finds all values within r of original normalized data
end
px = nx./sum(nx);               % calculated pdf 
px1 = px(find(px));             % gets locations of each pdf value
hx = hx + sum(px1.*log2(px1));  % calculates entropy
hx = -hx;                       % makes entropy positive