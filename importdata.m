FID = fopen('data.txt');
textdata = textscan(FID,'%f %f %f %f %s', 200, 'Delimiter',',');
data = cell2mat(textdata(:,1:4));
target=textdata{1,5};
[m,n] = size(target);
tr=[];
for k= 1:m
    a=target(k);
    if strcmp(a,'Iris-setosa')==1
        l=-1;
    elseif strcmp(a,'Iris-versicolor')==1
        l=0;
    else
        l=1;
    end
    tr=[tr;l];
end
clear a;

data=[data tr];
datatr=[data(1:1,:);data(51:51,:);data(101:101,:)];
datatst=[data(41:50,:);data(91:100,:);data(141:150,:)];


if max(abs(datatr(:,1:end-1)))> 1 
 
%Need to normalize
 
Norm_Input = datatr(:,1:end-1) / max(abs(datatr(:)));
 
else
 
Norm_Input = datatr;
 
end
bias=0;
wieghts=[0.1;0.1;0.1;0.1];
flag=0;
[sl,sw]=size(datatr);
while flag~=3
    wieghtsn=wieghts;
    for a=1:sl
        dtl=bias+wieghts(1)*datatr(a,1)+wieghts(2)*datatr(a,2)+wieghts(3)*datatr(a,3)+wieghts(4)*datatr(a,4);
        fn=sigmf(dtl,[1 0]);
        if fn~=datatr(a,5)
            bias=bias+0.6*(datatr(a,5)-fn);
            wieghts(1)=wieghts(1)+0.6*(datatr(a,5)-fn)*datatr(a,1);
            wieghts(2)=wieghts(2)+0.6*(datatr(a,5)-fn)*datatr(a,2);
            wieghts(3)=wieghts(3)+0.6*(datatr(a,5)-fn)*datatr(a,3);
            wieghts(4)=wieghts(4)+0.6*(datatr(a,5)-fn)*datatr(a,4);
            disp(wieghts)
        end
    end
    if wieghtsn==wieghts
        flag=flag+1;
    end
end