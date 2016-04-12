data = dlmread('yeast.matrix.txt');
bi=[];
%r is the value of the missing data in the dataset it is -1 for yeast data
r=-1;
%replacing missing data with random nos ranging from 0 to 800
for i=1:2884
    for j=1:17
        if data(i,j)==r
            data(i,j)=round(rand(1)*800);
        end
    end
end
%changing the value of n will give you the no of biclusters
n=1;
%alp is the threshold for multiple node deletion and is greater than 1
alp=1.1;
%del is the maximum accepted mean squared residue should be > or = to 0
del=300;
bicl=data;
bi_row=[];
bi_col=[];
bi_res=[];
for itr=1:n
%beginning of all the three algorithms
    bi_row=1:2884;
    bi_col=1:17;
    %1st algorithm starts........................................
    flag=0;
    while flag~=1
        %calculate H(I,J),aIj,aiJ,aij..................................
        %find row variance aiJ
        bir_n=bi_row;
        bic_n=bi_col;
        row_maj=[];
        col_maj=[];
        tot_maj=[];
        [~,m]=size(bi_row);
        [~,n]=size(bi_col);
        for i=bi_row
            mn=sum(bicl(i,:))/n;
            row_maj=[row_maj mn];
        end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj mn];
        end
       %find table variance aIJ
       tot_maj=0;
       for i=bi_row
           for j=bi_col
                tot_maj=tot_maj+bicl(i,j);
           end
       end
       tot_maj=tot_maj/(m*n);
       %find the residue a_{ij} - a_{i J} - a_{I j} + a_{I J}.........
       l1=0;
       l2=0;
       resd=[];
       
       for i=bi_row
           resdr=[];
           l2=0;
           for j=bi_col
               temp=(bicl(i,j)-(row_maj(1,i)+col_maj(1,j))+tot_maj);
               resdr=[resdr temp];
              l2=l2+1;
           end
           resd=[resd;resdr];
           l1=l1+1;
       end
       %Find H(I,J) mean squared residue(MSR) of the bicluster is..........
       MSR=0;
       for i=1:l1
           for j=1:l2
             MSR=MSR+resd(i,j)*resd(i,j);
           end
       end
       MSR=MSR/(m*n);

       if MSR<=del
          break; 
       end
       %remove multiple rows 1/J row[sum(residue^2)]> alp*MSR
       row_resd=[];
       %Find the row resd of each row
       for i=1:m
           row_res=0;
           for j=1:n
               row_res=row_res + resd(i,j)*resd(i,j);
           end
           row_res=row_res/n;
           row_resd=[row_resd;i row_res];
       end
       %delete rows
       i=1;
       fl=0;
       [m,~]=size(row_resd);
       while fl~=1   
           if row_resd(i,2)>(alp*MSR)
               bi_row=bi_row(1,bi_row(:,2)==row_resd(i,1));  
               row_resd(i,:)=[];
               m=m-1;
               i=i-1;
               brb=0;
           end
           
           i=i+1;
           if i==m
               fl=1;
           end
       end
       
       %calculate H(I,J),aIj,aiJ,aij..................................
        %find row variance aiJ
        bir_n=bi_row;
        bic_n=bi_col;
        row_maj=[];
        col_maj=[];
        tot_maj=[];
        [~,m]=size(bi_row);
        [~,n]=size(bi_col);
        for i=bi_row
            mn=sum(bicl(i,:))/n;
            row_maj=[row_maj mn];
        end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj mn];
        end
       %find table variance aIJ
       tot_maj=0;
       for i=bi_row
           for j=bi_col
                tot_maj=tot_maj+bicl(i,j);
           end
       end
       tot_maj=tot_maj/(m*n);
       %find the residue a_{ij} - a_{i J} - a_{I j} + a_{I J}.........
       l1=0;
       l2=0;
       resd=[];
       
       for i=bi_row
           resdr=[];
           l2=0;
           for j=bi_col
               temp=(bicl(i,j)-(row_maj(1,i)+col_maj(1,j))+tot_maj);
               resdr=[resdr temp];
              l2=l2+1;
           end
           resd=[resd;resdr];
           l1=l1+1;
       end
       %Find H(I,J) mean squared residue(MSR) of the bicluster is..........
       MSR=0;
       for i=1:l1
           for j=1:l2
             MSR=MSR+resd(i,j)*resd(i,j);
           end
       end
       MSR=MSR/(m*n);

       if MSR<=del
          break; 
       end
       %remove multiple columns 1/J col[sum(residue^2)]> alp*MSR
       col_resd=[];
       for i=1:n
           col_res=0;
           for j=1:m
               col_res=col_res+ resd(j,i)*resd(j,i);
           end
           col_res=col_res/m;
           col_resd=[col_resd;i col_res];
       end
       
        %delete columns
       i=1;
       fl=0;
       [m,~]=size(row_resd);
       while fl~=1
           
           if col_resd(i,2)>(alp*MSR)
               bi_col=bi_col(1,bi_col(:,2)==col_resd(i,1));  
               col_resd(i,:)=[];
               m=m-1;
               i=i-1;
               brb=0;
           end
           i=i+1;
           if i==m
               fl=1;
           end
       end
       if brb==1
         flag=1;
       end
    end
    %1st algorithm ends
    
    %2nd algorithm starts
    %Algoritm 2 to be implemented delete single rows and columns to improve
   %the quality of the bicluster
    %recalculate H(I,J),aIj,aiJ,aij..................................
       
        %find row variance aiJ
        bir_n=bi_row;
        bic_n=bi_col;
        row_maj=[];
        col_maj=[];
        tot_maj=[];
        [~,m]=size(bi_row);
        [~,n]=size(bi_col);
        for i=bi_row
            mn=sum(bicl(i,:))/n;
            row_maj=[row_maj mn];
        end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj mn];
        end
       %find table variance aIJ
       tot_maj=0;
       for i=bi_row
           for j=bi_col
                tot_maj=tot_maj+bicl(i,j);
           end
       end
       tot_maj=tot_maj/(m*n);
       %find the residue a_{ij} - a_{i J} - a_{I j} + a_{I J}.........
       l1=0;
       l2=0;
       resd=[];
       
       for i=bi_row
           resdr=[];
           l2=0;
           for j=bi_col
               temp=(bicl(i,j)-(row_maj(1,i)+col_maj(1,j))+tot_maj);
               resdr=[resdr temp];
              l2=l2+1;
           end
           resd=[resd;resdr];
           l1=l1+1;
       end
       %Find H(I,J) mean squared residue(MSR) of the bicluster is..........
       MSR=0;
       for i=1:l1
           for j=1:l2
             MSR=MSR+(resd(i,j)*resd(i,j));
           end
       end
       MSR=MSR/(m*n);

       if MSR<=del
          break; 
       end
       %check residue value of the largest row and column and delete
       %accordingly....................
       flag=1;
       if flag==0;
       while MSR>del
           %find d(i) of all rows................................
           row_resd=[];
           for i=1:m
               row_res=0;
               for j=1:n
                   row_res=row_res + resd(i,j)*resd(i,j);
               end
               row_res=row_res/n;
               row_resd=[row_resd;i row_res];
           end
           %find d(j) of all columns..............................
           col_resd=[];
           for i=1:n
               col_res=0;
               for j=1:m
                   col_res=col_res+ resd(j,i)*resd(j,i);
               end
               col_res=col_res/m;
               col_resd=[col_resd;i col_res];
           end
           max_r=max(row_resd);
           max_c=max(col_resd);
           [~,posr]=max(row_resd);
           [~,posc]=max(col_resd);
           if max_r>max_c
              bicl(posr,:)=[]; 
           
           else
              bicl(:,posc)=[];    
           end
           [m,n]=size(bicl);
           %find row variance aiJ
            for i=1:m
                mn=sum(bicl(i,:))/n;
                row_maj=[row_maj mn];
            end
            %find column variance aIj
            for i=1:n
                mn=sum(bicl(i,:))/m;
                col_maj=[col_maj mn];
            end
           %find table vaeiance aIJ
           tot_maj=sum(sum(bicl));
           tot_maj=tot_maj/(m*n);
           %find the residue a_{ij} - a_{i J} - a_{I j} + a_{I J}.........
           for i=1:m
               for j=1:n
                 resd(i,j)=bicl(i,j)-(row_maj(i)+col_maj(j))+tot_maj;
               end
           end
           %Find H(I,J) mean squared residue(MSR) of the bicluster is......
           MSR=0;
           for i=1:m
               for j=1:n
                 MSR=MSR+resd(i,j)*resd(i,j);
               end
           end
           MSR=MSR/(m*n);
           if MSR<=del
              break; 
           end
       end
       end
    %2nd algorithm ends
    %3rd algorithm starts
    
%end of all the three algorithms
end

