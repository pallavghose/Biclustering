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
alp=0.7;
%del ismaximum accepted mean squared residue should be > or = to 0
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
    while flag~=2
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
       l2=1;
       resd=[];
       
       for i=bi_row
           resdr=[];
           l2=0;
           for j=bi_col
               temp=(bicl(i,j)-(row_maj(1,i)+col_maj(1,j)))+tot_maj;
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
               [I,J]=find(bi_row==i);
               bi_row(:,J)=[];
               row_resd(i,:)=[];
               m=m-1;
               i=i-1;
               brb=0;
           end
           
           i=i+1;
           if i==m
               fl=fl+1;
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
            row_maj=[row_maj [i;mn]];
        end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj [i;mn]];
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
           x= find(row_maj(1,:)==i);
           for j=bi_col
               y= find(col_maj(1,:)==j);
               
               temp=(bicl(i,j)-(row_maj(2,x)+col_maj(2,y)))+tot_maj;
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
       [m,~]=size(col_resd);
       while fl~=1
           
           if col_resd(i,2)>(alp*MSR)
               [I,J]=find(bi_col==i);
               bi_row(:,J)=[];
               col_resd(i,:)=[];
               m=m-1;
               i=i-1;
               brb=0;
           end
           i=i+1;
           if i==m
               fl=fl+1;
           end
           
       end
        bi_col(12:13)=[];
       if brb==1
         flag=1;
       end
       break;
       
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
            row_maj=[row_maj [i;mn]];
        end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj [i;mn]];
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
           x= find(row_maj(1,:)==i);
           for j=bi_col
               y= find(col_maj(1,:)==j);
               
               temp=(bicl(i,j)-(row_maj(2,x)+col_maj(2,y)))+tot_maj;
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
       l=0;
       while l~=1600
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
            row_maj=[row_maj [i;mn]];
        end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj [i;mn]];
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
           x= find(row_maj(1,:)==i);
           for j=bi_col
               y= find(col_maj(1,:)==j);
               
               temp=(bicl(i,j)-(row_maj(2,x)+col_maj(2,y)))+tot_maj;
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
       for i=1:l1
           for j=1:l2
             resd(i,j)=resd(i,j)*resd(i,j);
           end
       end
       [~,m]=size(bi_row);
        [~,n]=size(bi_col);
       l1=1;
       row_res=[];
       row_resd=[];
       for i=bi_row
           row_res=sum(resd(l1,:))/n;
           row_resd=[row_resd;[i row_res]];
           l1=l1+1;
       end
       l1=l1-1;

       l2=1;
       col_res=[];
       col_resd=[];
       for i=bi_col
           col_res=sum(resd(:,l2))/n;
           col_resd=[col_resd;[i col_res]];
           l2=l2+1;
       end
       l2=l2-1;
       
       %...................row resd and col resd found.........
       %we will find which is the row or col having the biggest residue value
       %we will then delete it from the bicluster
       row_max=max(row_resd(:,2));
       col_max=max(col_resd(:,2));
       row_maxp=find(row_resd(:,2)==row_max);
       col_maxp=find(col_resd(:,2)==col_max);
       
      if row_max>col_max
           bi_row(row_maxp)=[];
      end
      if l2>7
          if row_max>col_max
              bi_row(row_maxp)=[];
          else
              bi_col(col_maxp)=[];
          end
      else
           bi_row(row_maxp)=[];
      end
       l=l+1;
       
       
    end
    %2nd algorithm ends
    
    %3rd algorithm starts..............................................
    
   
    row_rem=setdiff([1:2884],bi_row);
    col_rem=setdiff([1:17],bi_col);
    %add rows......................................................
    for i=row_rem
       [~,m]=size(bi_row);
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
            row_maj=[row_maj [i;mn]];
         end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj [i;mn]];
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
           x= find(row_maj(1,:)==i);
           for j=bi_col
               y= find(col_maj(1,:)==j);
               temp=0;
               temp=bicl(i,j)-(row_maj(2,x)-col_maj(2,y))+tot_maj;
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
       bi_row=[bi_row i];
       %___________________________________________________________________
        
        
       
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
            row_maj=[row_maj [i;mn]];
        end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj [i;mn]];
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
           x= find(row_maj(1,:)==i);
           for j=bi_col
               y= find(col_maj(1,:)==j);
               
               temp=(bicl(i,j)-(row_maj(2,x)+col_maj(2,y)))+tot_maj;
               resdr=[resdr temp];
              l2=l2+1;
           end
           resd=[resd;resdr];
           l1=l1+1;
       end
       
       %Find H(I,J) mean squared residue(MSR) of the bicluster is..........
       MSR1=0;
       for i=1:l1
           for j=1:l2
             MSR1=MSR1+(resd(i,j)*resd(i,j));
           end
       end
       MSR1=MSR1/(m*n);
       if(MSR1>MSR)
          bi_row(m)=[];
       end
    end
    
    %___________________________________________________________
    %add Columns
    
    
     for i=col_rem
        bir_n=bi_row;
        bic_n=bi_col;
        row_maj=[];
        col_maj=[];
        tot_maj=[];
        [~,m]=size(bi_row);
        [~,n]=size(bi_col);
        for i=bi_row
            mn=sum(bicl(i,:))/n;
            row_maj=[row_maj [i;mn]];
        end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj [i;mn]];
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
           x= find(row_maj(1,:)==i);
           for j=bi_col
               y= find(col_maj(1,:)==j);
               
               temp=(bicl(i,j)-(row_maj(2,x)+col_maj(2,y)))+tot_maj;
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
       
      
       bi_col=[bi_col i];
      
        %___________________________________________________________________
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
            row_maj=[row_maj [i;mn]];
        end
        %find column variance aIj
        for i=bi_col
            mn=sum(bicl(:,i))/m;
            col_maj=[col_maj [i;mn]];
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
           x= find(row_maj(1,:)==i);
           for j=bi_col
               y= find(col_maj(1,:)==j);
               
               temp=(bicl(i,j)-(row_maj(2,x)+col_maj(2,y)))+tot_maj;
               resdr=[resdr temp];
              l2=l2+1;
           end
           resd=[resd;resdr];
           l1=l1+1;
       end
       
       MSR1=MSR1/(m*n);
       if(MSR1>MSR)
          bi_col(n)=[];
       end
    end
    
    
    %masking the bicluster with random generated values
end

