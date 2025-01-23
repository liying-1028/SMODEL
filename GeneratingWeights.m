function [w]=GeneratingWeights(poolsget,a,tau,highorder,lambda)
[n1,n2]=size(poolsget);
b=a;
if highorder == 1
    %Find the second-order similarity matrix (Gaussian kernel)
    [m1,~]=size(b);
    A=0;
    for i=1:m1
        A=A+norm(b(i,:)-mean(b,1),2)^2;
    end
    sigma2=A/(m1-1);
    d=zeros(n1,n1);
    for i=1:n1
        for j=1:n1
            d(i,j)=exp(-norm(b(i,:)-b(j,:),2)^2/(2*sigma2));
        end
    end
    %Weighted by a first-order similarity matrix
    [~,a2]=sort(a,2);
    [~,d2]=sort(d,2);
    c=a2(:,(n1-tau+1):n1);
    d=d2(:,(n1-tau+1):n1);
    s=zeros(n1,n2);
    t=zeros(n1,n2);
    for i=1:n2
        for j=1:n1
            for k=1:tau
                if poolsget(j,i)==poolsget(c(j,k),i)
                    s(j,i)=s(j,i)+1;
                else s(j,i)=s(j,i);
                end
                if poolsget(j,i)==poolsget(d(j,k),i)
                    t(j,i)=t(j,i)+1;
                else t(j,i)=t(j,i);
                end
            end
            s(j,i)=s(j,i)/(length(find(poolsget(:,i)==poolsget(j,i)))+tau-s(j,i));
            t(j,i)=t(j,i)/(length(find(poolsget(:,i)==poolsget(j,i)))+tau-t(j,i));
        end
    end
    %Normalization
    W1=s./max(s,[],2);
    W2=t./max(t,[],2);
    %Weighted integration
    w=lambda*W1+(1-lambda)*W2;
else
    %Weighted by a first-order similarity matrix
    [~,a2]=sort(a,2);
    c=a2(:,(n1-tau+1):n1);
    s=zeros(n1,n2);
    for i=1:n2
        for j=1:n1
            for k=1:tau
                if poolsget(j,i)==poolsget(c(j,k),i)
                    s(j,i)=s(j,i)+1;
                else s(j,i)=s(j,i);
                end
            end
            s(j,i)=s(j,i)/(length(find(poolsget(:,i)==poolsget(j,i)))+tau-s(j,i));
        end
    end
    %归一化
    W1=s./max(s,[],2);
    %W1=s;
    w=W1;
end

