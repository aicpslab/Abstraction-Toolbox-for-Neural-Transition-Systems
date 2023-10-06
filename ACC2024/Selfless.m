function flag=Selfless(demos,interval,Num)
flag = 0;
counter=zeros(size(demos,2),1);
for i = 1:size(demos,2)
    counter(i) = 0;
    for j = 1: size(demos{1,i}.pos,2)
        if (interval(1,1)<demos{1,i}.pos(1,j))&&(interval(1,2)>demos{1,i}.pos(1,j))&&(interval(2,1)<demos{1,i}.pos(2,j))&&(interval(2,2)>demos{1,i}.pos(2,j))
            counter(i)=counter(i)+1;
        end
    end
end


for i = 1:size(demos,2)
   if counter(i)>Num
       flag=1;
   end
end

end