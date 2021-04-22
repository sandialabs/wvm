function VR=pvl_WVM_compute_VR(dist,tmscales,CS)

%computes variablity reduction (VR), which is a funciton of the distances
%between sites (dist), the timescales, and the clouds speed (CS)
 fn=@(x)abs((x.^2-x)./2-length(dist));
[len_dist]=round(fminsearch(fn,sqrt(length(dist))));

for i=1:length(tmscales)
    r1=exp(-dist./(CS/2*tmscales(i)));
    rhosum=sum(r1);
   
    VR(i)=1./((2*rhosum+len_dist)./(len_dist.^2));
end
