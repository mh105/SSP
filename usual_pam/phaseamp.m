function pa=phaseamp(a,p,pbins)
pa=zeros(length(pbins)-1,1);
%For each of the phase bins
for i=1:length(pbins)-1
    %Get the mean of the amplitude in that phase bin
    if ~isempty(a(p>pbins(i) & p<= pbins(i+1)))
    pa(i)=nanmean(a(p>pbins(i) & p<= pbins(i+1)));
    end
    
end

%Normalize by the amplitude
%pa=pa./nanmean(a);
