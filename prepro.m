load("FTISxprt-20200310_flight2.mat")

fields = fieldnames(flightdata);
directory = "C:\Users\vladg\OneDrive\Documents\GitHub\Flight_Dyanmics_SVV\Data\";

for i = 1:numel(fields)
  current = flightdata.(fields{i});
  
  tab = current.data;
  units = current.units;
  desc = current.description;
  name = fields{i};
  %disp(current);
  
  filename = string(directory + string(fields{i})+".txt");
  
  fid = fopen(filename,"w");
 
  %fprintf(fid,string(description+"\n"));
  formatSpec = string(desc+"\n"+string(units)+"\n");
  fprintf(fid,formatSpec,"\n");
  
  for j=1:numel(tab)
      if isnan(tab(j))
          fprintf(fid,"0\n");
      else
          fprintf(fid,string(tab(j)+"\n"));
      end
  end
  fclose(fid);
  disp(j)
end

   
    
