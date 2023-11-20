rm(list=ls())


args=commandArgs(trailingOnly = TRUE)

path_script_files<-args[1]
path_simulation_files<-args[2]
path_matrix_files<-args[3]
path_lineage_files<-args[4]
#simulation_type<-args[5]


#print (path_simulation_files)

scripts<-list.files(path = path_script_files,pattern = "RF")
inp<-list.files(path=path_simulation_files,pattern ="\\.phen$")

for (j in scripts){
  
   for (i in inp){
	 simulation_type<-(strsplit(i, "[_]")[[1]][2])
  	
      #system(sprintf("Rscript %s %s %s %s",j,i,path_matrix_files,path_lineage_files))
	system(command=sprintf("Rscript %s %s%s %s %s %s",j,path_simulation_files,i,path_matrix_files,path_lineage_files,simulation_type))
      	#print(sprintf("Rscript %s %s%s %s %s %s",j,path_simulation_files,i,path_matrix_files,path_lineage_files,simulation_type))
	}
  }




