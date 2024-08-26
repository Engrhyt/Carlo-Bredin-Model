format long
pkg load signal
pkg load optim
pkg load statistics
pkg load matgeom
pkg load control
directory=['/home/engrhyt/Desktop/carlofitfast/process2']
workset=[directory,'/work']

savefileplace=[directory,'/Bredin']

save_analyze_place=[directory,'/parameter_data']
save_analyze_file_place=[directory,'/parameter_data/parameter_table/Traw12/']
save_analyze_pic_place=[directory,'/parameter_data/parameter_pic']

for table_index=12:15
short_index_start=5;
parameter_analysis_filterCsfit
endfor

%for table_index=1%:4%:-1:1
%parameter_analysis
%endfor

%for gas_select=1
%density_analysis
%endfor

