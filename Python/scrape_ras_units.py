import fnmatch
import os
import pandas as pd
import os.path
from os import path
from re import search
from pathlib import Path
import pyproj

# TODO: not sure if this is still being used but let's add any new code to the existing workflow
# TODO: remove any hardcoded paths ('<PATH>\inland_routing\OWP_ras_share') to a function arg
matches = pd.DataFrame(columns = ['hdf_path', 'model_unit','full_path'])
for root, dirs, files in os.walk('<PATH>\inland_routing\OWP_ras_share'):
    for gfile in fnmatch.filter(files, '*.g[0-9][0-9].hdf'):
        for proot, pdirs, pfiles in os.walk(os.path.join(root)):
            for pfile in fnmatch.filter(pfiles, '*'+str(os.path.basename(gfile)[:-8])+'.prj'):
                with open(os.path.join(root,pfile)) as f:
                    first_file_line = f.read()
                
                # skip projection files
                if any(x in first_file_line for x in ['PROJCS','GEOGCS','DATUM','PROJECTION']):
                    continue
                
                if search("SI Units", first_file_line):
                    unit_scrape ='meter'
                    #matches = matches.append({'hdf_path':gfile,'model_unit':'meter'},ignore_index = True)
                    
                elif search("English Units", first_file_line):
                    unit_scrape ='foot'
                    #matches = matches.append({'hdf_path':gfile,'model_unit':'foot'},ignore_index = True)
                else:
                    unit_scrape ='manual'
                    matches = matches.append({'hdf_path':gfile,'model_unit':'manual'},ignore_index = True)
                    
                matches = matches.append({'hdf_path':gfile,'model_unit':unit_scrape,'full_path':os.path.join(root,gfile)},ignore_index = True)
                #matches = matches.append({'hdf_path':gfile,'proj':os.path.join(proot, pfile)},ignore_index = True)
                

manual_db = pd.read_csv ('<PATH>\inland_routing\data\survey_db.csv')
export = matches.merge(manual_db, on='hdf_path', how='left')

proj_paths = pd.read_csv('<PATH>\inland_routing\data\models_with_proj_data_edited.csv')
proj_paths = pd.read_excel('<PATH>\inland_routing\data\models_with_proj_data_edited.xlsx')
proj_paths['hdf_path'] = proj_paths[['model_name', 'model_extension']].apply(lambda x: '.'.join(x), axis=1)

export = export.merge(proj_paths, on='hdf_path', how='left')
#print(export)
#print('pre proj str:')
#print(export[export.proj.dropna()])

def return_proj_str(path_to_shp):
    try:
        shp = open(Path(path_to_shp),"r").readline()
    except:
        #print("could not open"+str(path))
        return False
        str_out = str(pyproj.CRS(shp))
    return str_out

#export['proj'] = export.apply(lambda row : return_proj_str(row['shapefile_path']), axis = 1)

#export.to_csv('<PATH>\inland_routing\data\survey_db_units.csv', index=False)  

survey_db = pd.read_excel('<PATH>\inland_routing\data\survey_db_units.xlsx')
out_dir = Path(r'<PATH>\inland_routing\data\20220309')

# apply to every valid row in DF

if not os.path.exists(os.path.join(out_dir,'tmp')):
    os.mkdir(os.path.join(out_dir,'tmp'))

all_files = os.listdir(os.path.join(out_dir,"tmp"))
csv_files = list(filter(lambda f: f.endswith('.csv'), all_files))
combined_csv = pd.concat([pd.read_csv(f) for f in csv_files])
combined_csv.to_csv(os.path.join(out_dir,"ras_xyz.csv"), index=False, encoding='utf-8-sig')

