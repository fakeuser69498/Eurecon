#################### ModelNet40-Eurecon Preprocessing and Augmentation Tool For Train/Test Dataset Splits ####################
# Python ENV: EL   

import os.path as path
import sys
import os
import glob
import shutil
#from termcolor import colored
from tqdm import tqdm
from time import sleep
#from torch_geometric.io import read_off, write_off
#from torch_geometric.data import Data
#from torch_geometric.transforms import NormalizeScale
import datetime

rmsd = 0.5
base_dir = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(base_dir)

now = datetime.datetime.now()
#file = open('/home/apm/Desktop/aug_20_rmsd_' + str(rmsd) + '.txt', 'w')
#file = open('/home/apm/Desktop/aug_rand_20_' + str(rmsd) + '_' + str(now) + '.txt', 'w')
file = open('/beegfs/gpfs0/a.morozov/EureconAugmentedDatasets/MN/Non-scaled/Eurepackage/aug_20_relmin302_recentered.txt', 'w')
#file = open('/home/apm/Desktop/EL/PointNet/nbs/aug_20_relative_min3_0.01_' + str(now) + '.txt', 'w')
#input_directory = '/home/apm/Desktop/EL/PointNet/nbs/ModelNet40A20RMSD'+ str(rmsd).replace(".", "")
input_directory = '/beegfs/gpfs0/a.morozov/EureconAugmentedDatasets/MN/Non-scaled/Eurepackage/2MN40RF'
path_list = os.listdir(input_directory)
result1 = [] 
result2 = []
counter1 = 0
counter2 = 0
for path in path_list:
    result1.append(input_directory + '/' + path + '/train')
    result2.append(input_directory + '/' + path + '/test')

counter_after1 = 0
for object_path1 in result1:
    files1 = [f for f in glob.glob(object_path1 + "/*.off", recursive=True)]
    counter1 += len(files1)

#################### File Fixing in Train Folders ####################

    for i in tqdm(files1, desc = 'Progress (Train Files)'):
        # print(colored("> Reading file: ", 'red', attrs = ['bold']) + i);
        # in_file = open(i, "r");
        # all_lines = in_file.readlines()
        # first_line = all_lines[0];
        # tokens = first_line.split();
        # # should have only one token: OFF
        # if(len(tokens) > 1):
        #     # Get the value with the 'OFF'
        #     tokens[0] = tokens[0].split("OFF")[1];
        #     edited_line = tokens[0] + " " + tokens[1] + " " + tokens[2] + "\n";
        #     all_lines[0] = edited_line;
        #     all_lines = ["OFF\n"] + all_lines;
        #     in_file.close();
        #     in_file = open(i, "w");
        #     in_file.writelines(all_lines);
        #     in_file.close();
        # else:
        #     in_file.close();
        # print(colored(">> Files fixed!", 'green', attrs = ['bold']))
        # fix_flag_train = True
#################### File Fixing in Train Folders ####################


#################### Normalization Scaling in Train Folders ####################
        # parsed_file = read_off(i)
        # sc = NormalizeScale()
        # scaled_pos = sc(parsed_file)
        # write_off(scaled_pos, i)
        # print(colored(">>> Rescaling done!\n", 'blue', attrs = ['bold']))
        # counter_after1 += 1 
        # print(str(counter_after1) + " out of " + str(counter1))
        # sleep(0.2)
        # os.system('clear')
        # scale_flag_train = True
        out_string_train = str('python /beegfs/gpfs0/a.morozov/EureconAugmentedDatasets/MN/Non-scaled/Eurepackage/eurecon_rmsd_relmin/examples/eurerun.py -r ' + str(rmsd) + ' -p 1 -a /beegfs/gpfs0/a.morozov/EureconAugmentedDatasets/MN/Non-scaled/Eurepackage/eurecon_rmsd_relmin/test_data/tess/tesselation_vertices_layer_0.txt -rr 1 -rrr 0.2 -i ' + str(i) + ' -o ' + str(object_path1) + '\n')
        # out_string_train = str('python /home/apm/Desktop/All_Eurecon/Eurepoint/Eurecon/eurecon/examples/eurerun.py -r ' + str(rmsd) + ' -p 1 -a /home/apm/Desktop/All_Eurecon/Eurepoint/Eurecon/eurecon/test_data/tess/tesselation_vertices_layer_0.txt -rr 1 -rrr 1 -i ' + str(i) + ' -o ' + str(object_path1) + '\n')
        file.write(out_string_train)
#################### Normalization Scaling in Train Folders ####################

counter_after2 = 0
for object_path2 in result2:
    files2 = [f for f in glob.glob(object_path2 + "/*.off", recursive=True)]
    counter2 += len(files2)

#################### File Fixing in Test Folders #####################
    
    for j in tqdm(files2, desc = 'Progress (Test Files)'):
        # print(colored("> Reading file: ", 'red', attrs = ['bold']) + j);
        # in_file = open(j, "r");
        # all_lines = in_file.readlines()
        # first_line = all_lines[0];
        # tokens = first_line.split();
        # # should have only one token: OFF
        # if(len(tokens) > 1):
        #     # Get the value with the 'OFF'
        #     tokens[0] = tokens[0].split("OFF")[1];
        #     edited_line = tokens[0] + " " + tokens[1] + " " + tokens[2] + "\n";
        #     all_lines[0] = edited_line;
        #     all_lines = ["OFF\n"] + all_lines;
        #     in_file.close();
        #     in_file = open(j, "w");
        #     in_file.writelines(all_lines);
        #     in_file.close();
        # else:
        #     in_file.close();
        # print(colored(">> Errors fixed!", 'green', attrs = ['bold']))
        # fix_flag_test = True

#################### File Fixing in Test Folders ####################


#################### Normalization Scaling in Test Folders ####################
        # parsed_file = read_off(j)
        # sc = NormalizeScale()
        # scaled_pos = sc(parsed_file)
        # write_off(scaled_pos, j)
        # print(colored(">>> Rescaling done!\n", 'green', attrs = ['bold']))
        # counter_after2 += 1 
        # print(str(counter_after2) + " out of " + str(counter2))
        # sleep(0.2)
        # os.system('clear')
        # scale_flag_test = True
        out_string_test = str('python /beegfs/gpfs0/a.morozov/EureconAugmentedDatasets/MN/Non-scaled/Eurepackage/eurecon_rmsd_relmin/examples/eurerun.py -r ' + str(rmsd) + ' -p 1 -a /beegfs/gpfs0/a.morozov/EureconAugmentedDatasets/MN/Non-scaled/Eurepackage/eurecon_rmsd_relmin/test_data/tess/tesselation_vertices_layer_0.txt -rr 1 -rrr 0.2 -i ' + str(i) + ' -o ' + str(object_path1) + '\n')
        file.write(out_string_test)
#################### Normalization Scaling in Test Folders ####################

# if fix_flag_train & fix_flag_test & scale_flag_train & scale_flag_test:
#     os.rename(input_directory, input_directory + 'SCRF')
            		
            
