#!/usr/bin/python
# AnGST


######################## add by zdh
import os
import time
import sys 

current_path = '/data2/zdh/software/angst-master/angst_lib/'
tree_lib_dir = os.path.join(current_path, '..', 'tree_lib')
tree_lib_dir = os.path.normpath(tree_lib_dir)
sys.path.append(tree_lib_dir)
sys.path.append(current_path)
######################## add by zdh end


# python libraries
import sys
import time
import pdb
import PyVM

sys.setrecursionlimit(1000000)
from AnGSTHelper import RunAnGST
from AnGSTInput import input_obj


# initiate variable for measuring running time and memory consumtion
start_time = time.time()
mem_str = []
mem_str.append(PyVM.MemoryUpdate("init",'return'))

# load inputs #
print "* read input"
input_info = input_obj(sys.argv[1])
input_dict = {}
input_dict['input_info'] = input_info
input_dict['mem_str'] = mem_str
input_dict['write_out'] = False
angst_inputs = [input_dict]

# run angst one more time at that scaling
input_dict['write_out'] = True
input_dict['start_time'] = start_time
RunAnGST(input_dict)


