#!/usr/bin/env python3

import os
import subprocess

# Set the Util path based on the current directory.
base_util_path = os.path.join(os.getcwd(), "Util")

# Check all subitems in the Util folder
for folder_name in os.listdir(base_util_path):
    folder_path = os.path.join(base_util_path, folder_name)

    # It is a directory, its name contains “Util,” and it is not “Util_monomer_modify_4.”
    if (os.path.isdir(folder_path) and 
        "Util" in folder_name and 
        "Util_monomer_modify_4" not in folder_name):
        
        print(f"Compiling: {folder_name}")
        try:
            # Go to the directory and execute the compilation command.
            subprocess.run("rm -rf *.o", shell=True, check=True, cwd=folder_path)
            subprocess.run("make", shell=True, check=True, cwd=folder_path)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred in {folder_name}: {e}")
