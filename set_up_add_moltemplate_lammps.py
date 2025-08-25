#!/usr/bin/env python3

import os
import subprocess

def run_command(command, cwd=None, description=""):
    """Execute a command and handle errors"""
    try:
        print(f"Running: {description if description else command}")
        subprocess.run(command, shell=True, check=True, cwd=cwd)
        print(f"✓ Success: {description if description else command}")
    except subprocess.CalledProcessError as e:
        print(f"✗ Error in {description if description else command}: {e}")
        return False
    return True

def main():
    print("Starting setup process...")
    
    # Set the Util path based on the current directory
    base_util_path = os.path.join(os.getcwd(), "Util")
    
    # Step 1: Compile Util folders
    print("\n=== Step 1: Compiling Util folders ===")
    for folder_name in os.listdir(base_util_path):
        folder_path = os.path.join(base_util_path, folder_name)

        # It is a directory, its name contains "Util," and it is not "Util_monomer_modify_4."
        if (os.path.isdir(folder_path) and 
            "Util" in folder_name and 
            "Util_monomer_modify_4" not in folder_name):
            
            print(f"\nCompiling: {folder_name}")
            
            # Clean previous builds
            if not run_command("rm -rf *.o", cwd=folder_path, description=f"Cleaning {folder_name}"):
                continue
                
            # Compile
            if not run_command("make", cwd=folder_path, description=f"Making {folder_name}"):
                continue
    
    # Step 2: Install moltemplate
    print("\n=== Step 2: Installing moltemplate ===")
    moltemplate_path = os.path.join(base_util_path, "moltemplate-master")
    
    if os.path.exists(moltemplate_path):
        print(f"Found moltemplate directory: {moltemplate_path}")
        run_command("pip3 install .", cwd=moltemplate_path, description="Installing moltemplate")
    else:
        print(f"Warning: moltemplate directory not found at {moltemplate_path}")
    
    # Step 3: Compile LAMMPS
    print("\n=== Step 3: Compiling LAMMPS ===")
    lammps_src_path = os.path.join(base_util_path, "lammps-2Aug2023", "src")
    
    if os.path.exists(lammps_src_path):
        print(f"Found LAMMPS src directory: {lammps_src_path}")
        
        # Compile serial version
        print("\nCompiling LAMMPS serial version...")
        run_command("make serial", cwd=lammps_src_path, description="Making LAMMPS serial")
        
        # Compile MPI version
        print("\nCompiling LAMMPS MPI version...")
        run_command("make mpi", cwd=lammps_src_path, description="Making LAMMPS MPI")
        
    else:
        print(f"Warning: LAMMPS src directory not found at {lammps_src_path}")
    
    print("\n=== Setup completed ===")

if __name__ == "__main__":
    main()
