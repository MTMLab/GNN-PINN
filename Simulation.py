#!/usr/bin/env python3

from rdkit import Chem
from rdkit import RDLogger
import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

# Hide RDKit warnings
RDLogger.DisableLog('rdApp.*')

def get_polymer_configuration():
    """
    Get polymer configuration from user input (similar to make_cpp.py)
    
    Returns:
        dict: Polymer configuration dictionary
    """
    print("\n=========================================")
    print("    Polymer Configuration")
    print("=========================================")
    
    # Select polymer type
    while True:
        polymer_type = input("Select polymer type (binary/ternary): ").strip().lower()
        if polymer_type in ["binary", "ternary"]:
            break
        print("Invalid input. Please enter 'binary' or 'ternary'.")
    
    # Input PPTA count
    while True:
        try:
            ppta_count = int(input("Enter the number of PPTA[TPC+PPD]: ").strip())
            if ppta_count <= 0:
                print("Number of PPTA must be positive.")
                continue
            break
        except ValueError:
            print("Enter a number.")
    
    # Input 34ODA count (for ternary)
    oda_count = None
    if polymer_type == "ternary":
        while True:
            try:
                oda_count = int(input("Enter the number of [TPC+34ODA]: ").strip())
                if oda_count < 0:
                    print("Number of [TPC+34ODA] must be zero or positive.")
                    continue
                break
            except ValueError:
                print("Enter a number.")
    else:
        oda_count = 0
    
    # Input cation count
    while True:
        try:
            cation_count = int(input("Enter the number of [DCA+PPD]: ").strip())
            if cation_count <= 0:
                print("Number of [DCA+PPD] must be positive.")
                continue
            break
        except ValueError:
            print("Enter a number.")
    
    return {
        'polymer_type': polymer_type,
        'ppta_count': ppta_count,
        'oda_count': oda_count if oda_count is not None else 0,
        'cation_count': cation_count
    }

def save_polymer_config(config):
    """
    Save polymer configuration to JSON file
    
    Args:
        config (dict): Polymer configuration dictionary
        
    Returns:
        bool: True if successful, False otherwise
    """
    config_path = os.path.join(os.getcwd(), "set", "polymer_config.json")
    
    try:
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, 'w', encoding='utf-8') as f:
            json.dump(config, f, indent=2, ensure_ascii=False)
        print(f"Polymer configuration saved to {config_path}")
        return True
    except Exception as e:
        print(f"Error saving polymer configuration: {e}")
        return False

def load_polymer_config():
    """
    Load polymer configuration from JSON file
    
    Returns:
        dict or None: Polymer configuration dictionary, or None if file not found
    """
    config_path = os.path.join(os.getcwd(), "set", "polymer_config.json")
    
    try:
        if os.path.exists(config_path):
            with open(config_path, 'r', encoding='utf-8') as f:
                config = json.load(f)
            return config
        else:
            print(f"Warning: Polymer configuration file not found at {config_path}")
            return None
    except Exception as e:
        print(f"Error reading polymer configuration: {e}")
        return None

def format_polymer_info(config):
    """
    Format polymer configuration into a readable string for result.txt
    
    Args:
        config (dict): Polymer configuration
        
    Returns:
        str: Formatted polymer information string
    """
    if config is None:
        return "N/A N/A N/A N/A"
    
    polymer_type = config.get('polymer_type', 'N/A')
    ppta_count = config.get('ppta_count', 'N/A')
    oda_count = config.get('oda_count', 0) if config.get('oda_count') is not None else 0
    cation_count = config.get('cation_count', 'N/A')
    
    return f"{polymer_type} {ppta_count} {oda_count} {cation_count}"

def check_existing_result(monomer_smiles, solvent_smiles, polymer_info):
    """
    Check if the simulation has already been done by looking in result.txt
    Compares canonical SMILES to handle different representations of the same molecule
    
    Args:
        monomer_smiles (str): Canonical monomer SMILES
        solvent_smiles (str): Canonical solvent SMILES
        polymer_info (str): Formatted polymer configuration string
        
    Returns:
        dict or None: Dictionary with existing results if found, None otherwise
    """
    result_file_path = os.path.join(os.getcwd(), "result.txt")
    
    if not os.path.exists(result_file_path):
        return None
    
    try:
        with open(result_file_path, 'r') as f:
            lines = f.readlines()
        
        # Parse each line to find matching system
        for line_num, line in enumerate(lines, 1):
            parts = line.strip().split()
            if len(parts) >= 8:  # Minimum expected parts
                # Extract components from the line
                # Format: monomer_SMILES solvent_SMILES polymer_type ppta_count oda_count cation_count stretched_interE solution_interE
                saved_monomer = parts[0]
                saved_solvent = parts[1]
                saved_polymer_info = ' '.join(parts[2:6])  # polymer_type ppta_count oda_count cation_count
                saved_stretched = parts[6]
                saved_solution = parts[7] if len(parts) > 7 else "N/A"
                
                # Canonicalize the saved SMILES for comparison
                try:
                    saved_monomer_canonical = canonicalize_smiles(saved_monomer)
                    saved_solvent_canonical = canonicalize_smiles(saved_solvent)
                    
                    # Skip if canonicalization fails (invalid SMILES in result.txt)
                    if saved_monomer_canonical is None or saved_solvent_canonical is None:
                        print(f"Warning: Invalid SMILES found in result.txt line {line_num}, skipping...")
                        continue
                    
                    # Compare canonical forms
                    if (saved_monomer_canonical == monomer_smiles and 
                        saved_solvent_canonical == solvent_smiles and 
                        saved_polymer_info == polymer_info):
                        
                    #    print(f"\n📌 Match found! (Line {line_num} in result.txt)")
                    #    print(f"   Saved monomer: {saved_monomer} →  {saved_monomer_canonical}")
                    #    print(f"   Input monomer: {monomer_smiles}")
                    #    print(f"   Saved solvent: {saved_solvent} →  {saved_solvent_canonical}")
                    #    print(f"   Input solvent: {solvent_smiles}")
                        
                        return {
                            'monomer_smiles': saved_monomer_canonical,
                            'solvent_smiles': saved_solvent_canonical,
                            'polymer_info': saved_polymer_info,
                            'stretched_interE': saved_stretched,
                            'solution_interE': saved_solution,
                            'line': line.strip(),
                            'original_monomer': saved_monomer,
                            'original_solvent': saved_solvent
                        }
                        
                except Exception as e:
                    print(f"Warning: Error processing SMILES in line {line_num}: {e}")
                    continue
    
    except Exception as e:
        print(f"Error reading result.txt: {e}")
        return None
    
    return None

def display_existing_result(result_data):
    """
    Display the existing result in a formatted way
    
    Args:
        result_data (dict): Dictionary containing the existing result data
    """
    print("\n" + "="*60)
    print("🔍 FOUND EXISTING SIMULATION RESULT!")
    print("="*60)
    
    # Show both original and canonical SMILES if different
    if 'original_monomer' in result_data and result_data['original_monomer'] != result_data['monomer_smiles']:
        print(f"Monomer SMILES (saved): {result_data['original_monomer']}")
        print(f"Monomer SMILES (canonical): {result_data['monomer_smiles']}")
    else:
        print(f"Monomer SMILES: {result_data['monomer_smiles']}")
    
    if 'original_solvent' in result_data and result_data['original_solvent'] != result_data['solvent_smiles']:
        print(f"Solvent SMILES (saved): {result_data['original_solvent']}")
        print(f"Solvent SMILES (canonical): {result_data['solvent_smiles']}")
    else:
        print(f"Solvent SMILES: {result_data['solvent_smiles']}")
    
    print(f"Polymer Configuration: {result_data['polymer_info']}")
    print("-"*60)
    print(f"✓ Stretched interE: {result_data['stretched_interE']}")
    print(f"✓ Solution interE: {result_data['solution_interE']}")
    print("="*60)
    #print("ℹ️  Skipping simulation - using cached result from result.txt")
    #print("="*60 + "\n")
    
    # Check if plot exists and display path
    polymer_parts = result_data['polymer_info'].split()
    if len(polymer_parts) == 4:
        polymer_type, ppta_count, oda_count, cation_count = polymer_parts
        # Try to find plot with canonical SMILES first
        plot_filename = f"{result_data['monomer_smiles']}_{result_data['solvent_smiles']}_{polymer_type}_{ppta_count}_{oda_count}_{cation_count}.png"
        plot_path = os.path.join("./Result_plot/", plot_filename)
        
        # If not found, try with original SMILES
        if not os.path.exists(plot_path) and 'original_monomer' in result_data:
            plot_filename = f"{result_data['original_monomer']}_{result_data['original_solvent']}_{polymer_type}_{ppta_count}_{oda_count}_{cation_count}.png"
            plot_path = os.path.join("./Result_plot/", plot_filename)
        
        # If still not found, try with sanitized filename
        if not os.path.exists(plot_path):
            # sanitize_filename is already defined in this file, no import needed
            base_filename = f"{result_data['monomer_smiles']}_{result_data['solvent_smiles']}_{polymer_type}_{ppta_count}_{oda_count}_{cation_count}"
            sanitized_filename = sanitize_filename(base_filename)
            plot_filename = f"{sanitized_filename}.png"
            plot_path = os.path.join("./Result_plot/", plot_filename)
        
        if os.path.exists(plot_path):
            print(f"📊 Plot available at: {plot_path}")
        else:
            print("📊 No plot file found for this system")


def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)

def save_smiles(smiles, folder="set", filename="name.txt"):
    path = os.path.join(os.getcwd(), folder)
    os.makedirs(path, exist_ok=True)
    filepath = os.path.join(path, filename)
    try:
        with open(filepath, 'w') as f:
            f.write(smiles)
        print(f"Saved canonical SMILES to {filepath}")
        return filepath
    except Exception as e:
        print(f"Error saving SMILES: {e}")
        return None

def save_solvent_smiles(smiles, folder="Solution/solvent", filename="solvent.smi"):
    path = os.path.join(os.getcwd(), folder)
    os.makedirs(path, exist_ok=True)
    filepath = os.path.join(path, filename)
    try:
        with open(filepath, 'w') as f:
            f.write(smiles)
        print(f"Saved solvent SMILES to {filepath}")
    except Exception as e:
        print(f"Error saving solvent SMILES: {e}")

def run_make_cpp_with_config(polymer_config):
    """
    Run make_cpp.py with polymer configuration using subprocess to provide inputs
    
    Args:
        polymer_config (dict): Polymer configuration dictionary
        
    Returns:
        bool: True if successful, False otherwise
    """
    cpp_dir = os.path.join(os.getcwd(), "Util", "Util_monomer_modify_4_cpp_make")
    script_path = os.path.join(cpp_dir, "make_cpp.py")

    if not os.path.isfile(script_path):
        print(f"Error: make_cpp.py not found at {script_path}")
        return False

    try:
        # Prepare input string for make_cpp.py
        input_lines = [
            polymer_config['polymer_type'],
            str(polymer_config['ppta_count']),
        ]
        
        if polymer_config['polymer_type'] == 'ternary':
            input_lines.append(str(polymer_config['oda_count']))
        
        input_lines.append(str(polymer_config['cation_count']))
        input_lines.append('y')  # Confirm settings
        
        input_string = '\n'.join(input_lines) + '\n'
        
        # Run make_cpp.py with inputs
        process = subprocess.run(
            ["python3", "make_cpp.py"], 
            cwd=cpp_dir, 
            input=input_string,
            text=True,
            capture_output=True,
            check=True
        )
        
        print("make_cpp.py executed successfully with polymer configuration.")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"Error during make_cpp.py execution: {e}")
        if e.stdout:
            print(f"Output: {e.stdout}")
        if e.stderr:
            print(f"Error: {e.stderr}")
        return False

def run_simulation(name_file_path, monomer_smiles, solvent_smiles):
    set_dir = os.path.dirname(name_file_path)
    exe_path = os.path.abspath(os.path.join(set_dir, "../Util/Util_Polymer_combine/exe"))

    if not os.path.isfile(exe_path):
        print(f"Error: Executable not found at {exe_path}")
        return

    try:
        subprocess.run([exe_path], cwd=set_dir, check=True)
        print("Initial simulation executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during initial simulation execution: {e}")

    copy_structures_to_targets(set_dir, monomer_smiles, solvent_smiles)

def copy_structures_to_targets(set_dir, monomer_smiles, solvent_smiles):
    source = os.path.join(set_dir, "structures")
    if not os.path.exists(source):
        print(f"Source directory {source} does not exist. Skipping copy.")
        return

    base_dir = os.getcwd()
    targets = ["Stretched", "Solution"]
    for target in targets:
        target_struct_path = os.path.join(base_dir, target, "structures")
        name_file_path = os.path.join(base_dir, target, "name.txt")

        if os.path.exists(target_struct_path):
            try:
                shutil.rmtree(target_struct_path)
                print(f"Removed existing {target_struct_path}")
            except Exception as e:
                print(f"Failed to remove {target_struct_path}: {e}")
                continue

        try:
            shutil.copytree(source, target_struct_path)
            print(f"Copied {source} → {target_struct_path}")
        except Exception as e:
            print(f"Failed to copy to {target_struct_path}: {e}")
            continue

        try:
            structure_files = sorted(os.listdir(source))
            with open(name_file_path, 'w') as f:
                for fname in structure_files:
                    name_without_ext = os.path.splitext(fname)[0]
                    f.write(name_without_ext + '\n')
            print(f"Wrote structure base names to {name_file_path}")
        except Exception as e:
            print(f"Failed to write structure filenames: {e}")

    run_final_stretch(base_dir)
    run_final_solution(base_dir, monomer_smiles, solvent_smiles)

def run_final_stretch(base_dir):
    stretch_dir = os.path.join(base_dir, "Stretched")
    exe_path = os.path.abspath(os.path.join(stretch_dir, "../Util/Util_Polymer_run_4_final_strain/exe"))

    if not os.path.isfile(exe_path):
        print(f"Error: Final stretching executable not found at {exe_path}")
        return

    try:
        subprocess.run([exe_path], cwd=stretch_dir, check=True)
        print("Final stretching simulation executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during final stretching simulation: {e}")

def run_final_solution(base_dir, monomer_smiles, solvent_smiles):
    solution_dir = os.path.join(base_dir, "Solution")
    exe_path = os.path.abspath(os.path.join(solution_dir, "../Util/Util_Polymer_run_4_final_polymer_solvent/exe"))

    if not os.path.isfile(exe_path):
        print(f"Error: Final solution executable not found at {exe_path}")
        return

    try:
        subprocess.run([exe_path], cwd=solution_dir, check=True)
        print("Final solution simulation executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during final solution simulation: {e}")

    print_interaction_energies(base_dir, monomer_smiles, solvent_smiles)

# === Analysis Functions ===

def read_lammps_data(filename):
    """
    Read LAMMPS output file and return data
    
    Args:
        filename (str): Path to the file to read
    
    Returns:
        pandas.DataFrame: DataFrame containing TimeStep and Value columns
    """
    try:
        # Read data, skipping comment lines
        data = pd.read_csv(filename, sep=r'\s+', comment='#', 
                          names=['TimeStep', 'Value'])
        return data
    except FileNotFoundError:
        print(f"File not found: {filename}")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

def sanitize_filename(filename):
    """
    Sanitize filename by replacing invalid characters with safe alternatives
    
    Args:
        filename (str): Original filename
        
    Returns:
        str: Sanitized filename safe for filesystem
    """
    # Replace characters that are invalid in filenames
    replacements = {
        '/': '_',
        '\\': '_',
        ':': '_',
        '*': '_',
        '?': '_',
        '"': '_',
        '<': '_',
        '>': '_',
        '|': '_',
        '(': '[',
        ')': ']',
        '=': '-',
        '+': 'plus',
        '#': 'hash'
    }
    
    sanitized = filename
    for old_char, new_char in replacements.items():
        sanitized = sanitized.replace(old_char, new_char)
    
    # Remove any remaining problematic characters and limit length
    sanitized = ''.join(c for c in sanitized if c.isalnum() or c in '-_[].')
    
    # Limit filename length to avoid filesystem issues
    if len(sanitized) > 200:
        sanitized = sanitized[:200]
    
    return sanitized

def plot_combined_figure(stretched_data, solution_data, monomer_smiles=None, solvent_smiles=None, save_plot=True):
    """
    Display Stretched and Solution data in 2x1 subplots
    
    Args:
        stretched_data: dict containing timesteps_ps, inter_energy, density for stretched
        solution_data: dict containing timesteps_ps, inter_energy, density for solution
        monomer_smiles (str): Monomer SMILES string for filename
        solvent_smiles (str): Solvent SMILES string for filename
        save_plot (bool): Whether to save the plot to file
    """
    # Load polymer configuration for filename
    polymer_config = load_polymer_config()
    polymer_info = format_polymer_info(polymer_config)
    # Font settings
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.unicode_minus'] = False
    
    # Create 2x1 subplots (adjust height and spacing)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 14))
    
    # Adjust subplot spacing
    plt.subplots_adjust(hspace=0.4, top=0.92, bottom=0.08)
    
    # === Stretched plot (top) ===
    color1 = 'tab:blue'
    color2 = 'tab:red'
    
    # First subplot - Stretched
    ax1_twin = ax1.twinx()
    
    # Interaction energy (blue, left y-axis)
    line1 = ax1.plot(stretched_data['timesteps_ps'], stretched_data['inter_energy'], 
                     '-', color=color1, linewidth=2.5, alpha=0.8, 
                     label='Chain-Chain Interaction Energy')
    ax1.set_xlabel('Time [ps]', fontsize=13)
    ax1.set_ylabel('Chain-Chain interE [Kcal/mol·Å³]', color=color1, fontsize=12)
    ax1.tick_params(axis='y', labelcolor=color1, labelsize=10)
    ax1.tick_params(axis='x', labelsize=10)
    ax1.grid(True, alpha=0.3, linestyle='--')
    
    # Density (red, right y-axis)
    line2 = ax1_twin.plot(stretched_data['timesteps_ps'], stretched_data['density'], 
                          '-', color=color2, linewidth=2.5, alpha=0.8, 
                          label='Density')
    ax1_twin.set_ylabel('Density [g/cm³]', color=color2, fontsize=12)
    ax1_twin.tick_params(axis='y', labelcolor=color2, labelsize=10)
    
    # Title and legend - Stretched
    ax1.set_title('Stretched State', fontsize=15, fontweight='bold', pad=15)
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_twin.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)
    
    # === Solution plot (bottom) ===
    ax2_twin = ax2.twinx()
    
    # Interaction energy (blue, left y-axis)
    line3 = ax2.plot(solution_data['timesteps_ps'], solution_data['inter_energy'], 
                     '-', color=color1, linewidth=2.5, alpha=0.8, 
                     label='Chain-Chain Interaction Energy')
    ax2.set_xlabel('Time [ps]', fontsize=13)
    ax2.set_ylabel('Chain-Chain interE [Kcal/mol·Å³]', color=color1, fontsize=12)
    ax2.tick_params(axis='y', labelcolor=color1, labelsize=10)
    ax2.tick_params(axis='x', labelsize=10)
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    # Density (red, right y-axis)
    line4 = ax2_twin.plot(solution_data['timesteps_ps'], solution_data['density'], 
                          '-', color=color2, linewidth=2.5, alpha=0.8, 
                          label='Density')
    ax2_twin.set_ylabel('Density [g/cm³]', color=color2, fontsize=12)
    ax2_twin.tick_params(axis='y', labelcolor=color2, labelsize=10)
    
    # Title and legend - Solution
    ax2.set_title('Solution State', fontsize=15, fontweight='bold', pad=15)
    lines3, labels3 = ax2.get_legend_handles_labels()
    lines4, labels4 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines3 + lines4, labels3 + labels4, loc='upper right', fontsize=10)
    
    # Overall title
    fig.suptitle('LAMMPS Time Series: Stretched vs Solution States', 
                 fontsize=17, fontweight='bold')
    
    if save_plot:
        # Create Result_plot folder if it doesn't exist
        output_dir = './Result_plot/'
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate filename from SMILES strings and polymer config - same format as result.txt
        if monomer_smiles and solvent_smiles:
            # Use format: monomerSMILES_solvent_SMILES_polymer_type_ppta_count_oda_count_cation_count.png
            # Parse polymer_info which is in format "polymer_type ppta_count oda_count cation_count"
            polymer_parts = polymer_info.split()
            if len(polymer_parts) == 4:
                polymer_type, ppta_count, oda_count, cation_count = polymer_parts
                filename = f"{monomer_smiles}_{solvent_smiles}_{polymer_type}_{ppta_count}_{oda_count}_{cation_count}.png"
            else:
                # Fallback if polymer config is not available
                filename = f"{monomer_smiles}_{solvent_smiles}.png"
        else:
            # Fallback to default filename
            filename = 'stretched_vs_solution_combined.png'
        
        # Set file path
        output_path = os.path.join(output_dir, filename)
        
        try:
            plt.savefig(output_path, dpi=300, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            plt.show()  # Display the plot
            return output_path  # Return the saved file path
        except Exception as e:
            # If saving with original SMILES fails, fall back to sanitized version
            print(f"Warning: Could not save with original SMILES filename. Using sanitized version.")
            print(f"Error: {e}")
            
            # Sanitize filename as fallback
            base_filename = f"{monomer_smiles}_{solvent_smiles}"
            polymer_parts = polymer_info.split()
            if len(polymer_parts) == 4:
                polymer_type, ppta_count, oda_count, cation_count = polymer_parts
                base_filename = f"{monomer_smiles}_{solvent_smiles}_{polymer_type}_{ppta_count}_{oda_count}_{cation_count}"
            
            sanitized_filename = sanitize_filename(base_filename)
            filename = f"{sanitized_filename}.png"
            output_path = os.path.join(output_dir, filename)
            
            plt.savefig(output_path, dpi=300, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            plt.show()
            return output_path
    
    plt.show()
    return None

def load_and_process_data(data_dir_name):
    """
    Load and process LAMMPS data from specific directory
    
    Args:
        data_dir_name (str): Data directory name ("Stretched" or "Solution")
    
    Returns:
        dict: Dictionary containing processed data
    """
    data_dir = f"./{data_dir_name}/lammps/monomer_0/"
    output1_file = os.path.join(data_dir, "output1.txt")
    output5_file = os.path.join(data_dir, "output5.txt")
    
    # Read data
    inter_energy_data = read_lammps_data(output1_file)
    density_data = read_lammps_data(output5_file)
    
    if inter_energy_data is None or density_data is None:
        return None
    
    # Extract data and convert time units
    timesteps_raw = inter_energy_data['TimeStep'].values
    inter_energy = inter_energy_data['Value'].values
    density = density_data['Value'].values
    
    # Convert TimeStep to actual time (0.5 fs timestep, convert to ps)
    timesteps_ps = timesteps_raw * 0.5 / 1000  # ps units
    
    return {
        'timesteps_raw': timesteps_raw,
        'timesteps_ps': timesteps_ps,
        'inter_energy': inter_energy,
        'density': density,
        'name': data_dir_name
    }

def run_analysis(monomer_smiles=None, solvent_smiles=None):
    """
    Analyze simulation results and generate graphs
    
    Args:
        monomer_smiles (str): Monomer SMILES string for filename
        solvent_smiles (str): Solvent SMILES string for filename
        
    Returns:
        str or None: Path to saved plot file, or None if failed
    """
    # Load Stretched data
    stretched_data = load_and_process_data("Stretched")
    if stretched_data is None:
        return None
    
    # Load Solution data
    solution_data = load_and_process_data("Solution")
    if solution_data is None:
        return None
    
    # Generate combined graph
    try:
        plot_path = plot_combined_figure(stretched_data, solution_data, 
                                       monomer_smiles=monomer_smiles, 
                                       solvent_smiles=solvent_smiles, 
                                       save_plot=True)
        return plot_path
    except Exception as e:
        return None

def print_interaction_energies(base_dir, monomer_smiles, solvent_smiles):
    stretch_ehbond = os.path.join(base_dir, "Stretched", "E_h_bond.txt")
    solution_ehbond = os.path.join(base_dir, "Solution", "E_h_bond.txt")

    stretched_value = read_first_energy(stretch_ehbond)
    unstretched_value = read_first_energy(solution_ehbond)

    if stretched_value is not None:
        print(f"stretched interE: {stretched_value}")
    else:
        print("stretched interE: N/A")

    if unstretched_value is not None:
        print(f"Solution interE: {unstretched_value}")
    else:
        print("Solution interE: N/A")

    # Load polymer configuration
    polymer_config = load_polymer_config()
    polymer_info = format_polymer_info(polymer_config)

    # Save result with polymer configuration
    stretched_str = f"{stretched_value:.6f}" if stretched_value is not None else "N/A"
    solution_str = f"{unstretched_value:.6f}" if unstretched_value is not None else "N/A"

    result_line = f"{monomer_smiles} {solvent_smiles} {polymer_info} {stretched_str} {solution_str}"

    result_file_path = os.path.join(os.getcwd(), "result.txt")
    result_relative_path = "./result.txt"
    try:
        with open(result_file_path, "a") as f:
            f.write(result_line + "\n")
        print(f"Results saved to {result_relative_path}")
    except Exception as e:
        print(f"Error writing to result.txt: {e}")
    
    # Run analysis after saving results
    plot_path = run_analysis(monomer_smiles=monomer_smiles, 
                           solvent_smiles=solvent_smiles)
    
    if plot_path:
        print(f"Results plots saved to {plot_path}")
    else:
        print("Error: Could not generate analysis plots.")

def read_first_energy(filepath):
    try:
        with open(filepath, 'r') as f:
            line = f.readline()
            parts = line.strip().split()
            if len(parts) >= 2:
                return float(parts[1])
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return None

def main():
    while True:
        smiles_input = input("Enter a dicarboxylic acid(DCA) monomer SMILES string (or 'q' to quit): ")
        if smiles_input.lower() == 'q':
            break

        canonical = canonicalize_smiles(smiles_input)
        if canonical is None:
            print("Invalid monomer SMILES. Please check the input.")
            continue

        solvent_input = input("Enter a solvent SMILES string: ")
        solvent_canonical = canonicalize_smiles(solvent_input)
        if solvent_canonical is None:
            print("Invalid solvent SMILES. Please check the input.")
            continue

        # Get polymer configuration from user BEFORE checking for duplicates
        polymer_config = get_polymer_configuration()
        
        # Format polymer info for comparison
        polymer_info = format_polymer_info(polymer_config)
        
        # Display configuration summary
        print(f"\n=== Configuration Summary ===")
        print(f"Monomer SMILES: {canonical}")
        print(f"Solvent SMILES: {solvent_canonical}")
        print(f"Polymer Type: {polymer_config['polymer_type']}")
        print(f"PPTA Count: {polymer_config['ppta_count']}")
        print(f"ODA Count: {polymer_config['oda_count']}")
        print(f"Cation Count: {polymer_config['cation_count']}")
        print("="*30)
        
        # Check if this system has already been simulated
        existing_result = check_existing_result(canonical, solvent_canonical, polymer_info)
        
        if existing_result:
            # Display the existing result
            display_existing_result(existing_result)
            
            # Ask user if they want to re-run the simulation anyway
            user_choice = input("\nDo you want to re-run the simulation anyway? (y/n): ").lower()
            if user_choice != 'y':
                continue  # Skip to next iteration
            else:
                print("\n⚙️  Re-running simulation as requested...")
        
        # Proceed with simulation (either no existing result or user chose to re-run)
        # Save SMILES files
        name_file_path = save_smiles(canonical)
        if not name_file_path:
            continue

        save_solvent_smiles(solvent_canonical)
        
        # Save polymer configuration to JSON
        if not save_polymer_config(polymer_config):
            print("Failed to save polymer configuration. Continuing anyway...")
        
        # Run make_cpp with the polymer configuration
        if not run_make_cpp_with_config(polymer_config):
            print("Failed to run make_cpp.py. Aborting simulation.")
            continue
        
        # Run the simulation
        run_simulation(name_file_path, canonical, solvent_canonical)

if __name__ == "__main__":
    main()
