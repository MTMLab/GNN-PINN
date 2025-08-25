#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import random
import re
import sys
import subprocess
import shutil
import json

class IntegratedPolymerGenerator:
    def __init__(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.template_file = os.path.join(script_dir, "MAIN_template.cpp")
        
    def read_template(self):
        """Read template file"""
        try:
            if not os.path.exists(self.template_file):
                print(f"Error: '{self.template_file}' file does not exist in the current directory.")
                print(f"Current directory: {os.getcwd()}")
                print(f"Files in directory: {os.listdir('.')}")
                sys.exit(1)
                
            with open(self.template_file, 'r', encoding='utf-8') as f:
                content = f.read()
                return content
        except Exception as e:
            print(f"Template file reading error: {str(e)}")
            sys.exit(1)
    
    def write_output(self, content, filename):
        """Write output file"""
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                f.write(content)
                
            if os.path.exists(filename):
                print(f"{filename} file has been created.")
            else:
                print(f"Warning: {filename} file was created but cannot be verified.")
        except Exception as e:
            print(f"Output file writing error: {str(e)}")
            sys.exit(1)

class LinearPolymerGenerator(IntegratedPolymerGenerator):
    """Linear polymer structure generator"""
    
    def __init__(self):
        super().__init__()
        self.output_file = "MAIN_linear.cpp"
        
    def generate_random_positions(self, ppta_count, oda_count=0, cation_count=0):
        """Generate random monomer positions (linear structure)"""
        total_internal_monomers = ppta_count + oda_count + cation_count
        
        positions = list(range(total_internal_monomers))
        random.shuffle(positions)
        
        ppta_positions = sorted(positions[:ppta_count])
        
        if oda_count > 0:
            oda_positions = sorted(positions[ppta_count:ppta_count+oda_count])
            cation_positions = sorted(positions[ppta_count+oda_count:])
        else:
            oda_positions = []
            cation_positions = sorted(positions[ppta_count:])
        
        # Actual indices start from 1 (0 is H_head, last is H_tail)
        ppta_positions = [p + 1 for p in ppta_positions]
        oda_positions = [p + 1 for p in oda_positions]
        cation_positions = [p + 1 for p in cation_positions]
        
        return {
            'ppta': ppta_positions,
            'oda': oda_positions,
            'cation': cation_positions
        }
    
    def generate_import_statements(self, use_oda=False):
        """Generate import statements"""
        if use_oda:
            return (
                'outputFile40 << "import monomer_add.lt" << endl;\n'
                'outputFile40 << "import H_head.lt" << endl;\n'
                'outputFile40 << "import H_tail.lt" << endl;\n'
                'outputFile40 << "import 34ODA.lt" << endl;\n'
                'outputFile40 << "import PPTA.lt" << endl;\n'
            )
        else:
            return (
                'outputFile40 << "import monomer_add.lt" << endl;\n'
                'outputFile40 << "import H_head.lt" << endl;\n'
                'outputFile40 << "import H_tail.lt" << endl;\n'
                'outputFile40 << "import PPTA.lt" << endl;\n'
            )
    
    def generate_polymer_code(self, total_monomers, positions, use_oda=False):
        """Generate linear polymer code"""
        code = []
        
        # H_head (index 0)
        code.append('outputFile40 << "monomers[0] = new H_head" << endl;')
        code.append('d0 = 11;')
        
        # Internal monomers (index 1 ~ total_monomers-2)
        for i in range(1, total_monomers - 1):
            if i in positions['ppta']:
                monomer_type = 'PPTA'
                next_d0 = 'd0 + interval'
            elif use_oda and i in positions['oda']:
                monomer_type = '34ODA'
                next_d0 = 'd0 + 39.2185'
            else:
                monomer_type = 'cation'
                next_d0 = 'd0 + 11'
            
            code.append(f'outputFile40 << "monomers[{i}] = new {monomer_type}.move(" << d0 << ",0,0)" << endl;  // {monomer_type}')
            code.append(f'd0 = {next_d0};')
        
        # H_tail (last index)
        last_idx = total_monomers - 1
        code.append(f'outputFile40 << "monomers[{last_idx}] = new H_tail.move(" << d0 << ",0,0)" << endl;')
        code.append('d0 = d0 + interval;')
        
        return '\n'.join(code)
    
    def generate_bond_code(self, total_monomers):
        """Generate linear bond code"""
        code = []
        
        for i in range(total_monomers - 1):
            bond_idx = i + 1
            current_idx = i
            next_idx = i + 1
            
            # First bond (H_head and first internal monomer)
            if i == 0:
                code.append(f'outputFile40 << "    $bond:b{bond_idx}  $atom:monomers[{current_idx}]/H_head $atom:monomers[{next_idx}]/atom1" << endl;')
            # Last bond (last internal monomer and H_tail)
            elif i == total_monomers - 2:
                code.append(f'outputFile40 << "    $bond:b{bond_idx}  $atom:monomers[{current_idx}]/atom2 $atom:monomers[{next_idx}]/H_tail" << endl;')
            # Internal bonds
            else:
                code.append(f'outputFile40 << "    $bond:b{bond_idx}  $atom:monomers[{current_idx}]/atom2 $atom:monomers[{next_idx}]/atom1" << endl;')
        
        return '\n'.join(code)
    
    def create_cpp_file(self, ppta_count=5, oda_count=None, cation_count=5):
        """Create linear structure C++ file"""
        try:
            # Set default values and check types
            if oda_count is not None:
                total_monomers = ppta_count + oda_count + cation_count + 2
                use_oda = True
            else:
                oda_count = 0
                total_monomers = ppta_count + cation_count + 2
                use_oda = False
            
            # Generate random positions
            positions = self.generate_random_positions(ppta_count, oda_count, cation_count)
            
            # Read template code
            template = self.read_template()
            
            # Check marker existence
            if "// Polymer chain generation code start" not in template or "// Polymer chain generation code end" not in template:
                print("Error: Template file does not contain polymer chain generation code markers.")
                return False
                
            if "// Bond structure definition start" not in template or "// Bond structure definition end" not in template:
                print("Error: Template file does not contain bond structure definition markers.")
                return False
            
            # Modify import statements
            import_pattern = r'outputFile40 << "import monomer_add\.lt" << endl;\s*' \
                             r'outputFile40 << "import 34ODA\.lt" << endl;\s*' \
                             r'outputFile40 << "import PPTA\.lt" << endl;'
            alt_import_pattern = r'outputFile40 << "import monomer_add\.lt" << endl;\s*' \
                                 r'outputFile40 << "import PPTA\.lt" << endl;'
            
            import_statements = self.generate_import_statements(use_oda)
            
            # Replace import statements
            if re.search(import_pattern, template):
                template = re.sub(import_pattern, import_statements.strip(), template)
            elif re.search(alt_import_pattern, template):
                template = re.sub(alt_import_pattern, import_statements.strip(), template)
            
            # Generate monomer code
            polymer_code = self.generate_polymer_code(total_monomers, positions, use_oda)
            
            # Generate bond code
            bond_code = self.generate_bond_code(total_monomers)
            
            # Find and replace markers
            pattern = r"// Polymer chain generation code start(.*?)// Polymer chain generation code end"
            replacement = f"// Polymer chain generation code start\n{polymer_code}\n// Polymer chain generation code end"
            content = re.sub(pattern, replacement, template, flags=re.DOTALL)
            
            pattern = r"// Bond structure definition start(.*?)// Bond structure definition end"
            replacement = f"// Bond structure definition start\noutputFile40 << \"   write('Data Bond List') {{\" << endl;\n{bond_code}\n// Bond structure definition end"
            content = re.sub(pattern, replacement, content, flags=re.DOTALL)
            
            # Modify system_new.lt file boundary settings
            system_pattern = r'outputFile401 << "0\.0    " << d0 << "  xlo xhi" << endl;'
            system_replacement = 'outputFile401 << "0.0    " << d0+2 << "  xlo xhi" << endl;'
            content = re.sub(system_pattern, system_replacement, content)
            
            # Write output file
            self.write_output(content, self.output_file)
            return True
                
        except Exception as e:
            print(f"Error occurred during linear structure C++ file generation: {str(e)}")
            return False

class FiberPolymerGenerator(IntegratedPolymerGenerator):
    """Fiber structure polymer generator"""
    
    def __init__(self):
        super().__init__()
        self.output_file = "MAIN_fiber.cpp"
        
    def generate_random_positions(self, polymer_count, ppta_count, oda_count=0, cation_count=0):
        """Generate random monomer positions (last monomer fixed as PPTA)"""
        total_monomers = ppta_count + oda_count + cation_count
        
        polymer_positions = {}
        
        for polymer_idx in range(polymer_count):
            adjusted_ppta_count = ppta_count - 1
            
            positions = list(range(total_monomers - 1))
            random.shuffle(positions)
            
            ppta_positions = sorted(positions[:adjusted_ppta_count])
            
            if oda_count > 0:
                oda_positions = sorted(positions[adjusted_ppta_count:adjusted_ppta_count+oda_count])
                cation_positions = sorted(positions[adjusted_ppta_count+oda_count:])
            else:
                oda_positions = []
                cation_positions = sorted(positions[adjusted_ppta_count:])
            
            # Add last monomer position as PPTA
            ppta_positions.append(total_monomers - 1)
            
            # Calculate actual indices
            base_idx = polymer_idx * total_monomers
            ppta_positions = [p + base_idx for p in ppta_positions]
            oda_positions = [p + base_idx for p in oda_positions]
            cation_positions = [p + base_idx for p in cation_positions]
            
            polymer_positions[polymer_idx] = {
                'ppta': ppta_positions,
                'oda': oda_positions,
                'cation': cation_positions
            }
        
        return polymer_positions
    
    def generate_import_statements(self, use_oda=False):
        """Generate import statements"""
        if use_oda:
            return (
                'outputFile40 << "import monomer_add.lt" << endl;\n'
                'outputFile40 << "import 34ODA.lt" << endl;\n'
                'outputFile40 << "import PPTA.lt" << endl;\n'
            )
        else:
            return (
                'outputFile40 << "import monomer_add.lt" << endl;\n'
                'outputFile40 << "import PPTA.lt" << endl;\n'
            )
    
    def generate_polymer_code(self, polymer_count, total_monomers, positions, use_oda=False):
        """Generate polymer code"""
        code = []
        polymer_y_positions = [0, 40, 80, 120, 160]
        
        for polymer_idx in range(polymer_count):
            y = polymer_y_positions[polymer_idx]
            base_idx = polymer_idx * total_monomers
            first_idx = base_idx
            
            # Determine first monomer type
            first_type = None
            if first_idx in positions[polymer_idx]['ppta']:
                first_type = 'PPTA'
            elif use_oda and first_idx in positions[polymer_idx]['oda']:
                first_type = '34ODA'
            else:
                first_type = 'cation'
            
            # Add first monomer code
            code.append(f'outputFile40 << "monomers[{first_idx}] = new {first_type}.move(0,{y},0)" << endl;')
            
            if first_type == 'PPTA':
                code.append(f'd0 = 12;')
            elif first_type == '34ODA':
                code.append(f'd0 = 39.2185;')
            else:
                code.append(f'd0 = interval;')
            
            # Add remaining monomers
            for i in range(1, total_monomers):
                current_idx = base_idx + i
                
                monomer_type = None
                if current_idx in positions[polymer_idx]['ppta']:
                    monomer_type = 'PPTA'
                    next_d0 = 'd0 + 12'
                elif use_oda and current_idx in positions[polymer_idx]['oda']:
                    monomer_type = '34ODA'
                    next_d0 = 'd0 + 39.2185'
                else:
                    monomer_type = 'cation'
                    next_d0 = 'd0 + interval'
                
                code.append(f'outputFile40 << "monomers[{current_idx}] = new {monomer_type}.move(" << d0 << ",{y},0)" << endl;  // {monomer_type}')
                code.append(f'd0 = {next_d0};')
        
        return '\n'.join(code)
    
    def generate_bond_code(self, polymer_count, total_monomers):
        """Generate bond code"""
        code = []
        
        for polymer_idx in range(polymer_count):
            base_idx = polymer_idx * total_monomers
            
            # Add comments for each polymer chain
            chain_names = ["First", "Second", "Third", "Fourth", "Fifth"]
            if polymer_idx < len(chain_names):
                max_idx = base_idx + total_monomers - 1
                code.append(f"// {chain_names[polymer_idx]} chain bonds ({base_idx}-{max_idx})")
            
            # Define bonds between each monomer
            for i in range(total_monomers):
                current_idx = base_idx + i
                next_idx = base_idx + ((i + 1) % total_monomers)
                bond_idx = base_idx + i + 1
                
                code.append(f'outputFile40 << "    $bond:b{bond_idx}  $atom:monomers[{current_idx}]/atom2 $atom:monomers[{next_idx}]/atom1" << endl;')
        
        return '\n'.join(code)
    
    def create_cpp_file(self, polymer_count=5, ppta_count=None, oda_count=None, cation_count=None):
        """Create fiber structure C++ file"""
        try:
            # Set default values
            if ppta_count is None:
                ppta_count = 10
            
            # Calculate total monomers
            if oda_count is not None:
                total_monomers = ppta_count + oda_count + cation_count
                use_oda = True
            else:
                oda_count = 0
                total_monomers = ppta_count + cation_count
                use_oda = False
            
            # Generate random positions
            positions = self.generate_random_positions(polymer_count, ppta_count, oda_count, cation_count)
            
            # Read template code
            template = self.read_template()
            
            # Check marker existence
            if "// Polymer chain generation code start" not in template or "// Polymer chain generation code end" not in template:
                print("Error: Template file does not contain polymer chain generation code markers.")
                return False
                
            if "// Bond structure definition start" not in template or "// Bond structure definition end" not in template:
                print("Error: Template file does not contain bond structure definition markers.")
                return False
            
            # Modify import statements
            import_pattern = r'outputFile40 << "import monomer_add\.lt" << endl;\s*' \
                             r'outputFile40 << "import 34ODA\.lt" << endl;\s*' \
                             r'outputFile40 << "import PPTA\.lt" << endl;'
            alt_import_pattern = r'outputFile40 << "import monomer_add\.lt" << endl;\s*' \
                                 r'outputFile40 << "import PPTA\.lt" << endl;'
            
            import_statements = self.generate_import_statements(use_oda)
            
            # Replace import statements
            if re.search(import_pattern, template):
                template = re.sub(import_pattern, import_statements.strip(), template)
            elif re.search(alt_import_pattern, template):
                template = re.sub(alt_import_pattern, import_statements.strip(), template)
            
            # Generate monomer code
            polymer_code = self.generate_polymer_code(polymer_count, total_monomers, positions, use_oda)
            
            # Generate bond code
            bond_code = self.generate_bond_code(polymer_count, total_monomers)
            
            # Find and replace markers
            pattern = r"// Polymer chain generation code start(.*?)// Polymer chain generation code end"
            replacement = f"// Polymer chain generation code start\n{polymer_code}\n// Polymer chain generation code end"
            content = re.sub(pattern, replacement, template, flags=re.DOTALL)
            
            pattern = r"// Bond structure definition start(.*?)// Bond structure definition end"
            replacement = f"// Bond structure definition start\noutputFile40 << \"   write('Data Bond List') {{\" << endl;\n{bond_code}\n// Bond structure definition end"
            content = re.sub(pattern, replacement, content, flags=re.DOTALL)
            
            # Write output file
            self.write_output(content, self.output_file)
            return True
                
        except Exception as e:
            print(f"Error occurred during fiber structure C++ file generation: {str(e)}")
            return False

def compile_linear_cpp(filename):
    """Move linear structure C++ file to different location and compile"""
    try:
        target_dir = "../Util_monomer_modify_4_polymer_end_H_make"
        current_dir = os.getcwd()
        
        # Check target directory existence
        if not os.path.exists(target_dir):
            print(f"Error: Target directory '{target_dir}' does not exist.")
            return False
        
        target_path = os.path.join(target_dir, "MAIN.cpp")
               
        # Copy file to target location and rename to MAIN.cpp
        shutil.copy(filename, target_path)
        print(f"Moved {filename} to {target_path}.")
        
        # Move to target directory and compile
        os.chdir(target_dir)
        print(f"\nCompiling in {target_dir}...")
        
        compile_process = subprocess.run("make", shell=True, capture_output=True, text=True)
        
        if compile_process.returncode != 0:
            print(f"Compilation error occurred:")
            print(f"Error message:\n{compile_process.stderr}")
            if compile_process.stdout:
                print(f"Standard output:\n{compile_process.stdout}")
            return False
        else:
            print(f"Linear structure compilation successful!")
            
            # Copy compiled executable to original location (optional)
            if os.path.exists("main"):
                original_main_path = os.path.join(current_dir, "main_linear")
                shutil.copy("main", original_main_path)
                print(f"Copied executable to {original_main_path}.")
            
            return True
            
    except Exception as e:
        print(f"Error occurred during linear structure compilation process: {str(e)}")
        return False
    finally:
        # Return to original directory
        os.chdir(current_dir)

def compile_fiber_cpp(filename):
    """Move fiber structure C++ file to different location and compile"""
    try:
        target_dir = "../Util_monomer_modify_4_new_make"
        current_dir = os.getcwd()
        
        # Check target directory existence
        if not os.path.exists(target_dir):
            print(f"Error: Target directory '{target_dir}' does not exist.")
            return False
        
        target_path = os.path.join(target_dir, "MAIN.cpp")
                
        # Copy file to target location and rename to MAIN.cpp
        shutil.copy(filename, target_path)
        print(f"Moved {filename} to {target_path}.")
        
        # Move to target directory and compile
        os.chdir(target_dir)
        print(f"\nCompiling in {target_dir}...")
        
        compile_process = subprocess.run("make", shell=True, capture_output=True, text=True)
        
        if compile_process.returncode != 0:
            print(f"Compilation error occurred:")
            print(f"Error message:\n{compile_process.stderr}")
            if compile_process.stdout:
                print(f"Standard output:\n{compile_process.stdout}")
            return False
        else:
            print(f"Fiber structure compilation successful!")
            
            # Copy compiled executable to original location (optional)
            if os.path.exists("main"):
                original_main_path = os.path.join(current_dir, "main_fiber")
                shutil.copy("main", original_main_path)
                print(f"Copied executable to {original_main_path}.")
            
            return True
            
    except Exception as e:
        print(f"Error occurred during fiber structure compilation process: {str(e)}")
        return False
    finally:
        # Return to original directory
        os.chdir(current_dir)

def save_polymer_config(config):
    """Save polymer configuration to JSON file"""
    try:
        # Save to the set directory where Simulation.py creates other files
        config_dir = "../../set"  # Go up two levels from Util/Util_monomer_modify_4_cpp_make/ and into set/
        os.makedirs(config_dir, exist_ok=True)  # Create set directory if it doesn't exist
        config_path = os.path.join(config_dir, "polymer_config.json")
        
        with open(config_path, 'w', encoding='utf-8') as f:
            json.dump(config, f, indent=2, ensure_ascii=False)
        
        print(f"Polymer configuration saved to {config_path}")
        return True
        
    except Exception as e:
        print(f"Error saving polymer configuration: {str(e)}")
        return False

def get_user_inputs():
    """Get user inputs and return configuration"""
    print("\n=========================================")
    print("    Polymer Structure Generator")
    print("=========================================")
    print("This program generates both fiber and single linear structures.")
    
    # Select polymer type
    while True:
        polymer_type = input("\nSelect polymer type (binary/ternary): ").strip().lower()
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
    
    # Input cation count
    while True:
        try:
            cation_count = int(input("Enter the number of [DCA+PPD]: ").strip())
            if cation_count <= 0:
                print("Number of [DCA+PPD] must be zero or positive.")
                continue
            break
        except ValueError:
            print("Enter a number.")
    
    return {
        'polymer_type': polymer_type,
        'ppta_count': ppta_count,
        'oda_count': oda_count,
        'cation_count': cation_count
    }

def run_integrated_interface():
    """Start Integrated Interactive Interface"""
    while True:  # Main loop for restart functionality
        # Get user inputs
        config = get_user_inputs()
        
        # Fixed polymer chain count for fiber structure
        polymer_count = 5
        
        # Print structure summary
        print(f"\n=== Structure Generation Summary ({config['polymer_type'].upper()}) ===")
        print(f"1. Fiber Structure :")
        print(f"   - Each polymer: PPTA({config['ppta_count']} units)")
        if config['polymer_type'] == "ternary":
            print(f"   + [TPC+34ODA]({config['oda_count']} units)")
        print(f"   + [DCA+PPD]({config['cation_count']} units)")
        
        print("\n2. Single Linear Polymer :")
        print(f"   - H_head → PPTA[TPC+PPD]({config['ppta_count']} units)")
        if config['polymer_type'] == "ternary":
            print(f"   + [TPC+34ODA]({config['oda_count']} units)")
        print(f"   + [DCA+PPD]({config['cation_count']} units) → H_tail")
        
        # Confirmation
        confirm = input("\nContinue with these settings? (y/n): ").strip().lower()
        if confirm == 'y':
            break  # Exit the loop and continue with generation
        elif confirm == 'n':
            print("\nRestarting configuration...")
            continue  # Go back to the beginning of the loop
        else:
            print("Invalid input. Please enter 'y' or 'n'.")
            # Ask again without restarting the whole loop
            while True:
                confirm = input("Continue with these settings? (y/n): ").strip().lower()
                if confirm in ['y', 'n']:
                    break
                print("Invalid input. Please enter 'y' or 'n'.")
            if confirm == 'y':
                break
            else:
                print("\nRestarting configuration...")
                continue
    
    # Save configuration to JSON file
    save_polymer_config(config)
    
    print("\nGenerating files...")

    # Generate fiber structure
    print("\n1. Generating fiber structure...")
    fiber_generator = FiberPolymerGenerator()
    if not fiber_generator.create_cpp_file(
        polymer_count=polymer_count,
        ppta_count=config['ppta_count'],
        oda_count=config['oda_count'],
        cation_count=config['cation_count']
    ):
        print("Failed to generate fiber structure")
        return
            
    # Generate linear structure
    print("\n2. Generating single structure...")
    linear_generator = LinearPolymerGenerator()
    if not linear_generator.create_cpp_file(
        ppta_count=config['ppta_count'],
        oda_count=config['oda_count'],
        cation_count=config['cation_count']
    ):
        print("Failed to generate single structure")
        return
    
    # Automatic compilation (no user input required)
    print("\n=== Compilation in progress ===")

    # Compile fiber structure (at different location)
    print("\n1. Compiling fiber structure...")
    fiber_success = compile_fiber_cpp("MAIN_fiber.cpp")
            
    # Compile linear structure (at different location)
    print("\n2. Compiling single structure...")
    linear_success = compile_linear_cpp("MAIN_linear.cpp")
    
    if linear_success and fiber_success:
        print("\nAll structure compilations completed!")
        print("Generated files:")
        print("  - MAIN_linear.cpp (linear structure source code)")
        print("  - MAIN_fiber.cpp (fiber structure source code)")
        print("  - main_linear (linear structure executable)")
        print("  - main_fiber (fiber structure executable)")
        print("  - ../Util_monomer_modify_4_polymer_end_H_make/MAIN.cpp (linear structure)")
        print("  - ../Util_monomer_modify_4_polymer_end_H_make/main (linear structure executable)")
        print("  - ../Util_monomer_modify_4_new_make/MAIN.cpp (fiber structure)")
        print("  - ../Util_monomer_modify_4_new_make/main (fiber structure executable)")
    elif linear_success:
        print("\nLinear structure compilation succeeded but fiber structure compilation failed.")
    elif fiber_success:
        print("\nFiber structure compilation succeeded but linear structure compilation failed.")
    
    print("\nCompleted!")

if __name__ == "__main__":
    run_integrated_interface()
