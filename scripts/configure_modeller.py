import os
import sys
import glob

def configure_modeller(license_key='MODELIRANJE'):
    print("Configuring Modeller License...")
    dest_config = None
    
    # 1. Try importing modeller to find path
    try:
        import modeller
        modeller_path = os.path.dirname(modeller.__file__)
        candidate = os.path.join(modeller_path, "config.py")
        if os.path.exists(candidate):
            dest_config = candidate
            print(f"Found config via import: {dest_config}")
    except Exception as e:
        print(f"Could not import modeller to find path: {e}")

    # 2. Fallback search
    if not dest_config:
        print("Searching for config.py in common locations...")
        possible_paths = [
            f"{sys.prefix}/lib/modeller-*/modlib/modeller/config.py",
            f"{sys.prefix}/lib/python*/site-packages/modeller/config.py",
            # Add specific conda env path if sys.prefix isn't enough
            "/usr/local/envs/FrankPEPstein/lib/modeller-*/modlib/modeller/config.py" 
        ]
        for pattern in possible_paths:
            found = glob.glob(pattern)
            if found:
                dest_config = found[0]
                print(f"Found config via search: {dest_config}")
                break

    if dest_config:
        try:
            with open(dest_config, 'r') as f:
                lines = f.readlines()
            
            new_lines = []
            changed = False
            for line in lines:
                if line.strip().startswith("license ="):
                    print(f"Found license line: {line.strip()}")
                    new_lines.append(f"license = '{license_key}'\n")
                    changed = True
                else:
                    new_lines.append(line)
            
            if changed:
                with open(dest_config, 'w') as f:
                    f.writelines(new_lines)
                print(f"✅ Successfully updated license in {dest_config}")
            else:
                print(f"⚠️ 'license =' line not found in {dest_config}. Appending it.")
                with open(dest_config, 'a') as f:
                    f.write(f"\nlicense = '{license_key}'\n")
                print(f"✅ Successfully appended license to {dest_config}")
                
        except Exception as e:
            print(f"❌ Error writing to config file: {e}")
    else:
        print("❌ Could not locate modeller/config.py")

if __name__ == "__main__":
    configure_modeller()
