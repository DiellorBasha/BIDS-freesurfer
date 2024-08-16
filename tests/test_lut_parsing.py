# test_lut_parsing.py

def parse_lut_file(lut_path):
    lut_dict = {}
    with open(lut_path, 'r') as file:
        for line in file:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) > 1:
                # Extract the number and region name
                region_number = parts[0]
                region_name = parts[1].replace('-', '_')
                lut_dict[region_name] = region_number
    return lut_dict

def get_suvr_ref_numbers(lut_dict, suvr_ref):
    left_ref = f"Left_{suvr_ref.replace('-', '_')}"
    right_ref = f"Right_{suvr_ref.replace('-', '_')}"
    return lut_dict.get(left_ref), lut_dict.get(right_ref)

# Test the functions
if __name__ == "__main__":
    lut_path = "/export01/data/toolboxes/freesurfer/FreeSurferColorLUT.txt"  # Replace with the actual path
    lut_dict = parse_lut_file(lut_path)
    print("LUT Dictionary:", lut_dict)  # This will print the dictionary

    suvr_ref = "Cerebellum-Cortex"
    left, right = get_suvr_ref_numbers(lut_dict, suvr_ref)
    print(f"SUVR Reference Numbers for {suvr_ref}: Left = {left}, Right = {right}")
