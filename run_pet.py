#!/usr/bin/env python3

import argparse
import os
import subprocess
from glob import glob
from subprocess import Popen, PIPE

def run(command, env=None, ignore_errors=False):
    if env is None:
        env = {}
    merged_env = os.environ.copy()
    merged_env.update(env)
    merged_env.pop('DEBUG', None)  # Prevent large logs from FreeSurfer

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, env=merged_env)
    output = []
    while True:
        line = process.stdout.readline()
        if not line and process.poll() is not None:
            break
        line = line.decode('utf-8').strip()
        if line:
            print(line)
            output.append(line)
    
    if process.returncode != 0 and not ignore_errors:
        raise Exception(f"Command failed with return code {process.returncode}")
    
    return "\n".join(output)


# Function to import,align and register PET data to T1
def process_pet_file(pet_file, fs_pet_dir, t1_mgz, template_mgz, subject_label, n_cpus, smooth):
    subject_label = f"sub-{subject_label}"  # Ensure the correct prefix is added
    num_frames = int(run(f"mri_info {pet_file} --nframes"))

    for i in range(num_frames):
        frame_file = os.path.join(fs_pet_dir, f"frame_{i}.mgz")
        run(f"mri_convert {pet_file} --frame {i} {frame_file}")

        reg_file = os.path.join(fs_pet_dir, f"frame_{i}.reg.lta")
        run(f"mri_coreg --s {subject_label} --mov {frame_file} --reg {reg_file} --threads {n_cpus}")

        aligned_frame_file = os.path.join(fs_pet_dir, f"frame_reg_{i}.mgz")
        run(f"mri_vol2vol --mov {frame_file} --targ {t1_mgz} --lta {reg_file} --o {aligned_frame_file}")
    
    # Proceed with averaging frames and any further steps...
    template_mean = os.path.join(fs_pet_dir, 'template_mean.mgz')
    template_reg = os.path.join(fs_pet_dir, 'template.reg.lta')
    template_mgz = os.path.join(fs_pet_dir, 'template.mgz')
    
    # Collect all frame files for concatenation
    frame_reg_files = " ".join(glob(os.path.join(fs_pet_dir, 'frame_reg_*.mgz')))
    run(f"mri_concat --i {frame_reg_files} --o {template_mean} --mean")
    # Apply smoothing to averaged PET
    run(f"mri_convert --fwhm {smooth} {template_mean} {template_mgz}")

    # Register the template to generate a template.reg.lta file
    run(f"mri_coreg --s {subject_label} --mov {template_mgz} --reg {template_reg} --threads {n_cpus}")

    # Clean up the directory by removing all files that start with frame_
    frame_files = glob(os.path.join(fs_pet_dir, 'frame_*'))
    for frame_file in frame_files:
        os.remove(frame_file)

def process_gtmseg(subjects_dir, subject_label, env):
    """
    Check if gtmseg.mgz already exists for a subject. If it doesn't, run the gtmseg command.

    Parameters:
    - subjects_dir: The base directory where subjects' data is stored.
    - subject_label: The label of the subject (without the "sub-" prefix).
    - env: The environment variables to use when running the command.
    """
    # Ensure the subject label has the "sub-" prefix
    subject_label_prefixed = f"sub-{subject_label}"

    # Define the path to the gtmseg output file
    gtmseg_output = os.path.join(subjects_dir, subject_label_prefixed, "mri", "gtmseg.mgz")

    # Check if gtmseg.mgz already exists
    if os.path.exists(gtmseg_output):
        print(f"gtmseg.mgz already exists for subject {subject_label_prefixed}. Skipping gtmseg processing.")
    else:
        # Run the gtmseg command
        run(f"gtmseg --s {subject_label_prefixed} --xcerseg", env=env)

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
    left_region = f"Left_{suvr_ref.replace('-', '_')}"
    right_region = f"Right_{suvr_ref.replace('-', '_')}"
    if left_region in lut_dict and right_region in lut_dict:
        return lut_dict[left_region], lut_dict[right_region]
    else:
        raise ValueError(f"Invalid SUVR reference region: {suvr_ref}")
    
def process_pvc(template_mgz, fs_pet_dir, fwhm, pvc_dir, surf_dir, tracer, env, rescale_arg, gtmseg_file):
    run(f"mri_gtmpvc --i {template_mgz} "
        f"--reg {os.path.join(fs_pet_dir, 'template.reg.lta')} "
        f"--psf {fwhm} "
        f"--seg {gtmseg_file} "
        f"--default-seg-merge "
        f"--auto-mask 1 .05 "
        f"--mgx .25 "
        f"{rescale_arg} "  # Add the rescale argument here
        f"--save-input "
        f"--o {pvc_dir}", env=env)

    # Map the PVC result to the cortical surface
    run(f"mri_vol2surf --src {os.path.join(pvc_dir, 'mgx.ctxgm.nii.gz')} "
        f"--srcreg {os.path.join(pvc_dir, 'aux/bbpet2anat.lta')} "
        f"--hemi lh --projfrac 1 "
        f"--o {os.path.join(surf_dir, f'lh_pvc_{tracer}.nii.gz')} --no-reshape", env=env)

    run(f"mri_vol2surf --src {os.path.join(pvc_dir, 'mgx.ctxgm.nii.gz')} "
        f"--srcreg {os.path.join(pvc_dir, 'aux/bbpet2anat.lta')} "
        f"--hemi rh --projfrac 1 "
        f"--o {os.path.join(surf_dir, f'rh_pvc_{tracer}.nii.gz')} --no-reshape", env=env)

def process_nopvc(template_mgz, fs_pet_dir, fwhm, nopvc_dir, surf_dir, tracer, env, rescale_arg, gtmseg_file):
    run(f"mri_gtmpvc --i {template_mgz} "
        f"--reg {os.path.join(fs_pet_dir, 'template.reg.lta')} "
        f"--psf 0 "
        f"--seg {gtmseg_file} "
        f"--default-seg-merge "
        f"--auto-mask 1 .05 "
        f"--mgx .25 "
        f"{rescale_arg} "  # Add the rescale argument here
        f"--save-input "
        f"--o {nopvc_dir}", env=env)

    # Map the PVC result to the cortical surface
    run(f"mri_vol2surf --src {os.path.join(nopvc_dir, 'mgx.ctxgm.nii.gz')} "
        f"--srcreg {os.path.join(nopvc_dir, 'aux/bbpet2anat.lta')} "
        f"--hemi lh --projfrac 1 "
        f"--o {os.path.join(surf_dir, f'lh_nopvc_{tracer}.nii.gz')} --no-reshape", env=env)

    run(f"mri_vol2surf --src {os.path.join(nopvc_dir, 'mgx.ctxgm.nii.gz')} "
        f"--srcreg {os.path.join(nopvc_dir, 'aux/bbpet2anat.lta')} "
        f"--hemi rh --projfrac 1 "
        f"--o {os.path.join(surf_dir, f'rh_nopvc_{tracer}.nii.gz')} --no-reshape", env=env)


# Parse command-line arguments
parser = argparse.ArgumentParser(description='Set up environment and locate the desired PET file.')
parser.add_argument('--subjects_dir', required=True, help='Directory containing subject data, including recon-all outputs.')
parser.add_argument('--subject', required=True, help='Subject ID to process.')
parser.add_argument('--tracer', required=True, help='Tracer used in the PET dataset.')
parser.add_argument('--fs_license', required=True, help='Path to the FreeSurfer license key file.')
parser.add_argument('--pet_dir', required=True, help='Directory containing raw PET data formatted according to BIDS.')
parser.add_argument('--session_label', nargs='+', help='Session labels to analyze. If not provided, all sessions will be analyzed.')
parser.add_argument('--n_cpus', help='Number of CPUs to use for parallel processing.', type=int, default=1)
parser.add_argument('--smooth', help='Smoothing kernel size (FWHM) for mri_convert.', type=float, default=2.8)
parser.add_argument('--suvr_ref', type=str, default="Cerebellum-Cortex", help='Reference region for SUVR calculation. Specify the region name as it appears in FreeSurferColorLUT.txt (e.g., "Cerebellum-Cortex").')
parser.add_argument('--fwhm', type=float, default=2.8, help='FWHM value for PVC (default: 2.8)')
parser.add_argument('--skip_pet_processing', action='store_true', help='Skip the PET processing steps and start from GTMseg')


args = parser.parse_args()

# Determine the subject directory based on whether a PET session is specified
if args.session_label:  # This implies a PET session is specified
    pet_session_label = args.session_label[0]  # Assume the first label is for PET
    if len(args.session_label) > 1:
        mri_session_label = args.session_label[1]  # Assume the second label is for MRI
    else:
        raise ValueError("If a PET session is specified, an MRI session must also be specified.")
    
    # Construct subject directory with MRI session label
    subject_dir = os.path.join(args.subjects_dir, f"sub-{args.subject}_ses-{mri_session_label}")
else:
    # No PET session specified, assume subject directory is without session label
    subject_dir = os.path.join(args.subjects_dir, f"sub-{args.subject}")

# Define fs_pet_dir within the subject directory
fs_pet_dir = os.path.join(subject_dir, "pet", args.tracer)

# Ensure the fs_pet_dir exists
os.makedirs(fs_pet_dir, exist_ok=True)

# Define the T1 and surf directories within the subject directory
t1_mgz = os.path.join(args.subjects_dir, f"sub-{args.subject}", 'mri', 'T1.mgz')
surf_dir = os.path.join(args.subjects_dir, f"sub-{args.subject}", "surf")

# Define template_mgz, pvc_dir, and nopvc_dir relative to the tracer subfolder in fs_pet_dir
template_mgz = os.path.join(fs_pet_dir, "template.mgz")
pvc_dir = os.path.join(fs_pet_dir, "pvc")
nopvc_dir = os.path.join(fs_pet_dir, "nopvc")
gtmseg_file = os.path.join(subject_dir, "mri", "gtmseg.mgz")

# Create the PVC and No-PVC directories if they don't exist
os.makedirs(pvc_dir, exist_ok=True)
os.makedirs(nopvc_dir, exist_ok=True)

# Set up environment variables
os.environ['SUBJECTS_DIR'] = args.subjects_dir
os.environ['FREESURFER_HOME'] = '/export01/data/toolboxes/freesurfer'  # Replace with the actual path

# Validate FreeSurfer license
if not os.path.exists(args.fs_license):
    raise FileNotFoundError("Provided FreeSurfer license file does not exist")
os.environ['FS_LICENSE'] = os.path.abspath(args.fs_license)

# Define the environment variable dictionary
env = {
    'SUBJECTS_DIR': args.subjects_dir,
    'FREESURFER_HOME': os.environ['FREESURFER_HOME'],
    'FS_LICENSE': os.environ['FS_LICENSE']
}

# Load the FreeSurferColorLUT.txt file
freesurfer_home = os.environ.get('FREESURFER_HOME', '/export01/data/toolboxes/freesurfer')
lut_file = os.path.join(freesurfer_home, 'FreeSurferColorLUT.txt')

if not os.path.exists(lut_file):
    raise FileNotFoundError(f"FreeSurferColorLUT.txt not found at {lut_file}")

lut_dict = parse_lut_file(lut_file)

# Get the region numbers for the SUVR reference region
suvr_ref = args.suvr_ref if args.suvr_ref else "Cerebellum_Cortex"
left_ref, right_ref = get_suvr_ref_numbers(lut_dict, suvr_ref)

if left_ref and right_ref:
    rescale_arg = f"--rescale {left_ref} {right_ref}"
else:
    raise ValueError(f"Could not find SUVR reference regions for {suvr_ref} in FreeSurferColorLUT.txt")


# Determine if the study is longitudinal by checking for multiple valid sessions
subject_session_dirs = glob(os.path.join(args.pet_dir, f"sub-{args.subject}", "ses-*"))
sessions = [os.path.basename(session_dir).split("-")[-1] for session_dir in subject_session_dirs]

# If the user specifies session labels, use them; otherwise, process all sessions
if args.session_label:
    sessions_to_analyze = set(sessions).intersection(set(args.session_label))
else:
    sessions_to_analyze = sessions

if not sessions_to_analyze:
    raise Exception(f"No valid sessions found for subject {args.subject}.")

# Find and process PET files for each session
for session_label in sessions_to_analyze:
    # Construct the pattern for the PET file
    pet_file_pattern = os.path.join(args.pet_dir, f"sub-{args.subject}", f"ses-{session_label}", "pet", f"sub-{args.subject}_ses-{session_label}_trc-{args.tracer}_pet.nii.gz")

    # Find the PET file
    pet_files = glob(pet_file_pattern)

    if not pet_files:
        print(f"No PET file found for subject {args.subject}, session {session_label}. Skipping.")
        continue

    # Print found PET files for verification
    for pet_file in pet_files:
        print(f"Found PET file for subject {args.subject}, session {session_label}: {pet_file}")

  # Define the target PET file path in FS_PET_DIR
    target_pet_file = os.path.join(fs_pet_dir, f"{args.tracer}.mgz")

    # Copy and rename the PET file to FS_PET_DIR
    print(f"Copying and renaming PET file to: {target_pet_file}")
    run(f"mri_convert {pet_file} {target_pet_file}")

# Print the defined directories and files for verification
print(f"Defined directories and files for subject {args.subject}:")
print(f"FS_PET_DIR: {fs_pet_dir}")
print(f"T1: {t1_mgz}")
print(f"TEMPLATE: {template_mgz}")
print(f"SURF: {surf_dir}")
print(f"PVC: {pvc_dir}")
print(f"NOPVC: {nopvc_dir}")

print(f"Environment setup and PET file location completed for subject {args.subject}.")

if not args.skip_pet_processing:
    # Process the PET file
    process_pet_file(target_pet_file, fs_pet_dir, t1_mgz, template_mgz, args.subject, args.n_cpus, args.smooth)

# Create anatomical segementation for the geometric transfer matrix (GTM) using gtmseg
process_gtmseg(args.subjects_dir, args.subject, env)

# Process the PVC and No-PVC data
process_pvc(template_mgz, fs_pet_dir, args.psf, pvc_dir, surf_dir, args.tracer, env, rescale_arg, gtmseg_file)
process_nopvc(template_mgz, fs_pet_dir, args.psf, nopvc_dir, surf_dir, args.tracer, env, rescale_arg, gtmseg_file)