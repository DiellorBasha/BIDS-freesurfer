import os
import argparse
import subprocess
import logging

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    """Runs a shell command and logs its output."""
    logging.info(f"Running command: {command}")
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logging.error(f"Command failed with error: {result.stderr.decode('utf-8')}")
        raise Exception(f"Command failed: {command}")
    logging.info(result.stdout.decode('utf-8'))

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run PETSurfer pipeline for PET analysis.")
    parser.add_argument('--bids_dir', required=True, help='BIDS formatted directory with input data.')
    parser.add_argument('--output_dir', required=True, help='Directory to store output data.')
    parser.add_argument('--subject', required=True, help='Subject identifier (e.g., sub-01).')
    parser.add_argument('--pet_modality', default='pet', help='PET modality name in BIDS dataset (default: pet).')
    parser.add_argument('--tracer', required=True, help='Name of the tracer used (e.g., FDG, PIB).')
    parser.add_argument('--fs_license', required=True, help='Path to the FreeSurfer license file.')
    parser.add_argument('--run_recon', action='store_true', help='Flag to run recon-all if anatomical processing is needed.')
    return parser.parse_args()

def main():
    setup_logging()
    args = parse_arguments()

    # Set FreeSurfer environment variables
    os.environ['SUBJECTS_DIR'] = os.path.join(args.output_dir, 'freesurfer')
    os.environ['FS_LICENSE'] = args.fs_license

    subject_id = args.subject
    pet_dir = os.path.join(args.bids_dir, subject_id, args.pet_modality)
    mri_dir = os.path.join(args.bids_dir, subject_id, 'anat')

    # Ensure directories exist
    if not os.path.exists(pet_dir):
        raise Exception(f"PET directory does not exist: {pet_dir}")
    if not os.path.exists(mri_dir):
        raise Exception(f"MRI directory does not exist: {mri_dir}")

    # Step 1: Convert PET to mgz
    pet_file = os.path.join(pet_dir, f"{subject_id}_task-rest_acq-{args.tracer}_pet.nii.gz")
    pet_mgz = os.path.join(pet_dir, f"{subject_id}_task-rest_acq-{args.tracer}_pet.mgz")
    run_command(f"mri_convert {pet_file} {pet_mgz}")

    # Step 2: Run recon-all if needed
    if args.run_recon:
        recon_all_cmd = f"recon-all -s {subject_id} -all"
        run_command(recon_all_cmd)

    # Step 3: Register PET to anatomical T1
    registration_file = os.path.join(pet_dir, f"{subject_id}_task-rest_acq-{args.tracer}_register.dat")
    run_command(f"bbregister --s {subject_id} --mov {pet_mgz} --reg {registration_file} --init-fsl --t2")

    # Step 4: Quantify the PET data using PETSurfer
    quant_output = os.path.join(args.output_dir, 'petsurfer', subject_id)
    run_command(f"mri_gtmpvc --i {pet_mgz} --reg {registration_file} --psf 5.0 --seg {quant_output}")

    logging.info("PETSurfer pipeline completed successfully.")

if __name__ == "__main__":
    main()
