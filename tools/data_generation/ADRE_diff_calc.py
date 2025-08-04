import os
import fnmatch
from cdo import Cdo as CDO
os.environ['CDO'] = '/appl/spack/v018/install-tree/gcc-8.5.0/cdo-2.0.5-zpo6xz/bin/cdo'
cdo = CDO()

"""
by Atte Laakso / Aalto University

Calculates difference between two ADRE files using CDO and renames output variable to dADRE.
"""

# Path to directory with simulated results
path_to_goal_set = '../output_DRE'

def main():
    target_paths = []
    ref_files = []
    # Search all ADRE files
    for root, dirs, files in os.walk(path_to_goal_set):
        if any(skip in root for skip in ['old', 'temporary']):
            continue

        found_files = [f for f in files if fnmatch.fnmatch(f, "*ADRE*.nc")]
        full_paths = [os.path.join(root, f) for f in found_files]

        if len(found_files) == 12:
            target_paths.append((root, full_paths))
            if not any(tag in str(found_files).lower() for tag in ['rh', 'mmr', 'ref', 'clt']):
                ref_files.append((root, full_paths))

    print(f"Found {len(target_paths)} directories with valid target files")
    print(f"Found {len(ref_files)} directories with valid reference files")

    # Map reference files (files based on only vbs sensitivity data) by month
    ref_dir = {}
    for root, files in ref_files:
        for f in files:
            try:
                month = f.split("ADRE_")[1].split("2010")[1][:2]
                ref_dir[month] = (root, f)
            except IndexError:
                continue

    total_files = 0
    skipped = 0
    failed = 0
    # Calculate dADRE
    for dir_path, files in target_paths:
        for filename in files:
            if 'ADRE' not in filename:
                skipped += 1
                continue

            total_files += 1
            try:
                month = filename.split("ADRE_")[1].split("2010")[1][:2]
                if month not in ref_dir:
                    print(f"No reference file for month {month}")
                    failed += 1
                    continue

                ref_root, ref_name = ref_dir[month]
                ref_file = os.path.join(ref_root, ref_name)
                file_path = os.path.join(dir_path, filename)

                temp_file = os.path.join(dir_path, f"{month}_temp.nc")
                output_file = os.path.join(dir_path, filename.replace("ADRE", "DELTA"))

                print(f"Calculating difference for: {filename} using {ref_name}")

                # Subtract ADRE values
                cdo.sub(input=f"{ref_file} {file_path}", output=temp_file)

                # Rename variable
                cdo.chname("ADRE,dADRE", input=temp_file, output=output_file)

                # Cleanup
                if os.path.exists(temp_file):
                    os.remove(temp_file)

            except Exception as e:
                print(f"Failed to process {filename}: {e}")
                failed += 1

    print("\n Summary:")
    print(f"  Total processed: {total_files}")
    print(f"  Skipped:         {skipped}")
    print(f"  Failed:          {failed}")
    print("Done.")

if __name__ == "__main__":
    main()
            
