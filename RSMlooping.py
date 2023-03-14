"""
Looping over RMS dataset
"""
import click
import os
from subprocess import run

def get_data_files(data_dir, extension = ".txt"):

    res = {}
    for (dirpath, dirnames, filenames) in os.walk(data_dir):
        for filename in filenames:
            if filename.endswith(extension):
                file_path = os.sep.join([dirpath, filename])

                sampleName = filename.split("_")[1]
                res[sampleName] = file_path
    return res

@click.command()
@click.option(
    "--input-dir",
    "-input",
    prompt=True,
    help=(
        "Directory where all RSM txt files from Rigaku are stored"
    )
)
@click.option(
    "--output-dir",
    "-output",
    prompt=True,
    help=(
        "Output directory for calculated data"
    )
)

def run_looping(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    data = get_data_files(input_dir)
    
    for fileName, path in data.items():
        
        print("working on file : {}".format(fileName), flush=True)

        try:

            run([
                "python",
                "RSMfit.py",
                "-input",
                path,
                "-output",
                os.path.join(output_dir, "{}.json".format(fileName))
            ])

        
        except Exception as expt_msg:
            print("this calculation ({}) failed: {}".format(
                fileName,
                expt_msg), flush=True)
    
    


if __name__ == "__main__":
    run_looping()

