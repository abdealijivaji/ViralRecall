import argparse
from pathlib import Path
import shutil
import time
import requests
import progressbar
import numpy as np
from .utils import prep_hmm
import hashlib
import os


def check_file_hash(filepath: Path, expected_hash: str) -> bool:
    """Check the md5 hash of a file to ensure it matches the expected hash"""

    with open(filepath, "rb") as f:
        file_hash = hashlib.md5(f.read()).hexdigest()
    return file_hash == expected_hash


def download_file(url: str, filepath: Path) -> None:
    compressed_hash = "4a1926b1b4ac9d6d60e521a28c64c15d"
    if filepath.exists() and check_file_hash(filepath, compressed_hash):
        print("Database is already downloaded")
        return
    n_chunk = 1000
    r = requests.get(url, stream=True)

    # Estimates the number of bar updates
    block_size = 1024
    file_size = int(r.headers.get("Content-Length", 0))
    print(f"Downloading database of total size: {int(file_size / (block_size**2))} MB")
    num_bars = np.ceil(file_size / (n_chunk * block_size))
    bar = progressbar.ProgressBar(maxval=num_bars).start()
    with open(filepath, "wb") as f:
        for i, chunk in enumerate(r.iter_content(chunk_size=n_chunk * block_size)):
            f.write(chunk)
            bar.update(i + 1)
            time.sleep(0.05)
    if filepath.exists() and check_file_hash(filepath, compressed_hash):
        print("Finished downloading HMM database")
    else:
        print("Error in downloading database. Please try again.")
    return


def parse_args(argv=None):
    args_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Download database for Viralrecall",
        epilog="*******************************************************************\n\n*******************************************************************",
    )
    args_parser.add_argument(
        "-d",
        "--db_dir",
        required=False,
        default=Path.cwd(),
        help="Directory to download the HMM database. Default is current working directory",
    )
    args_parser.add_argument(
        "-n",
        "--name",
        required=False,
        default="VR_hmm_database",
        help='Name to give the database directory. Default is "VR_hmm_database"',
    )
    args_parser = args_parser.parse_args()

    return args_parser


def main():
    args_list = parse_args()
    out_dir = Path(args_list.db_dir).expanduser()
    project_name = args_list.name

    print(f"Downloading database to {out_dir}")
    dbsrc = "https://zenodo.org/records/20401448/files/VR_hmm_database.tar.gz"

    db_dir = out_dir / project_name
    db_file = db_dir.with_suffix(".tar.gz")
    download_file(dbsrc, db_file)

    if not db_dir.is_dir():
        shutil.unpack_archive(db_file, out_dir, format="gztar")
        print("Finished unpacking")
    os.remove(db_file)
    os.rename(out_dir / "VR_hmm_database", db_dir)
    print(f"Preparing HMM database: {db_dir}")
    prep_hmm(db_dir)


if __name__ == "__main__":
    main()
