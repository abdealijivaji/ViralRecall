import os
from pathlib import Path
from pyfaidx import Fasta, FastaIndexingError
from multiprocessing import cpu_count
from pyhmmer import hmmpress, plan7
from . import __version__


def load_genome(input: Path) -> Fasta:
    try:
        genome_file = Fasta(input)
        is_DNA(genome_file)
    except FastaIndexingError:
        raise FastaIndexingError(
            f"Input file {input.name} is not in Fasta format. Please check input file"
        )
    except ValueError:
        raise ValueError(
            f"{input.name} does not look like a valid DNA sequence. Please check input file"
        )

    return genome_file


def compress_hmm(hmmpath: Path) -> None:
    outname = hmmpath.with_suffix("")

    with plan7.HMMFile(hmmpath) as hf:
        print(f"Compressing {hmmpath.name} into binary format")
        hmmpress(hf, outname)
    print(f"Compressed {hmmpath.name}")


def prep_hmm(hmm_dir: Path) -> None:
    for file in hmm_dir.glob("*.hmm"):
        if not file.with_suffix(".h3m").exists():
            compress_hmm(file)


def check_directory_permissions(directory_path: Path) -> bool:
    if os.access(directory_path, os.R_OK | os.W_OK):
        return True
    else:
        return False


valid_bases = set("ATCGN")


def filt_fasta(phagesize: int, genome_file: Fasta) -> list[str]:
    filt_contig_list: list[str] = []

    for contig in genome_file.keys():
        if len(genome_file[contig]) >= phagesize:
            filt_contig_list.append(contig)

    return filt_contig_list


def is_DNA(genome: Fasta) -> None:
    seq_for_check = set(str(genome[0][:1000]).upper())
    if not set(seq_for_check).issubset(valid_bases):
        raise ValueError


def mp_cpu(cpu: int | None) -> int:
    """Returns number of parallel processes to spawn in batch mode"""
    if cpu is None:
        return cpu_count()
    else:
        return int(cpu) if int(cpu) <= cpu_count() else cpu_count()


def find_db(database: Path) -> Path:

    file_list = ["gvog_mirus_cat.hmm", "NCLDV_markers.hmm", "gvog_annotation.tsv"]

    if database.is_file() and all((database.parent / i).is_file() for i in file_list):
        print(
            f"Note: Using {database.parent} as database.\nFor future runs, please provide the database directory and not individual files."
        )
        return database.parent
    elif not database.is_dir():
        raise NotADirectoryError(
            f"{database} is not a directory. Please check the database path."
        )
    elif not all((database / i).is_file() for i in file_list):
        raise FileNotFoundError(
            f"{database} has missing files. Please check the database path or redownload using viralrecall_database script."
        )

    return database


def check_version(database: Path) -> bool:
    """Added a version check for the new database with tigtog model to ensure that the database version is compatible with the ViralRecall version.
    The database version (v2.0) is stored in a file called 'vr_db_version' in the database directory.
    If they are not compatible, it raises an error and prompts the user to update the database or use a compatible version of ViralRecall."""

    tigtog_files = [
        "Tigtog_db/names_imp_gvog_fam.csv",
        "Tigtog_db/names_imp_GVOGs_both_levels.csv",
        "Tigtog_db/names_imp_gvog_order.csv",
        "Tigtog_db/retrained_clf_Order.joblib",
    ]
    check_files = all((database / i).is_file() for i in tigtog_files)

    try:
        with open(database / "vr_db_version", "r") as f:
            db_version = f.read().strip()
        if db_version != "v2.0":
            raise ValueError
        elif not check_files:
            raise FileNotFoundError
        else:
            updated_ver = True
    except ValueError:
        updated_ver = False
        print(
            f"Database version {db_version} is not compatible with ViralRecall version {__version__}.\nPlease update the database using viralrecall_database script or use a compatible version of ViralRecall."
        )
    except FileNotFoundError:
        updated_ver = False
        print(
            f"{database} has missing files required for running TIGTOG model. Please check the database path or redownload using viralrecall_database script."
        )

    return updated_ver
