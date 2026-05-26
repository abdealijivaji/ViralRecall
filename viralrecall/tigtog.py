"""
This is a modified version of TIGTOG originally written by Dr. Anh D. Ha and is available on <https://github.com/anhd-ha/TIGTOG>
"""

import numpy as np
import pandas as pd
import joblib
import re

# from pyfaidx import Fasta
from pathlib import Path


def get_gc(vreg: str) -> float:
    """Calculate the GC content of a viral region from DNA sequence as string"""
    seq_len = len(vreg)
    cnt = vreg.count("G") + vreg.count("C")
    GC_perc = round(cnt * 100 / seq_len, 2)
    return GC_perc


# def density(table: pd.DataFrame, vreg: Fasta) -> float :
#     prot_lens = pd.Series
#     prot_lens = table["pend"] - table["pstart"]
#     sum_len = sum(prot_lens)
#     dens = round(100 * sum_len / len(vreg[0]), 2)
#     return dens

# def large_seq(vreg: Fasta) -> bool :
#     if len(vreg[0]) > 5_000_000 :
#         return True
#     else:
#         return False


def imp_names(list_file: Path) -> list:
    """Read the list of important GVOGs from a file and return a list of their names with the same format as in the HMM table (with .trim extension)"""
    list_g = list()
    g_names = open(list_file, "r")
    for line in g_names.readlines():
        g = line.rstrip()
        g2 = re.sub("hmm$", "trim", g)
        # g2 = re.sub("_", "" , g2)
        list_g.append(g2)
    return list_g


def parse_hits(table, imp_names) -> dict:
    """Parse the HMM hits from the table and create a dictionary with the count of hits for each important GVOG.
    The keys of the dictionary are the names of the important GVOGs with the same format as in the HMM table (with .trim extension)"""
    gvg_hits = {i: 0 for i in imp_names}
    for hmm, count in table["HMM_hit"].value_counts().items():
        hmm = re.sub("GVOGm", "GVOGm_", hmm)
        hmm = hmm + ".trim"
        if hmm in imp_names:
            gvg_hits[hmm] = count
    return gvg_hits


# hitdict = parse_hits(vreg_tab, imp_name_list)
# print(imp_name_list)


def create_df(gc_perc, hitdict):

    df = pd.DataFrame(
        hitdict,
        index=[
            "seq",
        ],
    )
    df.insert(0, "GC_content", gc_perc)
    # print(df)
    return df


# create_df(get_gc(vreg), hitdict)


def tax_predict(clf_file, input_df):
    """Load the trained classifier from a file and use it to predict the taxonomic classification of the viral region based on the input dataframe.
    The input dataframe should have the same format as the one used to train the classifier, with the same columns and order of features.
    'Sequence name as Index' and 'GC_content' as the first column, followed by the important GVOGs in the same order as in the training data.
    The function returns the predicted taxonomic classification and the confidence of the prediction as a percentage.
    """
    clf = joblib.load(clf_file)
    pred = clf.predict(input_df)[0]
    conf_pred = clf.predict_proba(input_df)
    # class_labels = clf.classes_
    for i, prediction in enumerate(conf_pred):
        max_index = np.argmax(prediction)
        confidence = prediction[max_index]

    return pred, confidence


def run_tigtog(vreg: str, table: pd.DataFrame, data_dir: Path) -> tuple[str, str]:

    GC_perc = get_gc(vreg)

    both_levels = imp_names(data_dir / "Tigtog_db/names_imp_GVOGs_both_levels.csv")
    order_levels = imp_names(data_dir / "Tigtog_db/names_imp_gvog_order.csv")
    # fam_levels = imp_names(data_dir / "Tigtog_db/names_imp_gvog_fam.csv")

    hitdict = parse_hits(table, both_levels)
    input_df = create_df(GC_perc, hitdict)
    Ord_df = input_df[order_levels]
    # Fam_df = input_df[fam_levels]

    ord_file = data_dir / "Tigtog_db/retrained_clf_Order.joblib"

    Ord_pred, Ord_prob = tax_predict(ord_file, Ord_df)

    return Ord_pred, f"{Ord_prob * 100:.2f}"
