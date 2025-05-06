import scipy.stats as stats
import pandas as pd


def fisher_exact_test(
    sample1_meth: int, sample1_unmeth: int, sample2_meth: int, sample2_unmeth: int
):
    """
    Perform Fisher's Exact Test to compare methylation counts between two samples.

    Parameters:
    sample1_meth (int): Methylated CpG count for sample 1.
    sample1_unmeth (int): Unmethylated CpG count for sample 1.
    sample2_meth (int): Methylated CpG count for sample 2.
    sample2_unmeth (int): Unmethylated CpG count for sample 2.

    Returns:
    dict: A dictionary containing:
        - 'odds_ratio': Odds Ratio from Fisher's Exact Test
        - 'fisher_p_value': p-value from Fisher's Exact Test
    """
    table = [
        [sample1_meth, sample1_unmeth],
        [sample2_meth, sample2_unmeth],
    ]

    odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")

    return {"odds_ratio": odds_ratio, "fisher_p_value": p_value}
