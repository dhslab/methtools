import scipy.stats as stats
import pandas as pd


def fisher_exact_test(
    group1_meth: int, group1_unmeth: int, group2_meth: int, group2_unmeth: int
):
    """
    Perform Fisher's Exact Test to compare methylation counts between two groups.

    Parameters:
    group1_meth (int): Methylated CpG count for group 1.
    group1_unmeth (int): Unmethylated CpG count for group 1.
    group2_meth (int): Methylated CpG count for group 2.
    group2_unmeth (int): Unmethylated CpG count for group 2.

    Returns:
    dict: A dictionary containing:
        - 'odds_ratio': Odds Ratio from Fisher's Exact Test
        - 'fisher_p_value': p-value from Fisher's Exact Test
    """
    table = [
        [group1_meth, group1_unmeth],
        [group2_meth, group2_unmeth],
    ]

    odds_ratio, p_value = stats.fisher_exact(table, alternative="two-sided")

    return {"odds_ratio": odds_ratio, "fisher_p_value": p_value}
