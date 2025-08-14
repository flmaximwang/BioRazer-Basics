import pandas as pd


def sort_metrics(df: pd.DataFrame, metric_prefix: str, reverse=True, inplace=False):

    metric_columns = [col for col in df.columns if col.startswith(metric_prefix)]
    if not metric_columns:
        raise ValueError(
            f"No columns found with prefix '{metric_prefix}' in the DataFrame."
        )
    if inplace:
        for i in df.index:
            df.loc[i, metric_columns] = sorted(
                df.loc[i, metric_columns], reverse=reverse
            )
        return None
    else:
        sorted_df = df.copy()
        for i in sorted_df.index:
            sorted_df.loc[i, metric_columns] = sorted(
                sorted_df.loc[i, metric_columns], reverse=reverse
            )
        return sorted_df
