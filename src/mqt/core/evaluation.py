"""Evaluating the json file generated by the benchmarking script."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

if TYPE_CHECKING:
    from os import PathLike

# Avoid output truncation
pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_rows", None)
pd.set_option("display.width", None)

sort_options = ["ratio", "algorithm"]
higher_better_metrics = ["hits", "hit_ratio"]


class Bcolors:
    """Class for colored output in the terminal."""

    OKGREEN = "\033[92m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"


def __flatten_dict(d: dict[Any, Any], parent_key: str = "", sep: str = ".") -> dict[str, Any]:
    """Flatten a nested dictionary. Every value only has one key which is the path to the value."""
    items = {}
    for key, value in d.items():
        new_key = f"{parent_key}{sep}{key}" if parent_key else key
        if isinstance(value, dict):
            items.update(__flatten_dict(value, new_key, sep=sep))
        else:
            items[new_key] = value
    return items


def __post_processing(key: str) -> dict[str, str]:
    """Postprocess the key of a flattened dictionary to get the metrics for the DataFrame columns."""
    metrics_divided = key.split(".")
    result_metrics = {}
    if len(metrics_divided) < 4:
        raise ValueError("Benchmark " + key + " is missing algorithm, task, number of qubits or metric!")
    result_metrics["algorithm"] = metrics_divided.pop(0)
    result_metrics["task"] = metrics_divided.pop(0)
    result_metrics["num_qubits"] = metrics_divided.pop(0)
    num_remaining_benchmarks = len(metrics_divided)
    if num_remaining_benchmarks == 1:
        result_metrics["component"] = ""
        result_metrics["metric"] = metrics_divided.pop(0)
    elif num_remaining_benchmarks == 2:
        if metrics_divided[0] == "dd":
            result_metrics["component"] = "" if metrics_divided[0] == "dd" else metrics_divided.pop(0)
            result_metrics["metric"] = metrics_divided[-1]
    else:
        separator = "_"
        # if the second-to-last element is not "total" then only the last element is the metric and the rest component
        if metrics_divided[-2] == "total":
            metric = separator.join(metrics_divided[-2:])
            result_metrics["metric"] = metric
            component = separator.join(metrics_divided[:-2])
            if component.startswith("dd_"):
                component = component[3:]
            result_metrics["component"] = component
        else:
            result_metrics["metric"] = metrics_divided[-1]
            component = separator.join(metrics_divided[:-1])
            if component.startswith("dd_"):
                component = component[3:]
            result_metrics["component"] = component

    return result_metrics


def __aggregate(baseline_filepath: str | PathLike[str], feature_filepath: str | PathLike[str]) -> pd.DataFrame:
    """Aggregate the data from the baseline and feature json files into one DataFrame for visualization."""
    base_path = Path(baseline_filepath)
    with base_path.open(mode="r", encoding="utf-8") as f:
        d = json.load(f)
    flattened_data = __flatten_dict(d)
    feature_path = Path(feature_filepath)
    with feature_path.open(mode="r", encoding="utf-8") as f:
        d_feature = json.load(f)
    flattened_feature = __flatten_dict(d_feature)

    for k, v in flattened_data.items():
        value = v
        if value == "unused":
            value = float("nan")
        if k in flattened_feature:
            ls = [value, flattened_feature[k]]
            flattened_data[k] = ls
            del flattened_feature[k]
        else:
            ls = [value, float("nan")]
            flattened_data[k] = ls
    # If a benchmark is in the feature file but not in the baseline file, it should be added with baseline marked as
    # "skipped"
    for k, v in flattened_feature.items():
        value = v
        if value == "unused":
            value = float("nan")
        ls = [float("nan"), value]
        flattened_data[k] = ls

    before_ls, after_ls, ratio_ls, algorithm_ls, task_ls, num_qubits_ls, component_ls, metric_ls = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )

    for k, v in flattened_data.items():
        before = v[0]
        after = v[1]
        if math.isnan(before) or math.isnan(after):
            ratio = float("nan")
        else:
            ratio = after / before if before != 0 else 1 if after == 0 else math.inf
        key = k
        if k.endswith(tuple(higher_better_metrics)):
            ratio = 1 / ratio
            key += "*"
        before_ls.append(round(before, 3) if isinstance(before, float) else before)
        after_ls.append(round(after, 3) if isinstance(after, float) else after)
        ratio_ls.append(round(ratio, 3))

        # postprocessing
        result_metrics = __post_processing(key)
        algorithm_ls.append(result_metrics["algorithm"])
        task_ls.append(result_metrics["task"])
        num_qubits_ls.append(result_metrics["num_qubits"])
        component_ls.append(result_metrics["component"])
        metric_ls.append(result_metrics["metric"])

    df_all = pd.DataFrame()
    df_all["before"] = before_ls
    df_all["after"] = after_ls
    df_all["ratio"] = ratio_ls

    df_all["algo"] = algorithm_ls
    df_all["task"] = task_ls
    df_all["n"] = num_qubits_ls
    df_all["component"] = component_ls
    df_all["metric"] = metric_ls
    df_all.index = pd.Index([""] * len(df_all.index))

    return df_all


def compare(
    baseline_filepath: str | PathLike[str],
    feature_filepath: str | PathLike[str],
    factor: float = 0.1,
    sort: str = "ratio",
    dd: bool = False,
    only_changed: bool = False,
    no_split: bool = False,
    algorithm: str | None = None,
    task: str | None = None,
    num_qubits: str | None = None,
) -> None:
    """Compare the results of two benchmarking runs from the generated json file.

    Args:
        baseline_filepath: Path to the baseline json file.
        feature_filepath: Path to the feature json file.
        factor: How much a result has to change to be considered significant.
        sort: Sort the table by this column. Valid options are "ratio" and "experiment".
        dd: Whether to show the detailed dd benchmarks.
        only_changed: Whether to only show results that changed significantly.
        no_split: Whether to merge all results together in one table or to separate the results into benchmarks that improved, stayed the same, or worsened.
        algorithm: Only show results for this algorithm.
        task: Only show results for this task.
        num_qubits: Only show results for this number of qubits. Can only be used if algorithm is also specified.

    Returns:
        None
    Raises:
        ValueError: If factor is negative or sort is invalid or if num_qubits is specified while algorithm is not.
        FileNotFoundError: If the baseline_filepath argument or the feature_filepath argument does not point to a valid file.
        JSONDecodeError: If the baseline_filepath argument or the feature_filepath argument points to a file that is not a valid JSON file.
    """
    if factor < 0:
        msg = "Factor must be positive!"
        raise ValueError(msg)
    if sort not in sort_options:
        msg = "Invalid sort option! Valid options are 'ratio' and 'algorithm'."
        raise ValueError(msg)
    if algorithm is None and num_qubits is not None:
        msg = "num_qubits can only be specified if algorithm is also specified!"
        raise ValueError(msg)

    df_all = __aggregate(baseline_filepath, feature_filepath)

    if task is not None:
        df_all = df_all[df_all["task"].str.contains(task, case=False)]
    if algorithm is not None:
        df_all = df_all[df_all["algo"].str.contains(algorithm, case=False)]
    if num_qubits is not None:
        df_all = df_all[df_all["n"] == num_qubits]

    df_runtime = df_all[df_all["metric"] == "runtime"]
    df_runtime = df_runtime.drop(columns=["component", "metric"])
    m1_runtime = df_runtime["ratio"] < 1 - factor
    m2_runtime = df_runtime["ratio"] > 1 + factor
    m3_runtime = (df_runtime["ratio"] != df_runtime["ratio"]) | (
        (1 - factor < df_runtime["ratio"]) & (df_runtime["ratio"] < 1 + factor)
    )

    print("\nDD runtimes:")
    if no_split:
        if only_changed:
            df_runtime = df_runtime[m1_runtime | m2_runtime]
            print("\nAll changed runtimes:\n")
        else:
            print("\nAll runtimes:\n")
        print(df_runtime.to_markdown(index=False, stralign="right"))
    else:
        print(f"\n{Bcolors.OKGREEN}Runtimes that have improved:{Bcolors.ENDC}\n")
        df_runtime_improved = df_runtime[m1_runtime]
        df_runtime_improved = (
            df_runtime_improved.sort_values(by=sort)
            if sort == "ratio"
            else df_runtime_improved.sort_values(["algo", "task", "n"])
        )
        print(df_runtime_improved.to_markdown(index=False, stralign="right"))

        print(f"\n{Bcolors.FAIL}Runtimes that have worsened:{Bcolors.ENDC}\n")
        df_runtime_worsened = df_runtime[m2_runtime]
        df_runtime_worsened = (
            df_runtime_worsened.sort_values(by=sort, ascending=False)
            if sort == "ratio"
            else df_runtime_worsened.sort_values(["algo", "task", "n"])
        )
        print(df_runtime_worsened.to_markdown(index=False, stralign="right"))

        if not only_changed:
            print("\nRuntimes that have stayed the same:\n")
            df_runtime_same = df_runtime[m3_runtime]
            df_runtime_same = (
                df_runtime_same.sort_values(by=sort)
                if sort == "ratio"
                else df_runtime_same.sort_values(["algo", "task", "n"])
            )
            print(df_runtime_same.to_markdown(index=False, stralign="right"))
    if dd:
        print("\nDD details:")
        df_all = df_all[df_all["metric"] != "runtime"]

        m1 = df_all["ratio"] < 1 - factor  # after significantly smaller than before
        m2 = df_all["ratio"] > 1 + factor  # after significantly larger than before
        m3 = (df_all["ratio"] != df_all["ratio"]) | ((1 - factor < df_all["ratio"]) & (df_all["ratio"] < 1 + factor))
        # ratio is NaN or after not significantly different from before

        if no_split:
            if only_changed:
                df_all = df_all[m1 | m2]
                print("\nAll changed DD benchmarks:\n")
            else:
                print("\nAll DD benchmarks:\n")

            df_all = (
                df_all.sort_values(by=sort)
                if sort == "ratio"
                else df_all.sort_values(["algo", "task", "n", "component", "metric"])
            )
            print(df_all.to_markdown(index=False, stralign="right"))
        else:
            print(f"\n{Bcolors.OKGREEN}DD Benchmarks that have improved:{Bcolors.ENDC}\n")
            df_improved = df_all[m1]
            df_improved = (
                df_improved.sort_values(by=sort)
                if sort == "ratio"
                else df_improved.sort_values(["algo", "task", "n", "component", "metric"])
            )
            print(df_improved.to_markdown(index=False, stralign="right"))

            print(f"\n{Bcolors.FAIL}DD Benchmarks that have worsened:{Bcolors.ENDC}\n")
            df_worsened = df_all[m2]
            df_worsened = (
                df_worsened.sort_values(by=sort, ascending=False)
                if sort == "ratio"
                else df_worsened.sort_values(["algo", "task", "n", "component", "metric"])
            )
            print(df_worsened.to_markdown(index=False, stralign="right"))

            if not only_changed:
                print("\nDD Benchmarks that have stayed the same:\n")
                df_same = df_all[m3]
                df_same = (
                    df_same.sort_values(by=sort)
                    if sort == "ratio"
                    else df_same.sort_values(["algo", "task", "n", "component", "metric"])
                )
                print(df_same.to_markdown(index=False, stralign="right"))


def main() -> None:
    """Main function for the command line interface."""
    parser = argparse.ArgumentParser(
        description="Compare the results of two benchmarking runs from the generated json files."
    )
    parser.add_argument("baseline_filepath", type=str, help="Path to the baseline json file.")
    parser.add_argument("feature_filepath", type=str, help="Path to the feature json file.")
    parser.add_argument(
        "--factor", type=float, default=0.1, help="How much a result has to change to be considered significant."
    )
    parser.add_argument(
        "--sort",
        type=str,
        default="ratio",
        help="Sort the table by this column. Valid options are 'ratio' and 'algorithm'.",
    )
    parser.add_argument("--dd", action="store_true", help="Whether to show the detailed dd benchmarks.")
    parser.add_argument(
        "--only_changed", action="store_true", help="Whether to only show results that changed significantly."
    )
    parser.add_argument(
        "--no_split",
        action="store_true",
        help="Whether to merge all results together in one table or to separate the results into "
        "benchmarks that improved, stayed the same, or worsened.",
    )
    parser.add_argument("--algorithm", type=str, help="Only show results for this algorithm.")
    parser.add_argument("--task", type=str, help="Only show results for this task.")
    parser.add_argument(
        "--num_qubits",
        type=str,
        help="Only show results for this number of qubits. Can only be used " "if algorithm is also specified.",
    )
    args = parser.parse_args()
    assert args is not None
    compare(
        args.baseline_filepath,
        args.feature_filepath,
        args.factor,
        args.sort,
        args.dd,
        args.only_changed,
        args.no_split,
        args.algorithm,
        args.task,
        args.num_qubits,
    )
