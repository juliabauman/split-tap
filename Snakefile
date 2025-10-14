configfile: "config/config.yaml"

from copy import deepcopy
from itertools import product
from pathlib import Path
import string


def _iter_variable_values(spec):
    """Expand a variable specification into a list of concrete values."""

    if isinstance(spec, (list, tuple)):
        return list(spec)

    if isinstance(spec, dict):
        if "values" in spec:
            return list(spec["values"])

        start = spec.get("start")
        stop = spec.get("stop")

        if start is None or stop is None:
            raise ValueError("range-based variables require 'start' and 'stop' keys")

        step = spec.get("step", 1)
        inclusive = spec.get("inclusive", True)
        fmt = spec.get("format")
        prefix = spec.get("prefix", "")
        suffix = spec.get("suffix", "")

        stop_value = stop + (1 if inclusive else 0)
        values = range(start, stop_value, step)

        if fmt is not None:
            values = [fmt.format(value) for value in values]
        else:
            values = [str(value) for value in values]

        if prefix or suffix:
            values = [f"{prefix}{value}{suffix}" for value in values]

        return values

    return [spec]


_FORMATTER = string.Formatter()


def _format_string(value, mapping):
    """Safely format a single string with the provided mapping."""

    has_fields = False
    for _, field_name, _, _ in _FORMATTER.parse(value):
        if field_name is None:
            continue
        has_fields = True
        if field_name and field_name not in mapping:
            return value

    if not has_fields:
        return value

    return value.format(**mapping)


def _format_structure(structure, mapping):
    """Recursively format strings within a nested structure."""

    if isinstance(structure, str):
        return _format_string(structure, mapping)

    if isinstance(structure, list):
        return [_format_structure(item, mapping) for item in structure]

    if isinstance(structure, dict):
        return {key: _format_structure(value, mapping) for key, value in structure.items()}

    return structure


def _deep_merge(base, update):
    """Deeply merge ``update`` into ``base`` and return the merged dict."""

    result = deepcopy(base)
    for key, value in update.items():
        if (
            key in result
            and isinstance(result[key], dict)
            and isinstance(value, dict)
        ):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result


def _expand_sample_sets(config_dict):
    """Combine explicit samples with any templated sample sets."""

    samples = deepcopy(config_dict.get("samples", {}))
    sample_sets = config_dict.get("sample_sets", []) or []

    for entry in sample_sets:
        template = entry.get("template")
        if not template:
            continue

        sample_key = entry.get("sample_key", "sample")
        variables = entry.get("variables", {})

        if sample_key not in template:
            raise ValueError(
                f"sample set template is missing the '{sample_key}' field: {template}"
            )

        variable_names = list(variables.keys())
        variable_values = [
            _iter_variable_values(variables[name]) for name in variable_names
        ]

        if variable_names:
            combinations = [
                dict(zip(variable_names, values)) for values in product(*variable_values)
            ]
        else:
            combinations = [{}]

        for combination in combinations:
            sample_name_template = template[sample_key]
            if not isinstance(sample_name_template, str):
                raise ValueError(
                    f"sample name template for '{sample_key}' must be a string"
                )

            sample_name = sample_name_template.format(**combination)
            mapping = {**combination, sample_key: sample_name}
            if sample_key != "sample":
                mapping.setdefault("sample", sample_name)

            remaining_template = {
                key: value for key, value in template.items() if key != sample_key
            }

            formatted_template = _format_structure(remaining_template, mapping)

            if sample_name in samples:
                samples[sample_name] = _deep_merge(samples[sample_name], formatted_template)
            else:
                samples[sample_name] = deepcopy(formatted_template)

    return samples


SAMPLES = _expand_sample_sets(config)
OUTPUT_DIR = Path(config.get("paths", {}).get("output_dir", "Data"))

CDNA_SAMPLES = [sample for sample, cfg in SAMPLES.items() if cfg.get("cdna") and cfg.get("cdna", {}).get("enabled", True)]
GUIDE_SAMPLES = [sample for sample, cfg in SAMPLES.items() if cfg.get("guides") and cfg.get("guides", {}).get("enabled", True)]
DUAL_GUIDE_SAMPLES = [sample for sample in GUIDE_SAMPLES if SAMPLES[sample]["guides"].get("mode", "dual").lower() == "dual"]
SINGLE_GUIDE_SAMPLES = [sample for sample in GUIDE_SAMPLES if SAMPLES[sample]["guides"].get("mode", "dual").lower() != "dual"]


def cdna_counts_path(sample):
    return str(OUTPUT_DIR / "cdna" / sample / "counts.tsv.gz")


def guide_counts_path(sample):
    return str(OUTPUT_DIR / "guides" / sample / "guide_counts.tsv")


def combined_cdna_counts_path():
    return str(OUTPUT_DIR / "cdna" / "combined_counts.tsv.gz")


def combined_guide_counts_path():
    return str(OUTPUT_DIR / "guides" / "combined_counts.tsv")


def histogram_path(sample):
    hist_dir = Path(config.get("guide_processing", {}).get("histogram_dir", OUTPUT_DIR / "guides" / "histograms"))
    return str(Path(hist_dir) / f"{sample}_tpt_histogram.tsv")


workflow.globals.update(
    {
        "OUTPUT_DIR": OUTPUT_DIR,
        "SAMPLES": SAMPLES,
        "CDNA_SAMPLES": CDNA_SAMPLES,
        "GUIDE_SAMPLES": GUIDE_SAMPLES,
        "DUAL_GUIDE_SAMPLES": DUAL_GUIDE_SAMPLES,
        "SINGLE_GUIDE_SAMPLES": SINGLE_GUIDE_SAMPLES,
        "cdna_counts_path": cdna_counts_path,
        "guide_counts_path": guide_counts_path,
        "combined_cdna_counts_path": combined_cdna_counts_path,
        "combined_guide_counts_path": combined_guide_counts_path,
        "histogram_path": histogram_path,
    }
)

include: "rules/cdna.smk"
include: "rules/guides.smk"


rule all:
    input:
        "../Data/cdna_counts.tsv.gz"

include: "rules/preprocess.smk"
include: "rules/count.smk"
        ([combined_cdna_counts_path()] if CDNA_SAMPLES else [])
        + [cdna_counts_path(sample) for sample in CDNA_SAMPLES]
        + ([combined_guide_counts_path()] if GUIDE_SAMPLES else [])
        + [guide_counts_path(sample) for sample in GUIDE_SAMPLES]
