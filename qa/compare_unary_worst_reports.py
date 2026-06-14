#!/usr/bin/env python3
import argparse
import csv
import math
import os
import re
from typing import Dict, List, Tuple

REQUIRED_COLUMNS = {
    "build_variant",
    "type",
    "function",
    "operation",
    "domain",
    "input_mpfr",
    "input_limbs",
    "relerr",
    "allowed",
    "pass",
    "condition_number",
    "seed",
    "ieee_add",
    "sloppy_mul",
    "sloppy_div",
    "fma",
}

UNIFIED_PRECISION_ORDER = ["dd", "td", "qd", "edd"]
_UNARY_REPORT_RE = re.compile(r"^oracle_unary_(.+)_worst\.csv$")


def to_float(text: str) -> float:
  try:
    return float(text)
  except (TypeError, ValueError):
    return float('nan')


def parse_report_path(spec: str) -> Tuple[str, str]:
  if "=" not in spec or spec.startswith("="):
    raise argparse.ArgumentTypeError(
        "each --report argument must be 'NAME=PATH' e.g. dd=/tmp/oracle_unary_dd_worst.csv")
  name, path = spec.split("=", 1)
  name = name.strip().lower()
  path = path.strip()
  if not name:
    raise argparse.ArgumentTypeError("empty report name in --report")
  if not path:
    raise argparse.ArgumentTypeError("empty report path in --report")
  return name, path


def read_rows(path: str) -> List[dict]:
  with open(path, newline='') as f:
    reader = csv.DictReader(f)
    fieldnames = set(reader.fieldnames or ())
    missing = sorted(REQUIRED_COLUMNS - fieldnames)
    if missing:
      raise ValueError(f"missing columns in {path}: {', '.join(missing)}")
    rows = list(reader)

  if not rows:
    raise SystemExit(f'report is empty: {path}')

  return rows


def build_index(rows_by_precision: Dict[str, List[dict]]) -> Tuple[List[Tuple[str, str, str]], Dict[Tuple[str, str, str], Dict[str, dict]]]:
  index: Dict[Tuple[str, str, str], Dict[str, dict]] = {}

  for precision, rows in rows_by_precision.items():
    for row in rows:
      key = (row.get('function', '').strip(), row.get('operation', '').strip(), row.get('domain', '').strip())
      if not any(key):
        continue
      relerr = abs(to_float(row.get('relerr', 'nan')))
      copy = dict(row)
      copy['__abs_relerr__'] = relerr
      index.setdefault(key, {})[precision] = copy

  sortable = []
  for key, precision_map in index.items():
    finite = [to_float(rec.get('__abs_relerr__')) for rec in precision_map.values() if math.isfinite(to_float(rec.get('__abs_relerr__')))]
    worst = max(finite) if finite else float('nan')
    sortable.append((key, worst))

  sortable.sort(key=lambda item: (math.isfinite(item[1]) is False, -item[1] if math.isfinite(item[1]) else 0.0, item[0]))

  return [key for key, _ in sortable], index


def format_float(value) -> str:
  try:
    f = float(value)
  except (TypeError, ValueError):
    return str(value) if value else '-'

  if math.isnan(f):
    return 'nan'
  if math.isinf(f):
    return 'inf' if f > 0 else '-inf'
  return f'{f:.16g}'


def stats_for_precision_rows(rows: Dict[str, dict], precisions: List[str]) -> Tuple[List[str], str, float]:
  vals: List[str] = []
  finite: List[Tuple[str, float]] = []

  for precision in precisions:
    row = rows.get(precision)
    rel = row.get('__abs_relerr__', float('nan')) if row else float('nan')
    vals.append(format_float(rel))
    f = to_float(rel)
    if math.isfinite(f):
      finite.append((precision, f))

  if finite:
    best_precision, best_rel = min(finite, key=lambda item: item[1])
    worst_rel = max(v for _, v in finite)
    return vals, f'{best_precision}:{format_float(best_rel)}', worst_rel

  return vals, '-', float('nan')


def print_markdown_table(sorted_keys: List[Tuple[str, str, str]], table: Dict[Tuple[str, str, str], Dict[str, dict]], precisions: List[str]):
  print('| function | operation | domain | ' + ' | '.join(precisions) + ' | best_precision | best_relerr | worst_relerr |')
  print('|---|---|---|' + '|'.join('---' for _ in precisions) + '|---|---|---|')

  for function_name, operation, domain in sorted_keys:
    vals, best, worst = stats_for_precision_rows(table.get((function_name, operation, domain), {}), precisions)
    print('| ' + ' | '.join([function_name, operation, domain] + vals + [best, format_float(worst)]) + ' |')


def print_tsv(sorted_keys: List[Tuple[str, str, str]], table: Dict[Tuple[str, str, str], Dict[str, dict]], precisions: List[str]):
  print('\t'.join(['function', 'operation', 'domain'] + precisions + ['best_precision', 'best_relerr', 'worst_relerr']))

  for function_name, operation, domain in sorted_keys:
    vals, best, worst = stats_for_precision_rows(table.get((function_name, operation, domain), {}), precisions)
    best_precision = '-'
    best_relerr = '-'
    if best != '-':
      best_precision, best_relerr = best.split(':', 1)
    print('\t'.join([function_name, operation, domain] + vals + [best_precision, best_relerr, format_float(worst)]))


def main() -> int:
  parser = argparse.ArgumentParser()
  parser.add_argument('--report', action='append', type=parse_report_path,
                      help='report in the form PREC=PATH (e.g. --report dd=/tmp/oracle_unary_dd_worst.csv)')
  parser.add_argument('--report-dir', help='directory that contains oracle_unary_<precision>_worst.csv files')
  parser.add_argument('--precision', action='append',
                      help='optional precision ordering override; repeatable, e.g. --precision dd --precision td')
  parser.add_argument('--format', choices=['table', 'tsv'], default='table',
                      help='output format (default: table)')
  args = parser.parse_args()

  report_map: Dict[str, str] = {}

  if args.report:
    for name, path in args.report:
      if name in report_map:
        raise SystemExit(f'duplicate report name: {name}')
      report_map[name] = path

  if args.report_dir:
    for name in os.listdir(args.report_dir):
      match = _UNARY_REPORT_RE.match(name)
      if not match:
        continue
      precision = match.group(1)
      if precision not in report_map:
        report_map[precision] = os.path.join(args.report_dir, name)

  if not report_map:
    raise SystemExit('No reports given. Use --report NAME=PATH or --report-dir')

  rows_by_precision: Dict[str, List[dict]] = {}
  for precision, path in sorted(report_map.items()):
    if not os.path.exists(path):
      raise SystemExit(f'report path not found: {path}')
    rows_by_precision[precision] = read_rows(path)

  if args.precision:
    requested = [p for p in args.precision if p]
    missing = [p for p in requested if p not in rows_by_precision]
    if missing:
      raise SystemExit("requested precision not found: %s" % ', '.join(missing))
    rows_by_precision = {p: rows_by_precision[p] for p in requested}
  else:
    ordered = [p for p in UNIFIED_PRECISION_ORDER if p in rows_by_precision]
    remaining = [p for p in sorted(rows_by_precision.keys()) if p not in ordered]
    rows_by_precision = {p: rows_by_precision[p] for p in ordered + remaining}

  sorted_keys, table = build_index(rows_by_precision)
  precisions = list(rows_by_precision.keys())

  if args.format == 'table':
    print_markdown_table(sorted_keys, table, precisions)
  else:
    print_tsv(sorted_keys, table, precisions)

  return 0


if __name__ == '__main__':
  raise SystemExit(main())
