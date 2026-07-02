#!/usr/bin/env python3
import argparse
import csv
import math
import os
from typing import Dict, Tuple


def to_float(text: str) -> float:
  try:
    return float(text)
  except (TypeError, ValueError):
    return float('nan')


def read_reports(report_dir: str):
  data = []
  for name in sorted(os.listdir(report_dir)):
    if not name.endswith('.rounding_corners.csv'):
      continue
    if name.startswith('default.') or name == 'default.rounding_corners.csv':
      # keep default config isolated for CI smoke coverage checks.
      continue
    path = os.path.join(report_dir, name)
    with open(path, newline='') as f:
      rows = csv.DictReader(f)
      required = {
          "build_variant", "type", "input_id", "operation", "relerr",
          "ieee_add", "sloppy_mul", "sloppy_div", "fma", "pass", "allowed",
          "seed"
      }
      if not required.issubset(set(rows.fieldnames or ())):
        raise RuntimeError(f"Missing columns for {path}. expected={sorted(required)} columns={rows.fieldnames}")
      for row in rows:
        data.append(row)
  return data


def normalized_input_id(row):
  input_id = row['input_id']
  if row['operation'] == 'add' and input_id.startswith(
      'variant_distinguishing_carry_ripple_cancellation_'):
    return 'variant_distinguishing_carry_ripple_cancellation'
  if row['operation'] == 'mul' and input_id.startswith(
      'qd_mul_variant_distinguishing_'):
    return 'qd_mul_variant_distinguishing'
  return input_id


def max_by_key(rows):
  add_map: Dict[Tuple[str, str, str, str, str, str], Dict[str, float]] = {}
  mul_map: Dict[Tuple[str, str, str, str, str, str], Dict[str, float]] = {}

  for row in rows:
    t = row['type']
    op = row['operation']
    input_id = normalized_input_id(row)
    ieee_add = row['ieee_add']
    sloppy_mul = row['sloppy_mul']
    sloppy_div = row['sloppy_div']
    fma = row['fma']
    relerr = abs(to_float(row['relerr']))

    if op == 'add':
      key = (t, input_id, sloppy_mul, sloppy_div, fma)
      if key not in add_map:
        add_map[key] = {}
      prev = add_map[key].get(ieee_add)
      if prev is None or relerr > prev:
        add_map[key][ieee_add] = relerr

    if op == 'mul' and t == 'qd':
      key = (input_id, sloppy_div, fma, ieee_add)
      if key not in mul_map:
        mul_map[key] = {}
      prev = mul_map[key].get(sloppy_mul)
      if prev is None or relerr > prev:
        mul_map[key][sloppy_mul] = relerr

  return add_map, mul_map


def compare_add_variants(add_map: Dict[Tuple[str, str, str, str, str], Dict[str, float]],
                        tol: float):
  ok = True
  mismatches = 0
  missing = 0

  for key, variants in sorted(add_map.items()):
    t, input_id, sloppy_mul, sloppy_div, fma = key
    if '0' not in variants or '1' not in variants:
      missing += 1
      ok = False
      continue
    sloppy = variants['0']
    ieee = variants['1']
    if math.isnan(ieee) or math.isnan(sloppy) or math.isinf(ieee) or math.isinf(sloppy):
      ok = False
      mismatches += 1
      continue
    if ieee > sloppy + tol:
      mismatches += 1
      ok = False
      print('add mismatch:', 'type=%s input_id=%s sloppy_mul=%s sloppy_div=%s fma=%s '
            'ieee=%.17g sloppy=%.17g' %
            (t, input_id, sloppy_mul, sloppy_div, fma, ieee, sloppy))

  print(f'add comparisons: {len(add_map)} total, {mismatches} violations, {missing} missing')
  return ok


def compare_qd_mul_variants(mul_map: Dict[Tuple[str, str, str, str], Dict[str, float]],
                           tol: float):
  ok = True
  mismatches = 0
  missing = 0

  for key, variants in sorted(mul_map.items()):
    input_id, sloppy_div, fma, ieee_add = key
    if '0' not in variants or '1' not in variants:
      missing += 1
      ok = False
      continue
    sloppy = variants['1']
    accurate = variants['0']
    if math.isnan(accurate) or math.isnan(sloppy) or math.isinf(accurate) or math.isinf(sloppy):
      ok = False
      mismatches += 1
      continue
    if accurate > sloppy + tol:
      mismatches += 1
      ok = False
      print('qd mul mismatch:',
            'input_id=%s sloppy_div=%s fma=%s ieee_add=%s '
            'accurate=%.17g sloppy=%.17g' %
            (input_id, sloppy_div, fma, ieee_add, accurate, sloppy))

  print(f'qd mul comparisons: {len(mul_map)} total, {mismatches} violations, {missing} missing')
  return ok


def main() -> int:
  parser = argparse.ArgumentParser()
  parser.add_argument('report_dir')
  parser.add_argument('--tolerance', type=float, default=0.0,
                      help='Tolerance to relax strict comparison (default: 0)')
  args = parser.parse_args()

  rows = read_reports(args.report_dir)
  if not rows:
    print(f'No rounding-corner reports found in {args.report_dir}')
    return 1

  add_map, mul_map = max_by_key(rows)

  add_ok = compare_add_variants(add_map, args.tolerance)
  mul_ok = compare_qd_mul_variants(mul_map, args.tolerance)

  return 0 if (add_ok and mul_ok) else 2


if __name__ == '__main__':
  raise SystemExit(main())
