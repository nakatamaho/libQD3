#!/usr/bin/env python3
"""Reformat old Oracle TAP logs to current diagnostics display policy."""

import argparse
import os
import re
import tempfile


def display_digits_for_test(name):
  if name.startswith("dd"):
    return 33
  if name.startswith("td"):
    return 49
  if name.startswith("qd"):
    return 65
  if name.startswith("edd"):
    return 39
  return 80


def should_shorten_value(name):
  return (
      name == "mpfr_reference"
      or name == "got_value"
      or name == "input_value"
      or name == "input_a_value"
      or name == "input_b_value"
      or name == "mpfr_sin"
      or name == "mpfr_cos"
      or name == "mpfr_tan"
      or name == "abs_error_mpfr"
      or name == "td_reduced_arg"
      or name == "td_reduced_ref"
  )


def trim_mantissa_digits(text, max_significant_digits):
  if max_significant_digits <= 0 or not text:
    return text

  e_pos = text.find("e")
  if e_pos < 0:
    e_pos = text.find("E")
  if e_pos < 0:
    return text

  dot_pos = text.find(".")
  if dot_pos < 0 or dot_pos >= e_pos:
    return text

  sign_pos = 1 if (text.startswith("-") or text.startswith("+")) else 0
  if sign_pos >= e_pos:
    return text

  current_digits = e_pos - sign_pos
  if current_digits <= max_significant_digits:
    return text

  max_fraction_digits = max_significant_digits - 1
  if max_fraction_digits <= 0:
    return text[:dot_pos] + text[e_pos:]

  max_len = dot_pos + 1 + max_fraction_digits
  if max_len >= e_pos:
    return text

  return text[:max_len] + text[e_pos:]


def yaml_unquote(text):
  if not (text.startswith("'") and text.endswith("'") and len(text) >= 2):
    return text
  return text[1:-1].replace("''", "'")


def yaml_quote(text):
  return "'" + text.replace("'", "''") + "'"


def reformat_lines(lines):
  test_name = ""
  digits = 80
  diag_mode = False

  test_name_re = re.compile(r"^(?:ok|not ok)\\s+\\d+\\s+-\\s+(.*)$")
  diag_re = re.compile(r"^  (\\S+): (.*)$")

  for line in lines:
    stripped = line.rstrip("\n")
    match = test_name_re.match(stripped)
    if match:
      test_name = match.group(1)
      digits = display_digits_for_test(test_name)
      diag_mode = False
      yield line
      continue

    if stripped == "  ---":
      diag_mode = True
      yield line
      continue
    if stripped == "  ...":
      diag_mode = False
      yield line
      continue

    if diag_mode:
      diag_match = diag_re.match(stripped)
      if diag_match:
        key = diag_match.group(1)
        value = diag_match.group(2)
        value = yaml_unquote(value)
        if should_shorten_value(key):
          value = trim_mantissa_digits(value, digits)
        prefix = "  " + key + ": "
        if key == "got_value":
          prefix += "     "
        yield f"{prefix}{yaml_quote(value)}\n"
        continue

    yield line


def process_file(path, in_place):
  with open(path, "r", encoding="utf-8") as inp:
    lines = list(inp)
  out_lines = reformat_lines(lines)

  if not in_place:
    for line in out_lines:
      print(line, end="")
    return

  fd, tmp = tempfile.mkstemp(prefix=".tmp_oracle_tap_", dir=os.path.dirname(path) or ".")
  with os.fdopen(fd, "w", encoding="utf-8") as out:
    for line in out_lines:
      out.write(line)
  os.replace(tmp, path)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("files", nargs="*")
  parser.add_argument("--in-place", action="store_true", help="Rewrite files in place")
  args = parser.parse_args()

  if args.in_place and not args.files:
    parser.error("--in-place requires one or more file arguments")

  if args.files:
    for path in args.files:
      process_file(path, args.in_place)
  else:
    for line in reformat_lines(__import__("sys").stdin):
      print(line, end="")


if __name__ == "__main__":
  main()
