#!/usr/bin/env python


"""Script to read in dump of Pythia8 settings, and comparing them"""


from __future__ import print_function

import sys
from collections import OrderedDict

def text_to_dict(text):
    settings = OrderedDict()
    for line in text:
        this_line = line.strip()
        if not this_line.startswith("|"):
            continue
        parts = [p.strip() for p in this_line.split("|") if p.strip()]
        if len(parts) < 3:
            continue
        name = parts[0]
        settings[name] = {"now": parts[1], "default": parts[2]}
        # last part can be default min max
        if " " in parts[-1]:
            last_parts = [p.strip() for p in parts[-1].split() if p.strip()]
            settings[name]['default'] = last_parts[0]
            settings[name]['min'] = last_parts[1]
            if len(last_parts) == 3:  # but doesn't always exist
                settings[name]['max'] = last_parts[2]
    return settings


def compare_dicts_common(dict1, dict2, only_diff_values=True):
    """Compare values for keys common to both dicts"""
    keys1 = set(dict1.keys())
    keys2 = set(dict2.keys())

    print("Common settings:")
    if only_diff_values:
        print("(with different values)")
    print("Setting name                 file1                 file2")
    print("-"*80)
    overlap = sorted(list(keys1 & keys2)) # in both
    for k in overlap:
        if only_diff_values and dict1[k]['now'] == dict2[k]['now']:
            continue
        print(k, dict1[k]['now'], dict2[k]['now'])
    print("-"*80)


def compare_dicts_diff(dict1, dict2):
    """Print out settings that are in one dict but not the other"""
    keys1 = set(dict1.keys())
    keys2 = set(dict2.keys())

    print("Only in file1:")
    print("-"*80)
    only_this = sorted(list(keys1 - keys2))
    for k in only_this:
        print(k, dict1[k]['now'])
    print("-"*80)

    print("Only in file2:")
    print("-"*80)
    only_this = sorted(list(keys2 - keys1))
    for k in only_this:
        print(k, dict2[k]['now'])
    print("-"*80)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage:")
        print("%s filename1.txt filename2.txt" % __file__)
        exit()

    filename1 = sys.argv[1]
    dict1 = None
    with open(filename1) as f:
        dict1 = text_to_dict(f)

    filename2 = sys.argv[2]
    dict2 = None
    with open(filename2) as f:
        dict2 = text_to_dict(f)

    compare_dicts_common(dict1, dict2, only_diff_values=True)
    compare_dicts_diff(dict1, dict2)
