#!/usr/bin/env python

# -*- coding: utf-8 -*-

# std import
import sys
import argparse
import os.path
from collections import defaultdict
import csv
import re


def main():
    """ The main function of vde no argument """

    enable_input = [name.split("2")[0] for name in globals().keys()
                    if name.endswith("2eva")]

    parser = argparse.ArgumentParser(
        prog="vde",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-e", "--experiment",
                        type=str,
                        help="File of experimente result.",
                        required=True)
    parser.add_argument("-t", "--truth",
                        type=str,
                        help="File of truth result.",
                        required=True)
    parser.add_argument("-d", "--delta",
                        type=int,
                        help="Acceptable diff betwen truth and experimente.",
                        default=5)
    parser.add_argument("-ef", "--experiment-format",
                        type=str,
                        help="Format of experiment file",
                        choices=enable_input,
                        default="eva")
    parser.add_argument("-tf", "--truth-format",
                        type=str,
                        help="Format of truth file",
                        choices=enable_input,
                        default="eva")

    # parsing cli argument
    argument = vars(parser.parse_args())

    expfunc = globals()[argument["experiment_format"]+"2eva"]
    truthfunc = globals()[argument["truth_format"]+"2eva"]

    experiment, count = expfunc(argument["experiment"])
    truth, count = truthfunc(argument["truth"])

    result = compare(experiment, truth, argument["delta"])

    result_printing(result, count)


def result_printing(result, count):
    """ Printing the result in csv format """

    head = ",".join(("type", "TP", "FP", "recall", "precision"))
    print(head)
    for gap in result.keys():
        total = result[gap]["TP"] + result[gap]["FP"]

        recall = 1 if total == 0 else result[gap]["TP"]/float(total)
        prec = 1 if count[gap] == 0 else result[gap]["TP"]/float(count[gap])

        print(",".join((str(gap),
                        str(result[gap]["TP"]),
                        str(result[gap]["FP"]),
                        str(recall),
                        str(prec))))


def compare(exp, truth, delta):
    """ Compare experimente and truth return TP FP precision and recall
    for each type """

    result = defaultdict(lambda: defaultdict(int))

    for exp_pos in exp.keys():
        for type_gap in exp[exp_pos]:
            find = False
            if exp_pos in truth and type_gap in truth[exp_pos]:
                __iterate_result(result, type_gap, "TP")
                find = True
            else:
                for i in range(1, delta+1):
                    prev_pos = str(int(exp_pos) - i)
                    next_pos = str(int(exp_pos) + i)

                    if prev_pos in truth and type_gap in truth[prev_pos]:
                        __iterate_result(result, type_gap, "TP")
                        find = True
                        break
                    if next_pos in truth and type_gap in truth[next_pos]:
                        __iterate_result(result, type_gap, "TP")
                        find = True
                        break

            if not find:
                __iterate_result(result, type_gap, "FP")

    return result


def eva2eva(filename):
    """ Read eva file and return value in dict
    position is key and type is value """

    __check_file_exist(filename)

    data = defaultdict(set)
    count = defaultdict(int)

    with open(filename) as csvfile:
        linereader = csv.reader(csvfile)
        {__add_in_data_count(val[0], val[1], data, count)
         for val in linereader}

    return data, count


def breakpoints2eva(filename):
    """ Read breakpoint file and return value in dict
    position is key  and type is value """

    __check_file_exist(filename)

    data = defaultdict(set)
    count = defaultdict(int)

    mtg2eva = {"HOM": "homo",
               "HET": "hete",
               "SNP": "snp",
               "MSNP": "multi_snp",
               "DEL": "deletion",
               "BACKUP": "backup"}

    findpos = re.compile(r'pos_(\d+)')
    findtype = re.compile(r'_([a-zA-Z]+)$')

    with open(filename) as filehand:
        for line in filehand:
            line = line.strip()
            if line.startswith(">"):
                __add_in_data_count(findpos.search(line).group(1),
                                    mtg2eva[findtype.search(line).group(1)],
                                    data, count)

    return data, count


def __iterate_result(result, type_gap, tpofp):
    """ If key is in dict iterate this else init this. """

    if type_gap in result.keys():
        result[type_gap][tpofp] += 1
    else:
        result[type_gap][tpofp] = 1


def __add_in_data_count(pos, type_gap, data, counter):
    """ Add value pos: type_gap in data and increment counter[data] """

    data[pos].add(type_gap)
    counter[type_gap] += 1


def __check_file_exist(filename):
    """ If file doesn't exist trow assert """

    assert os.path.isfile(filename), "Error when I try open " + filename


if(__name__ == '__main__'):
    main()
