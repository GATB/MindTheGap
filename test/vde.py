#!/usr/bin/env python

# -*- coding: utf-8 -*-

# std import
import sys
import argparse
import os.path
from collections import defaultdict
import csv
import re
import math


class Variant(object):

    def __init__(self, type, comment):
        self.type = type
        self.comment = comment

    def __eq__(self, other):
        if self.type != other.type:
            return False

        if self.comment != other.comment:
            return False

        return True

    def __str__(self):
        return "%s_%s" % (str(self.type), str(self.comment))

    def __repr__(self):
        return "<%s - %s>" % (self.type, self.comment)

    def __hash__(self):
        return hash(self.type + self.comment)


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

        prec = 1 if total == 0 else result[gap]["TP"]/float(total)
        recall = 1 if count[gap] == 0 else result[gap]["TP"]/float(count[gap])

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
        find = False
        if not __pos_in_truth(exp_pos, truth, exp, result):
            for swift in range(delta):
                if __pos_in_truth(exp_pos + swift, truth, exp, result):
                    find = True
                    break
                if __pos_in_truth(exp_pos - swift, truth, exp, result):
                    find = True
                    break
        else:
            find = True

        if not find:
            find = False
            for variant in exp[exp_pos]:
                __iterate_result(result, variant.type, "FP")

    return result


def eva2eva(filename):
    """ Read eva file and return value in dict
    position is key and type is value """

    __check_file_exist(filename)

    data = defaultdict(list)
    count = defaultdict(int)

    with open(filename) as csvfile:
        linereader = csv.reader(csvfile)
        for val in linereader:
            data[int(val[0])].append(Variant(val[1], val[2]))
            count[val[1]] += 1

    return data, count


def breakpoints2eva(filename):
    """ Read breakpoint file and return value in dict
    position is key  and type is value """

    __check_file_exist(filename)

    data = defaultdict(list)
    count = defaultdict(int)

    mtg2eva = {"HOM": "homo",
               "HET": "hete",
               "SNP": "snp",
               "MSNP": "multi_snp",
               "DEL": "deletion",
               "BACKUP": "backup"}

    findpos = re.compile(r'pos_([-\d]+)')
    findtype = re.compile(r'_([a-zA-Z]+)$')
    findcomment = re.compile(r'contig_\d+_(.+)_pos')

    with open(filename) as filehand:
        for line in filehand:
            line = line.strip()
            if line.startswith(">left_contig_"):
                data[int(findpos.search(line).group(1))].append(Variant(
                    mtg2eva[findtype.search(line).group(1)],
                    findcomment.search(line).group(1)))

                count[mtg2eva[findtype.search(line).group(1)]] += 1

    return data, count


def __pos_in_truth(pos, truth, exp, result):
    """If pos is in truth add in result, exp variant with good value"""
    if pos in truth.keys():
        set_exp = set(exp[pos])
        set_tru = set(truth[pos])
        for variant in set_exp & set_tru:
            # SNP exception hardcode
            if not variant in ["snp", "multi_snp"]:
                # Normal variant
                __iterate_result(result, variant.type, "TP")
            else if pos in exp.keys():
                # SNP variant in good pos
                __iterate_result(result, variant.type, "TP")
            else:
                # SNP variant not in good pos
                __iterate_result(result, variant.type, "FP")
            for variant in set_exp & set_exp - set_tru:
            __iterate_result(result, variant.type, "FP")

        return True

    return False


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
