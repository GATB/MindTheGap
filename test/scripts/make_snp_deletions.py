#!/usr/bin/env python

# -*- coding: utf-8 -*-

# base import
import argparse
import sys
import os

# specific import
import random
from collections import defaultdict

# log import
import logging


# logger configuration
def conf_logger():
    """ Set configuraion of root logger """

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(levelname)s : %(message)s')

    steam_handler = logging.StreamHandler()
    steam_handler.setLevel(logging.DEBUG)
    steam_handler.setFormatter(formatter)
    logger.addHandler(steam_handler)


# argparse control argument function
def file_exist(filename):
    """ Check if file exist """

    if not os.path.isfile(filename):
        raise argparse.ArgumentTypeError("we can't access to %s file"
                                         % filename)
    return str(filename)


def unsigned_int(numbre):
    """ Check if numbre is a unsigned integer"""

    inumber = int(numbre)
    if numbre < 0:
        raise argparse.ArgumentTypeError("%s isn't a positive int value"
                                         % numbre)
    return inumber


def file_with_extension_exist(filename):
    """ Check if file with vde or fasta doesn't exist """

    for ext in ["vde", "fasta"]:
        if os.path.isfile(str(filename)+"."+str(ext)):
            raise argparse.ArgumentTypeError(
                "we need %s.%s file not exist" % str(filename), str(ext))

    return str(filename)


# user argument check
def check_interval(first, second):
    """ If intervale isn't valid generate a warning and
     return the good interval"""

    if first > second:
        logging.getLogger().warning(
            "[%d, %d] isn't a valid interval we use [%d, %d]"
            % (first, second, second, first))
        return second, first

    if first == second:
        logging.getLogger().warning(
            "[%d, %d] interval content just one number"
            % (first, second))
        return first, second

    return first, second


# specific function
def pos_near(base_list, pos, min_dist):
    """ Check if pos isn't near a pos in base_list """
    for base_pos in base_list:
        if abs(base_pos - pos) < min_dist:
            return True

    return False


def generate_snp_del(seq, pos_del, pos_snp, del_size):
    """ Create a SNP and a deletion """

    nuc = ['A', 'C', 'T', 'G']

    nuc.remove(seq[pos_snp])
    seq = seq[:pos_snp] + nuc[random.randint(0, 2)] + seq[pos_snp+1:]

    seq = seq[:pos_del] + seq[pos_del+del_size:]


def write_vde(file_handler, pos, type, comment):
    """ Write information in vde format """

    file_handler.write("%s,%s,%s\n" % (pos, type, comment))


def write_seq(file_handler, comment, seq):
    """ Write information in fasta format """
    file_handler.write(">%s\n" % comment)
    file_handler.write("%s\n" % seq)


# Main programme function
def main():
    """ The main function of make_snp_deletions no argument """

    logger = logging.getLogger()
    parser = argparse.ArgumentParser(
        prog="make_snp_deletions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-g", "--genome", type=file_exist,
                        help="fasta file content the genome", required=True)
    parser.add_argument("-o", "--output", type=file_with_extension_exist,
                        help="""prefix of output file .fasta is genome with SNP
                         and deletion and .vde is the postion and type for each
                          'variant""", required=True)
    parser.add_argument("-n", "--number-del", type=unsigned_int,
                        help="number of deletions to generate", default=1)
    parser.add_argument("-m", "--min-size-del", type=unsigned_int,
                        help="minimal size of the deletions (in bp)",
                        default=100)
    parser.add_argument("-M", "--max-size-del", type=unsigned_int,
                        help="maximal size of the deletions (in bp)",
                        default=200)
    parser.add_argument("-s", "--min-dist-snp", type=unsigned_int,
                        help="minimal distance between snp and deletion (in bp)",
                        default=5)
    parser.add_argument("-S", "--max-dist-snp", type=unsigned_int,
                        help="maximal distance between snp and deletion (in bp)",
                        default=31)
    parser.add_argument("-d", "--variant-dist", type=unsigned_int,
                        help="distance minimal between two variant (in bp)",
                        default=232)

    # Parse cli argument
    arg = vars(parser.parse_args())

    # Check interval
    del_size_min, del_size_max = check_interval(
        arg["min_size_del"], arg["max_size_del"])
    dist_snp_min, dist_snp_max = check_interval(
        arg["min_dist_snp"], arg["max_dist_snp"])

    # Check variant distance
    if arg["variant_dist"] <= (arg["max_size_del"] + arg["max_dist_snp"]):
        logging.getLogger().warning(
            "variant distance is minus possible variant max size.")

    comment = ""
    comment2seq = defaultdict(str)
    genome_size = 0
    with open(arg["genome"]) as genome_file:
        for line in genome_file:
            if line.startswith(">"):
                comment = line.lstrip(">").rstrip()
            else:
                genome_size += len(line.rstrip())
                comment2seq[comment] += line.rstrip()

    nuc_per_del = genome_size / arg["number_del"]

    # output file openning
    vde_file = open(arg["output"]+".vde", "a")
    output_file = open(arg["output"]+".fasta", "a")

    # generate snp and deletion
    seq_del_cpt = 0
    del_cpt = 0
    list_pos = list()
    for comment in comment2seq.keys():
        while seq_del_cpt < (len(comment2seq[comment]) / nuc_per_del):

            del_cpt += 1
            seq_del_cpt += 1

            del_pos = random.randint(0, len(comment2seq[comment]))

            while pos_near(list_pos, del_pos, arg["variant_dist"]):
                del_pos = random.randint(0, len(comment2seq[comment]))
            list_pos.append(del_pos)

            snp_pos = del_pos - random.randint(dist_snp_min,
                                               dist_snp_max)

            generate_snp_del(comment2seq[comment], del_pos, snp_pos,
                             random.randint(del_size_min,
                                            del_size_max))

            write_vde(vde_file, snp_pos, "snp", comment)
            write_vde(vde_file, del_pos, "homo", comment)

        seq_del_cpt = 0
        write_seq(output_file, comment, comment2seq[comment])

    # output file closeing
    vde_file.close()
    output_file.close()


if __name__ == "__main__":
    conf_logger()
    main()
