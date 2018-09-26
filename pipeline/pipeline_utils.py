import argparse
import sys
import logging

def contig_stats(contigFile):
    file=open(contigFile,'r')
    nbContigs=0
    length=0
    for line in file:
        if line.startswith('>'):
            nbContigs+=1
        else:
            length+=len(line)
    logging.info("Number of contigs : "+str(nbContigs))
    logging.info("Assembled size : "+str(length))


class MtgParser(argparse.ArgumentParser):

    def error(self, message):
        print("")
        sys.stderr.write('error: %s\n' % message)
        print("")
        self.print_help()
        sys.exit(2)


class ArgumentFormatterMtg(argparse.HelpFormatter):


    #def _fill_text(self, text, width, indent):
    #    return ''.join([indent + line for line in text.splitlines(True)])
    def _split_lines(self, text, width):
        return text.splitlines()

    #remove default args layout
    def _format_args(self, action, default_metavar):
        result = ""
        return result

    #Remove "usage: ..." header
    def _format_usage(self, usage, actions, groups, prefix):
        return ""


    #Changed layout of each item
    def _get_help_string(self, action):

        text = ""

        if type(action) == argparse._StoreAction:
            text =  "(1 arg) :    " + action.help
        elif type(action) == argparse._StoreTrueAction:
            text =  "(0 arg) :    " + action.help

        if type(action) == argparse._StoreAction and action.default != None:
            text += " [Default: " + str(action.default) + "]"
        #print type(action), action
        #print action
        #return "-5-"
        #return action.help
        if text != "":
            return text

        return "__none__"

    #Hack for removing useless "optional arguments:" section
    def _join_parts(self, part_strings):
        #print part_strings
        return ''.join([part
                        for part in part_strings
                        if part and part is not argparse.SUPPRESS and not "optional arguments:" in part and not "__none__" in part and not "--help" in part])
