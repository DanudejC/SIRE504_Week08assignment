from bioseq.calculation.SeqCal import *
from bioseq.pattern.SeqPattern import *
from bioseq.manipulation.SeqMan import *
# print("in Main.py")

def argparserLocal():
    from argparse import ArgumentParser
    '''Argument parser for the commands'''
    parser = ArgumentParser(prog='myseq', description='Work with sequence')

    subparsers = parser.add_subparsers(
        title='commands', description='Please choose command below:',
        dest='command'
    )
    subparsers.required = True

    cgc_command = subparsers.add_parser('gcContent', help='Calculates the GC content of the input sequence')
    cgc_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")
    
    atc_command = subparsers.add_parser('atContent', help='Calculates the AT content of the input sequence')
    atc_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")    

    cb_command = subparsers.add_parser('countBase', help='Counts the total number of the specified base in the input sequence')
    cb_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")
    cb_command.add_argument("-b", "--base", type=str, default=None, dest='base',
                             help="Input your specified base here (Can only contain one letter : A/T/C/G)")
    cb_command.add_argument("-r", "--revcomp", action='store_true',
                             help="Reverse complements your input sequence prior to analysis")
    
    cbd_command = subparsers.add_parser('countBasesDict', help='Counts the total number of each base in the input sequence')
    cbd_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")
    cbd_command.add_argument("-r", "--revcomp", action='store_true',
                             help="Reverse complements your input sequence prior to analysis")

    cgs_command = subparsers.add_parser('cpgSearch', help='Search for all CpG sequences in the input sequence')
    cgs_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")
    cgs_command.add_argument("-r", "--revcomp", action='store_true',
                             help="Reverse complements your input sequence prior to analysis")

    ets_command = subparsers.add_parser('enzTargetsScan', help='Scans the input sequence for restriction sites')
    ets_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")
    ets_command.add_argument("-e", "--enz", type=str, default=None, dest='enz',
                             help="Input name of restriction enzyme here to scan for in the input sequence (enzymes available : EcoRI, BamHI, HindIII, AccB2I, AasI, AceI)")
    ets_command.add_argument("-r", "--revcomp", action='store_true',
                             help="Reverse complements your input sequence prior to analysis") 

    rs_command = subparsers.add_parser('reverseSeq', help='Returns the reverse of the input sequence')
    rs_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")

    cs_command = subparsers.add_parser('complementSeq', help='Returns the complement of the input sequence')
    cs_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")

    rcs_command = subparsers.add_parser('reverseComplementSeq', help='Returns the reverse complement of the input sequence')
    rcs_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")

    d2r_command = subparsers.add_parser('dna2rna', help="Transcribes the DNA nucleotide 'T' into the RNA nucleotide 'U' in the input sequence")
    d2r_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")
    d2r_command.add_argument("-r", "--revcomp", action='store_true',
                             help="Reverse complements your input sequence prior to analysis")

    d2p_command = subparsers.add_parser('dna2protein', help="Transcribes and translates the input DNA sequence into a protein sequence")
    d2p_command.add_argument("-s", "--seq", type=str, default=None, dest='seq',
                             help="Input your sequence here for analysis (Can only contain A/T/C/G)")
    d2p_command.add_argument("-r", "--revcomp", action='store_true',
                             help="Reverse complements your input sequence prior to analysis")

    subparsers.add_parser('loadCodons', help="Returns the DNA codon table")          

    # parser.print_help()
    return parser

def test():
    # Input
    # parser = argparserLocal()
    # args = parser.parse_args(["cpgScan","-s","AAATTTCCCGGGCGGGGG"])
    # print(args)
    # print(cpgSearch(args.seq))
    seq = 'ATGGGccGTAGAATTCTTGCaaGCCCGT'
    seq = seq.upper()
    print("Transcription: ", dna2rna(seq))
    print("Transcription-revcomp: ", dna2rna(reverseComplementSeq(seq)))
    print("Translation: ", dna2protein(seq))
    print("Translation-revcomp: ", dna2protein(reverseComplementSeq(seq)))
    print("GC Content:", gcContent(seq))
    print("Count Bases: ", countBasesDict(seq))
    print("Count Bases-revcomp: ", countBasesDict(reverseComplementSeq(seq)))
    print("Search EcoRI: ", enzTargetsScan(seq, 'EcoRI'))
    print("Search EcoRI-revcomp: ", enzTargetsScan(reverseComplementSeq(seq), 'EcoRI'))

    
def main():
    parser = argparserLocal()
    args = parser.parse_args()
    # print(args)
    # print(args.command, args.seq)
    # print("Input",args.seq,"\nGC content =", gcContent(args.seq) )

    # if args.seq == None:
    #     print("------\nError: You did not input your sequence\n------\n")
    # else:
    #     seq = args.seq.upper()
    
    # Input
    # seq = 'ATGGGccGTAGAATTCTTGCaaGCCCGT'

    if args.command == 'gcContent':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['gcContent','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        print("Input :",args.seq.upper(),"\nGC content =", gcContent(seq) )
        
    elif args.command == 'atContent':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['atContent','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        print("Input :",args.seq.upper(),"\nAT content =", atContent(seq) )    

    elif args.command == 'countBase':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['countBase','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        if args.base == None:
            print("------\nError: You did not specify the base to search for\n------\n")
            exit(parser.parse_args(['countBase','-h']))
        elif bool(len(args.base) > 1) == True:
            print("------\nError: You input more than one letter in base search\n------\n")
            exit()
        elif bool(re.search(r'[^ATCGatcg]', args.base)) == True:
            print("------\nError: You input other letters than A/T/C/G in base search\n------\n")
            exit()
        base = args.base.upper()
        if args.revcomp == True:
            seq = complementSeq(reverseSeq(seq))
            print('The input sequence was reverse complemented prior to analysis')
        print("Input :",args.seq.upper(),"\nBase specified :",args.base.upper(), "\nCount =", countBase(seq, base) )

    elif args.command == 'countBasesDict':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['countBasesDict','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        if args.revcomp == True:
            seq = complementSeq(reverseSeq(seq))
            print('The input sequence was reverse complemented prior to analysis')
        print("Input :",args.seq.upper(),"\nAll bases count =", countBasesDict(seq) )    

    elif args.command == 'cpgSearch':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['cpgSearch','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        if args.revcomp == True:
            seq = complementSeq(reverseSeq(seq))
            print('The input sequence was reverse complemented prior to analysis')
        print("Input :",args.seq.upper(),"\nCpG sites found :", cpgSearch(seq) )    

    elif args.command == 'enzTargetsScan':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['enzTargetsScan','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        if args.enz == None:
            print("------\nError: You did not specify a restriction enzyme to scan for\n------\n")
            exit(parser.parse_args(['enzTargetsScan','-h']))
        elif bool(args.enz == "EcoRI" or args.enz == "BamHI" or args.enz == "HindIII" or args.enz == "AccB2I" or args.enz == "AasI" or args.enz == "AceI") == False:
            print("------\nError: Restriction enzyme input is incorrect\n------\n")
            exit()
        if args.revcomp == True:
            seq = complementSeq(reverseSeq(seq))
            print('The input sequence was reverse complemented prior to analysis')
        print("Input :",args.seq.upper(),"\nRestriction enzyme chosen :",args.enz,"\n",args.enz,"sites found :", enzTargetsScan(seq, args.enz ) )

    elif args.command == 'reverseSeq':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['reverseSeq','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        print("Input : ",args.seq.upper(),"\nOutput :", reverseSeq(seq) )    

    elif args.command == 'complementSeq':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['complementSeq','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        print("Input : ",args.seq.upper(),"\nOutput :", complementSeq(seq) )    

    elif args.command == 'reverseComplementSeq':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['reverseComplementSeq','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        print("Input : ",args.seq.upper(),"\nOutput :", reverseComplementSeq(seq) )    

    elif args.command == 'dna2rna':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['dna2rna','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        if args.revcomp == True:
            seq = complementSeq(reverseSeq(seq))
            print('The input sequence was reverse complemented prior to analysis')
        print("Input : ",args.seq.upper(),"\nOutput :", dna2rna(seq) )    

    elif args.command == 'dna2protein':
        if args.seq == None:
            print("------\nError: You did not input your sequence\n------\n")
            exit(parser.parse_args(['dna2protein','-h']))
        elif bool(re.search(r'[^ATCGatcg]', args.seq)) == True:
            print("------\nError: Your sequence contains other letters than A/T/C/G\n------\n")
            exit()
        seq = args.seq.upper()
        if args.revcomp == True:
            seq = complementSeq(reverseSeq(seq))
            print('The input sequence was reverse complemented prior to analysis')
        print("Input : ",args.seq.upper(),"\nOutput :", dna2protein(seq) )    

    elif args.command == 'loadCodons':
        print(loadCodons())

# print(__name__)
if __name__ == "__main__":
    test()
    # main()


