import sys
import re
import os
import argparse
import dna_jellyfish as jf
from timer import Timer

def main(db_path,query_path):
    try:
        t = Timer()
        qf  = jf.QueryMerFile(db_path)
        seq_dict = parse_fasta(query_path)
        #print("parser ends")                                                  \
                                                                                
        wrong_bases_dict = {}
        for seqname,seq in seq_dict.items():
            t.start()
            rare_occurance = 0
            wrong_bases_list = []
            good = -1 #index of the last guaranteed good base
            backtracked = False
            i = 0 #first k mer at position 0
            while i < len(seq)-25+1: #k is 25 :          
                #make kmers and check occurances                                
                match = re.match("^[ACTGactg]*$",seq[i])
                if match is None:
                    rare_occurance = 0
                    #skip ambigious chars 
                    i += 1                                    
                    continue
                mer_string = seq[i:25+i]
                mer = jf.MerDNA(mer_string).get_canonical()
                occurrance = qf[mer]
                if occurrance <= 0:
                    rare_occurance += 1
                    if backtracked == True:
                        i+=1
                        continue  
                    j = i-1
                    backtracked = True
                    while qf[jf.MerDNA(seq[j:25+j]).get_canonical()] <= 0 and j>=0:
                        j = j-1
                        rare_occurance += 1
                    good = j+25-1 #the end base of a good kmer
                    if good >= i: #if the end base of a good kmer is after i, then the bases between i and good must be good too
                        i = good+1
                        rare_occurance = 0
                        continue
                    if j == -1: #even the first kmer is bad
                        good = -1
                    i +=1                                               
                else: #good kmer
                    if rare_occurance >=25: #there is a bad kmer before it
                        wrong_bases_list.extend([*range(good+1,i)]) #the end of the last bad kmer = bad base
                        #print(seq[])
                    backtracked = False
                    good = i+24 #the end base of a good kmer
                    i += 25
                    rare_occurance = 0
 
                
            t.stop()
            print("put into dict")
            t.start()
            wrong_bases_dict[seqname] = wrong_bases_list
            print(seqname)
            t.stop()
        with open("error_bases_try4.txt",'w') as of:
            for seqname in wrong_bases_dict.keys():
                of.write("{}: ".format(seqname))
                of.write(', '.join(str(ind) for ind in wrong_bases_dict[seqname]) + '\n\n')
        return
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \
                                                                                
         sys.exit(1)

def parse_fasta(query_file):
    f = open(query_file,"r")
    temp = ''
    seq = {}
    name = "placeholder"
    #print("parser running")                                                   \
                                                                                
    for line in f:
        if line.startswith(">"):
            seq[name] = temp
            name = line.split()[0][1:] #save the sequence name/number                           \
                                                                                                 
            temp = ''
            #print(name)                                                                         
        else:
            temp += line.replace('\n','') #remove whitespaces in the sequence                   \
                                                                                                 
            #seq[name] += temp                                                                  \
                                                                                                 
    seq[name] = temp
    seq.pop("placeholder") #remove the first key put into the dict as a placeholder              
    #print("closing the file")                                                                   
    f.close()
    return seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", help="The path to the .jf  database file")
    parser.add_argument("--query", help = "The path to the .fasta query file")
    args = parser.parse_args()
    main(args.db,args.query)
