import sys
import re
import os
import argparse
import math
import dna_jellyfish as jf
#from timer import Timer

def main(db_path,query_path,threshold,k,test,fix,fout,tout):
    try:
        #t = Timer()
        qf  = jf.QueryMerFile(db_path)
        seq_dict = parse_fasta(query_path)                                                                                                                 
        fixed_bases_dict = {}
        wrong_kmers_dict = {}
        seqs = []
        for seqname,seq in seq_dict.items():
            #t.start()
            print(seqname+":")
            rare_occurance = 0
            fixed_bases_list = []
            wrong_kmers_list = []
            good_before = -1 #index of the last guaranteed good base before the mismatch                                                                                     
            backtracked = False
            i = 0 #first k mer at position 0                                                                                                       
            while i < len(seq)-k+1: #k is 25 :                                                                                                    
                #make kmers and check occurances   
                match = re.match("^[ACTGactg]*$",seq[i])
                if match is None:
                    rare_occurance = 0
                    #skip ambigious chars                                                                                                          
                    i += 1
                    continue
                mer_string = seq[i:k+i]
                mer = jf.MerDNA(mer_string).get_canonical()
                occurrance = qf[mer]
                if occurrance <= threshold:
                    rare_occurance += 1
                    if backtracked == True:
                        i+=1
                        continue
                    j = i-1
                    backtracked = True
                    while qf[jf.MerDNA(seq[j:k+j]).get_canonical()] <= threshold and j>=0:
                        j = j-1
                        rare_occurance += 1
                    good_before = j+k-1 #the end base of a good kmer  
                    if j == -1: #even the first kmer is bad                                                                                        
                        good_before = -1 
                                                                                  
                    i = good_before +1 #skip to the first possible bad base
                    while qf[jf.MerDNA(seq[i:k+i]).get_canonical()] <= threshold:
                        i+=1
                    good_after = i #the first base of the first good kmer after the mismatch
                    num_bad_kmers = good_after-good_before-2+k
                    to_be_fixed = seq[max(0,good_before-k+2):good_after+k-1]
                    wrong_kmers_list.extend([*range(max(0,good_before-k+2),good_after)])
                    #print(good_before+2-k,good_after-k)
                    if fix == True:
                        seq = fixing(seq,to_be_fixed,k,threshold,qf,num_bad_kmers,good_before,good_after)
                        #seqs.append(seq)
                    
                else: #good kmer                                                                                                                   
                    #if rare_occurance >=25: #there is a bad kmer before it                                                                                                                                                                                   
                    backtracked = False
                    good_before = i+24 #the end base of a good kmer                                                                                       
                    i += 25
                    #rare_occurance = 0
            
                    
            seqs.append(seq)
            wrong_kmers_dict[seqname] = wrong_kmers_list
            if test == True:
                p_good = 1-len(wrong_kmers_list)/len(seq)
                e = 1-p_good**(1/k)
                print("Bad kmers:{}".format(wrong_kmers_list))
                print("Total num of bad kmers: "+str(len(wrong_kmers_list)))
                if e != 0:
                    Q = -10*math.log(e,10)
                else:
                    Q = "Inf"
                #print("{}: Q = {}".format(seqname,Q))
                print("Q = {}\n".format(Q))

        if test == True:
            with open(tout,'w') as of:
                for seqname in wrong_kmers_dict.keys():
                    of.write("{}: ".format(seqname))
                    of.write(', '.join(str(ind+1) for ind in wrong_kmers_dict[seqname]) + '\n\n') #output index starts with 1 instead of 0
        
        
        if fix == True:
            i = 0
            with open(fout,'w') as of:
                for seqname in wrong_kmers_dict.keys():
                    of.write("{}: ".format(seqname))
                    of.write(seqs[i])
                    of.write("\n")
                    i+=1
        return
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \                                                                   

         sys.exit(1)


def fixing(seq,to_be_fixed,k,threshold,qf,num_bad_kmers,good_before,good_after):
    if num_bad_kmers == k: #substitution or insertion
        b = fix_sub(to_be_fixed,k,threshold,qf)
        if b !=  None:
            print("Index {} should be {} instead of {}".format(good_after-1,b,seq[good_after-1]))
            seq = seq[:good_after-1] + b + seq[good_after:]
        else:
            b = fix_insert(to_be_fixed,k,threshold,qf)
            if b != None:
                print("{} was inserted as index {}, now removed".format(b,good_after-1))
                seq = seq[:good_after-1] + seq[good_after:]
                                
    elif num_bad_kmers == k-1: #deletion by one or more bases
       removed_bases = fix_del(to_be_fixed,k,threshold,qf)
       if removed_bases != None:
           seq = seq[:good_after]+removed_bases+seq[good_after:]
           print("{} was lost after index {}".format(removed_bases,good_after-1))
       elif (seq[good_before] == seq[good_before+1]) and (fix_same_base_insertion(to_be_fixed,k,threshold,qf) == True):
           print("{}, the same base as the good base after it, was inserted. Now removed".format(seq[good_after-1]))
           seq = seq[:good_after-1] + seq[good_after:] 
                            
    else:
       print("No available fixes for bad kmers from index {} to {}".format(good_before+1,good_after-1))
         
    return seq
                


def fix_sub(seq_to_be_fixed,k,threshold,qf):
    bad_base = seq_to_be_fixed[k-1]
    for b in 'ACTG':
        trial = seq_to_be_fixed
        if b == bad_base:
            continue
        else:
            trial = trial[:k-1]+b+trial[k:]
            fixed = True
            for i in range(len(trial)-k+1):
                if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
                    fixed  = False
            if fixed == True:
                return b
    return None


def fix_insert(seq_to_be_fixed,k,threshold,qf):
    #original_seq = seq_to_be_fixed
    ind_to_be_removed = k-1
    base_to_be_removed = seq_to_be_fixed[ind_to_be_removed]
    seq_to_be_fixed = seq_to_be_fixed[:ind_to_be_removed] + seq_to_be_fixed[ind_to_be_removed+1:]
    fixed = True
    for i in range(len(seq_to_be_fixed)-k+1):
        if qf[jf.MerDNA(seq_to_be_fixed[i:k+i]).get_canonical()] < threshold:
            fixed  = False
    if fixed == True:
            return base_to_be_removed
    return None



def fix_del(seq_to_be_fixed,k,threshold,qf): #assume at most three deletions
    for x in 'NACTG':
        for y in 'NACTG':
            for z in 'NACTG':
                trial = seq_to_be_fixed
                added_bases = (x+y+z).replace("N", "")
                trial = seq_to_be_fixed[:k-1]+added_bases+seq_to_be_fixed[k-1:]
                fixed = True
                for i in  range(len(trial)-k+1):
                    if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
                        fixed  = False
                        break
                if fixed == True:
                    return added_bases
    return None



def fix_same_base_insertion(seq_to_be_fixed,k,threshold,qf):
    ind_to_be_removed = k-1
    seq_to_be_fixed = seq_to_be_fixed[:ind_to_be_removed] + seq_to_be_fixed[ind_to_be_removed+1:]
    fixed=True
    for i in  range(len(seq_to_be_fixed)-k+1):
        if qf[jf.MerDNA(seq_to_be_fixed[i:k+i]).get_canonical()] < threshold:
            fixed  = False
    if fixed == True:
        return True
    return None



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
    parser.add_argument("-q","--query", help = "The path to the .fasta query file")
    parser.add_argument("-thr","--threshold", type=int,help = "The threshold for a bad kmer")
    parser.add_argument("-k","--ksize", type=int,help = "The kmer size")
    parser.add_argument("--test", action='store_true',help = "Print loc of bad kmers, total num of bad kmers, and estimate for Q")
    parser.add_argument("--fix", action='store_true', help="Output the index of fixed bases and output the new sequence")
    parser.add_argument("--fout",default = "fout.fasta", help = "The output file containing the fixed sequence")
    parser.add_argument("--tout", default = "fout.fasta", help = "The output file containing the locations of bad kmers")
    args = parser.parse_args()
    main(args.db,args.query,args.threshold,args.ksize,args.test,args.fix,args.fout,args.tout)
