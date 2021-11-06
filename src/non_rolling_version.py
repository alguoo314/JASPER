import sys
import re
import os
import argparse
import math
import csv
import dna_jellyfish as jf
#from timer import Timer

def main(contigs,query_path,k,test,fix,fout,tout,fixedout,database,thre):
    try:
        #t = Timer()
        threshold = thre
        db = database
        if contigs != None and database == None:
            threshold,db = jellyfish(contigs,k,thre)
        if (contigs == None and (database == None or thre == None)) or (contigs != None and database != None): 
            sys.stderr.write("Wrong arguments. One and only one between the contigs and database argument should be given. And if contigs is not given, database and threshold must both be given.")
            return
        print("Threshold =  {}".format(threshold))
        qf  = jf.QueryMerFile(db)
        seq_dict = parse_fasta(query_path)                                                                                                                 
        wrong_kmers_dict = {}
        seqs = []
        fixed_bases_list = []
        total_wrong_kmers = 0
        total_kmers = 0

        for seqname,seq in seq_dict.items():
            #t.start()
            print(seqname+":")
            rare_occurance = 0
            good_before = -1 #index of the last guaranteed good base before the mismatch                                                                                     
            backtracked = False
            i = 0 #first k mer at position 0                                                                                                       
            while i < len(seq)-k+1: #k is 25 : 
                wrong_kmers_list = []                                                                                                   
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
                
                if occurrance < threshold:
                    rare_occurance += 1
                    #count_record.append(occurrance)
                    if backtracked == True:
                        i+=1
                        continue
                    j = i-1
                    backtracked = True
                    occurrance = qf[jf.MerDNA(seq[j:k+j]).get_canonical()] 
                    while occurrance < threshold and j>=0:
                        #count_record.insert(0,occurrance)
                        j = j-1
                        rare_occurance += 1
                        occurrance = qf[jf.MerDNA(seq[j:k+j]).get_canonical()]
                    good_before = j+k-1 #the end base of a good kmer
                    print("previous good kmers:")
                    prev_good_count = qf[jf.MerDNA(seq[j:k+j]).get_canonical()]
                    print(prev_good_count, end = " ")
                    print(jf.MerDNA(seq[j:k+j]).get_canonical()) 
                    kmer_count = qf[jf.MerDNA(seq[i:k+i]).get_canonical()]                                                            
                    #go forward back to i
                    if j == -1: #even the first kmer is bad                                                                                                                                 
                        good_before = -1
                    while kmer_count < threshold:
                        i+=1
                        kmer_count = qf[jf.MerDNA(seq[i:k+i]).get_canonical()]
                        
                    good_after = i #the first base of the first good kmer after the mismatch

                    #A kmer is bad if it is below the threshold AND its count is less than 1/2 of the moving average of rolling_num good k-mers before it.
                    
                    while  qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()] > prev_good_count/2 and good_before-k+1 < good_after:
                        if good_before == -1:
                            break
                        prev_good_count = qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()]
                        print(prev_good_count,end=" ") #testing purpose
                        print(jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical())
                        good_before +=1
                    print() #testing purpose
                    
                    num_below_thres_kmers = good_after-good_before-2+k
                    to_be_fixed = seq[max(0,good_before-k+2):good_after+k-1]
                        
                    wrong_kmers_list.extend([*range(max(0,good_before-k+2),good_after)])
                    print("Bad kmers:")
                    for ind in range(max(0,good_before-k+2),good_after):
                         print(qf[jf.MerDNA(seq[ind:k+ind]).get_canonical()],end = " ") #testing purpose
                    print()
                    print()
                    print("Next iteration")
                    if fix == True:
                        seq,fixed_base = fixing_sid(seq,to_be_fixed,k,threshold,qf,num_below_thres_kmers,good_before,good_after) #fix simple sub/insert/del cases
                        if fixed_base != "nN":
                           fixed_bases_list.append([seqname,good_after-1,fixed_base]) 
                    
                else: #good kmer                                                                                                                   
                    #if rare_occurance >=25: #there is a bad kmer before it                                                                                                                                                                                   
                    backtracked = False
                    good_before = i+k-1 #the end base of a good kmer                                                                                       
                    i += k
                    #rare_occurance = 0
            
                    
            seqs.append(seq)
            wrong_kmers_dict[seqname] = wrong_kmers_list
            total_wrong_kmers += len(wrong_kmers_list)
            total_kmers += len(seq)-k+1


        if test == True:
            p_good = 1-total_wrong_kmers/total_kmers
            e = 1-p_good**(1/k)
            if e != 0:
                Q = round(-10*math.log(e,10),2)
            else:
                Q = "Inf"
            print("Q value = {}, # of bad kmers = {}, error rate = {}, # of total bases in the fasta file = {}".format(Q,total_wrong_kmers,e,total_kmers))
            kmers_fields = ['Contig','Kmer Coord']
            with open(tout,'w') as csvt:
                csvwriter = csv.writer(csvt,delimiter=' ')
                csvwriter.writerow(kmers_fields)
                for seqname in wrong_kmers_dict.keys():
                    for ind in wrong_kmers_dict[seqname]:
                        csvwriter.writerow([seqname,ind+1]) #output index starts with 1 instead of 0
        
        
        if fix == True:
            i = 0
            with open(fixedout,'w') as of:
                for seqname in wrong_kmers_dict.keys():
                    of.write(">{}\n".format(seqname))
                    of.write(seqs[i])
                    of.write("\n")
                    i+=1
            base_fields = ['Contig', 'Base_coord', 'Operation']
            with open(fout,'w') as csvf:
                csvwriter = csv.writer(csvf,delimiter=' ')
                csvwriter.writerow(base_fields)
                csvwriter.writerows(fixed_bases_list)
                
        return
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \                                                                   

         sys.exit(1)



def jellyfish(contigs,k,thre):
    count = math.inf
    threshold = 0
    db_name = os.path.splitext(os.path.basename(contigs[0]))[0]
    contigs = ' '.join(contig_file for contig_file in contigs)
    #print(contigs)
    os.system("jellyfish count -s 300000000 -t 32 -m {} -C -o {} {}".format(k,db_name+".jf",contigs))
    os.system("jellyfish histo -t 32 {}> {}".format(db_name+".jf",db_name+".csv"))
    with open(db_name+".csv",'r') as histo:
        if thre != None:
            return thre, db_name+".csv"
        csvreader = csv.reader(histo,delimiter=' ')
        for row in csvreader:
            if count >= int(row[-1]):
                count = int(row[-1])
                threshold = int(int(row[0])*math.log(2))
            else: #found the local min
                return threshold,db_name+".jf"
                



def fixing_sid(seq,to_be_fixed,k,threshold,qf,num_below_thres_kmers,good_before,good_after):
    fixed_base = "nN"
    if num_below_thres_kmers == k: #substitution or insertion
        b = fix_sub(to_be_fixed,k,threshold,qf)
        if b !=  None:
            #print("Index {} should be {} instead of {}".format(good_after-1,b,seq[good_after-1]))
            seq = seq[:good_after-1] + b + seq[good_after:]
            fixed_base = "s"+b
        else:
            b = fix_insert(to_be_fixed,k,threshold,qf)
            if b != None:
                #print("{} was inserted as index {}, now removed".format(b,good_after-1))
                seq = seq[:good_after-1] + seq[good_after:]
                
                                
    elif num_below_thres_kmers == k-1: #deletion by one or more bases
       removed_bases = fix_del(to_be_fixed,k,threshold,qf)
       if removed_bases != None:
           seq = seq[:good_after]+removed_bases+seq[good_after:]
           #print("{} was lost after index {}".format(removed_bases,good_after-1))
           fixed_base = "d"+removed_bases
       elif (seq[good_before] == seq[good_before+1]) and (fix_same_base_insertion(to_be_fixed,k,threshold,qf) == True):
           #print("{}, the same base as the good base after it, was inserted. Now removed".format(seq[good_after-1]))
           seq = seq[:good_after-1] + seq[good_after:]     
           fixed_base = "i"+seq[good_after-1]
    #else:
       #print("No available fixes for bad kmers from index {} to {}".format(good_before+1,good_after-1))
         
    return seq,fixed_base
                


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
    parser.add_argument("--db", default = None, help="The path to the .jf  database file. Not needed if --contigs is given.")
    parser.add_argument("--contigs",nargs='+',default = None, help="The path to the .fasta file(s） containing the contigs to build the jellyfish database. Not needed if --db is provided")
    parser.add_argument("-q","--query", help = "The path to the .fasta query file")
    parser.add_argument("-thr","--threshold", type=int, default = None, help = "The threshold for a bad kmer.")
    parser.add_argument("-k","--ksize", type=int,help = "The kmer size")
    parser.add_argument("--test", action='store_true',help = "Print loc of bad kmers, total num of bad kmers, and estimate for Q")
    parser.add_argument("--fix", action='store_true', help="Output the index of fixed bases and output the new sequence")
    parser.add_argument("--fout",default = "fout.csv",help = "Output the index of the fixed bases (index based on the original unfix seq)" )
    parser.add_argument("-fo","--fixedout",default = "fixed_seq.fasta",help = "The output file containing the fixed sequence")
    parser.add_argument("--tout", default = "tout.csv", help = "The output file containing the locations of bad kmers")
    args = parser.parse_args()
    main(args.contigs,args.query,args.ksize,args.test,args.fix,args.fout,args.tout,args.fixedout,args.db,args.threshold)
