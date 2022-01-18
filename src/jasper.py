import sys
import re
import os
import argparse
import math
import csv
import dna_jellyfish as jf
import textwrap

def main(contigs,query_path,k,test,fix,fout,tout,fixedout,database,thre,rep_thre):
    try:
        threshold = thre
        rep_region_threshold = rep_thre
        db = database
        if ((contigs == None and database == None) or (contigs != None and database != None)): 
            sys.stderr.write("Wrong arguments. One and only one between the contigs and database argument should be given. ")
            return
        
        threshold,rep_region_threshold,db = jellyfish(contigs, database, k,thre,rep_thre)
        
        print("Threshold =  {}".format(threshold))
        print("Threshold for Repetitive Region=  {}".format(rep_region_threshold))
        qf  = jf.QueryMerFile(db)
        seq_dict = parse_fasta(query_path)
        wrong_kmers_dict = {}
        wrong_kmers_list = []
        seqs = []
        fixed_bases_list = []
        total_wrong_kmers = 0
        total_kmers = 0

        for seqname,seq in seq_dict.items():
            total_kmers += len(seq)-k+1
            print(seqname+":")
            rare_occurance = 0
            good_before = -1 #index of the last guaranteed good base before the mismatch                                                                                     
            backtracked = False
            i = 0 #first k mer at position 0
            wrong_kmers_list = []                                                                                                       
            while i < len(seq)-k+1: #k is 25 :  
                mer_string = seq[i:k+i]
                #match = re.match("^[ACTGactg]*$",mer_string)
                #if match is None:
                N = mer_string.find('N')
                if N >= 0:
                    i+=(N+1)
                    rare_occurance = 0
                    continue
                    
                mer = jf.MerDNA(mer_string).get_canonical()
                occurrance = qf[mer]
                
                if occurrance < threshold:
                    rare_occurance += 1
                    if backtracked == True:
                        i+=1
                        continue
                    j = i-1
                    backtracked = True
                    occurrance = qf[jf.MerDNA(seq[j:k+j]).get_canonical()] 
                    while occurrance < threshold and j>=0:
                        j = j-1
                        rare_occurance += 1
                        occurrance = qf[jf.MerDNA(seq[j:k+j]).get_canonical()]
                    good_before = j+k-1 #the end base of a good kmer
                    prev_good_count = qf[jf.MerDNA(seq[j:k+j]).get_canonical()]
                    kmer_count = qf[jf.MerDNA(seq[i:k+i]).get_canonical()]                                                            
                    #go forward back to i

                    if j == -1: #even the first kmer is bad                                                                                                                                 
                        good_before = -1
                    while kmer_count < threshold and i < len(seq)-k+1:
                        i+=1
                        kmer_count = qf[jf.MerDNA(seq[i:k+i]).get_canonical()]
                        
                    good_after = i #the first base of the first good kmer after the mismatch

                    #A kmer is bad if it is below the threshold AND its count is less than 1/2 of the good k-mers before it AND its previous good kmer is not from a repetitive region (special case apply).
                    if prev_good_count > rep_region_threshold:
                        #repetitive region
                        #i=good_after
                        #continue
                        #jan 14 revision, directly consider it as diploid case
                        if len([*range(max(0,good_before-k+2),good_after)]) < k-1:
                            shortcut_flag = True
                        else:
                            continue

                    else:  #nonrepetitive region case
                        shortcut_flag = False
                        while  qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()] >= prev_good_count/2 and good_before-k+1 < good_after:
                            if good_before == -1:
                                break
                            prev_good_count = qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()]
                            good_before +=1
                        
                        if good_before >= len(seq)-1:
                            break
                    to_be_fixed = seq[max(0,good_before-k+2):good_after+k-1]
                        
                    wrong_kmers_list.extend([*range(max(0,good_before-k+2),good_after)])
                    
                    if fix == True:
                        seq,fixed_base,original,fixed_ind = fixing_sid(shortcut_flag,seq,to_be_fixed,k,threshold,qf,len([*range(max(0,good_before-k+2),good_after)]),good_before,good_after) #fix simple sub/insert/del cases
                        if fixed_base != "nN":
                            if len(fixed_ind) == 1:
                                fixed_bases_list.append([seqname,fixed_ind[0],fixed_base,original])
                            else:
                                fixed_bases_list.append([seqname,fixed_ind[0],fixed_base[0],original[0]])
                                fixed_bases_list.append([seqname,fixed_ind[1],fixed_base[1],original[1]])
                    
                else: #good kmer                                                                                                                                                      
                    backtracked = False
                    good_before = i+k-1 #the end base of a good kmer                                                                                       
                    i += k
                    #rare_occurance = 0
            
                    
            seqs.append(seq)
            wrong_kmers_dict[seqname] = wrong_kmers_list
            total_wrong_kmers += len(wrong_kmers_list)


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
        
        #print(one)
        if fix == True:
            base_fields = ['Contig', 'Base_coord', 'Original',"Mutated"]
            with open(fout,'w') as csvf:
                csvwriter = csv.writer(csvf,delimiter=' ')
                csvwriter.writerow(base_fields)
                csvwriter.writerows(fixed_bases_list)

            i = 0
            with open(fixedout,'w') as of:
                for seqname in wrong_kmers_dict.keys():
                    of.write(">{}\n".format(seqname))
                    new_seq = split_output(seqs[i],60)
                    for l in new_seq:
                        of.write(l+"\n")
                    i+=1
        return
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \                                                                   

         sys.exit(1)




def split_output(seq, num_per_line=60): #make a new line after num_per_line bases
    lines = math.ceil(len(seq)/num_per_line)
    output = []
    for i in range(lines):
        output.append(seq[num_per_line*i:num_per_line*(i+1)])
    return output

         

def jellyfish(contigs,database,k,thre,rep_thre):
    count = math.inf
    threshold = 0
    if database != None:
        db_name = database
        base_name = os.path.splitext(db_name)[0]
    if contigs != None:
        db_name = os.path.splitext(os.path.basename(contigs[0]))[0]+".jf"
        base_name = os.path.splitext(db_name)[0]
        contigs = ' '.join(contig_file for contig_file in contigs)
        #print(contigs)
        os.system("jellyfish count -s 300000000 -t 32 -m {} -C -o {} {}".format(k,db_name,contigs))
    os.system("jellyfish histo -t 32 {}> {}".format(db_name,base_name+".csv"))
    found_thres = False
    with open(base_name+".csv",'r') as histo:
        if thre != None and rep_thre!=None:
            return thre,rep_thre,db_name
        csvreader = csv.reader(histo,delimiter=' ')
        for row in csvreader:
            if count >= int(row[-1]) and found_thres == False:
                count = int(row[-1])
                threshold = int(int(row[0])/2)
            elif found_thres == False: #found the local min, start to find next local max
                found_thres = True
                count = int(row[-1])
            elif count >= int(row[-1]) and found_thres == True:
                repetitive_seq_thr = int(int(row[0])*2)
                if thre == None and rep_thre!=None:
                    return threshold,rep_thre,db_name
                elif thre != None and rep_thre==None:
                    return thre,repetitive_seq_thr,db_name
                elif thre == None and rep_thre==None:
                    return threshold,repetitive_seq_thr,db_name
                else:
                    print("Error")
                    return(None)
            else:
                count = int(row[-1])
                



def fixing_sid(shortcut_flag,seq,to_be_fixed,k,threshold,qf,num_below_thres_kmers,good_before,good_after):
    try:
        fixed_base = "nN"
        original = '-'
        fixed_ind = None

        if num_below_thres_kmers == k and shortcut_flag == False: #substitution or insertion
            b = fix_sub(to_be_fixed,k,threshold,qf)
            if b !=  None:
                original = "s"+seq[good_after-1]
                temp  = seq[:good_after-1] + b + seq[good_after:]                                                           
                seq = temp
                fixed_base = b
                fixed_ind = [good_after-1]
            else:
                b = fix_insert(to_be_fixed,k,threshold,qf)
                if b != None:
                    original = "i"+seq[good_after-1]
                    temp = seq[:good_after-1] + seq[good_after:]
                    seq = temp
                    fixed_base = '-'
                    fixed_ind = [good_after-2]
                                    
        if num_below_thres_kmers == k-1 and shortcut_flag == False: 
           removed_base = fix_del(to_be_fixed,k,threshold,qf)
           if removed_base != None: #deletion of a base
               original = "d-"
               fixed_ind = [good_after]
               temp = seq[:good_after]+removed_base+seq[good_after:]
               seq = temp
               fixed_base = removed_base 
           elif (seq[good_before] == seq[good_before+1]) and (fix_same_base_insertion(to_be_fixed,k,threshold,qf) == True):
               original = "i"+seq[good_before]
               temp = seq[:good_before] + seq[good_before+1:]     
               seq = temp
               fixed_base = "-"
               fixed_ind = [good_before] #insertion after this index
        if num_below_thres_kmers < k-1 and len(to_be_fixed)>k:#skip the good_before = -1 base. Repetitive case is addressed here (ie shortcut is true)
            start,end,s_or_e =  fixhetero(to_be_fixed,k,threshold,qf)
            if s_or_e !=  None: #diploidy
                if s_or_e == "b": #b stands for both bases are changed
                    original = ["s"+seq[good_after-1],"s"+seq[good_before+1]]
                    fixed_base = [str(start),str(end)]
                    fixed_ind = [good_after-1,good_before+1]
                elif s_or_e == "s": #starting base is changed
                    original = "s"+seq[good_after-1]
                    fixed_base = str(start)
                    fixed_ind = [good_after-1]
                else:#the ending base is changed
                    original = "s"+seq[good_before+1]
                    fixed_base = str(end)
                    fixed_ind = [good_before+1]
                temp = seq[:good_after-1]+start+seq[good_after:good_before+1] + end + seq[good_before+2:]
                seq =  temp
            else: #check same base deletion
                removed_base = fix_sb_del(to_be_fixed,k,threshold,qf,num_below_thres_kmers)
                if removed_base != None: #deletion of a base
                    original = "d-"
                    fixed_ind = [good_after]
                    temp = seq[:good_after]+removed_base+seq[good_after:]
                    seq = temp
                    fixed_base = removed_base 
        return seq,fixed_base, original, fixed_ind
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \                                                                   
         sys.exit(1)         
            

        

         
                
def fixhetero(seq_to_be_fixed,k,threshold,qf):
    try:
        left_bad = seq_to_be_fixed[len(seq_to_be_fixed)-k]
        right_bad = seq_to_be_fixed[k-1]
        start = left_bad
        end = right_bad
        for x in 'ACTG':
            for y in 'ACTG':
                if x ==left_bad and y==right_bad:
                    continue
                trial = seq_to_be_fixed[:len(seq_to_be_fixed)-k]+x+seq_to_be_fixed[len(seq_to_be_fixed)-k+1:k-1]+y+seq_to_be_fixed[k:]
                fixed = True
                for i in  range(len(trial)-k+1):
                    if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
                        fixed  = False
                        break
                if fixed == True:
                    start = x
                    end = y
                    if x == left_bad:
                        s_or_e = "e" 
                        
                    elif y == right_bad:
                        s_or_e = "s"
                        
                    else: #both changed
                        s_or_e = "b"
                        
     
                    return(start,end,s_or_e)
        return start,end,None
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \                                                                   
         sys.exit(1)    
    


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



def fix_del(seq_to_be_fixed,k,threshold,qf): 
    for alt in 'ATCG':
        trial = seq_to_be_fixed[:k-1]+alt+seq_to_be_fixed[k-1:]
        fixed = True
        for i in  range(len(trial)-k+1):
            if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
                fixed  = False
                break
        if fixed == True:
            return alt
                        
    return None

def fix_sb_del(seq_to_be_fixed,k,threshold,qf,num_below_thres_kmers):
    sb = seq_to_be_fixed[k-2] #sb stands for same base
    rr = k-1-num_below_thres_kmers #rr stands for remaining repeat (base)
    fixed = False
    if seq_to_be_fixed[k-1-rr:k-1] == sb * rr:
        trial = seq_to_be_fixed[:k-1]+sb+seq_to_be_fixed[k-1:]
        fixed = True
        for i in  range(len(trial)-k+1):
            if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
                fixed  = False
                break
    if fixed == True:
        return sb
    else:
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
    for line in f:
        if line.startswith(">"):
            seq[name] = temp
            name = line.split()[0][1:] #save the sequence name/number
            temp = ''
        else:
            temp += line.replace('\n','') #remove whitespaces in the sequence

    seq[name] = temp
    seq.pop("placeholder") #remove the first key put into the dict as a placeholder                                                                
    f.close()
    return seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default = None, help="The path to the .jf  database file. Not needed if --contigs is given.")
    parser.add_argument("--reads",nargs='+',default = None, help="The path to the .fasta file(s）containing the contigs to build the jellyfish database. Not needed if --db is provided")
    parser.add_argument("-q","--query", help = "The path to the .fasta query file")
    parser.add_argument("-thre","--threshold", type=int, default = None, help = "The threshold for a bad kmer.")
    parser.add_argument("-rep_thre", type=int, default = None, help = "The threshold  of occurrance for a kmer at a repeitive region.")
    parser.add_argument("-k","--ksize", type=int,help = "The kmer size")
    parser.add_argument("--test", action='store_true',help = "Ouput the indexes of bad kmers, total num of bad kmers, and an estimation for Q value")
    parser.add_argument("--fix", action='store_true', help="Output the index of fixed bases and output the new sequence")
    parser.add_argument("--fout",default = "fout.csv",help = "The path to output the index of the fixed bases." )
    parser.add_argument("-ff","--fixedfasta",default = "fixed_seq.fasta",help = "The path to output the fixed assembly sequences")
    parser.add_argument("--tout", default = "tout.csv", help = "The output file containing the locations of bad kmers")
    args = parser.parse_args()
    main(args.reads,args.query,args.ksize,args.test,args.fix,args.fout,args.tout,args.fixedfasta,args.db,args.threshold,args.rep_thre)

