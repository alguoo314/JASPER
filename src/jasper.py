#!/usr/bin/env python
import sys
import re
import os
import argparse
import math
import csv
import difflib
import dna_jellyfish as jf
#import textwrap

def main(contigs,query_path,k,test,fix,fout,fixedout,database,thre,num_iter):
    try:
        db = database
        if ((contigs == None and database == None) or (contigs != None and database != None)): 
            sys.stderr.write("Wrong arguments. One and only one between the contigs and database argument should be given. ")
            return
        
        threshold,db = jellyfish(contigs, database, k,thre)        

        qf  = jf.QueryMerFile(db)
        for ite in range(num_iter+1): #num_iter rounds of fixing plus one more round to find the final q value
            query_path = iteration(num_iter,ite,qf,query_path,k,test,fix,fout,fixedout,database,threshold)
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \                                                                   
         sys.exit(1)         
            

def iteration(num_iter,ite,qf,query_path,k,test,fix,fout,fixedout,database,threshold):   
    try:
        if ite == num_iter:
            fix=False
        seq_dict = parse_fasta(query_path)
        wrong_kmers_list = []
        seqs = []
        fixed_bases_list = []
        total_wrong_kmers = 0
        total_kmers = 0
        fout = os.path.split(fout)
        fout = fout[0]+"_iter"+str(ite)+"_"+fout[1]
        fixedout = os.path.split(fixedout)
        fixedout = fixedout[0]+"_iter"+str(ite)+"_"+fixedout[1]
        seq_names=[]
        for seqname,seq in seq_dict.items():
            seq_names.append(seqname)
            total_kmers += len(seq)-k+1
            good_before = -1 #index of the last guaranteed good base before the mismatch                                                                                     
            backtracked = False
            i = 0 #first k mer at position 0
            wrong_kmers_list = []                                                                                                       
            while i < len(seq)-k+1:
                mer_string = seq[i:k+i]
                N = mer_string.find('N') #ignore all kmers containing non acgt bases
                if N >= 0:
                    i+=(N+1)
                    continue
                n = mer_string.find('n')
                if n >= 0:
                    i+=(n+1)
                    continue
                match = re.match("^[ACTGactg]*$",mer_string) #other invalid characters other than N or n
                if match is None:
                    i +=1
                    continue
                    
                mer = jf.MerDNA(mer_string).get_canonical()
                occurrance = qf[mer]
                if occurrance < threshold:
                    if backtracked == True:
                        i+=1
                        continue
                    j = i-1
                    backtracked = True
                    occurrance = qf[jf.MerDNA(seq[j:k+j]).get_canonical()] 
                    while occurrance < threshold and j>=0:
                        j = j-1
                        occurrance = qf[jf.MerDNA(seq[j:k+j]).get_canonical()]
                    good_before = j+k-1 #the right base of a good kmer
                    prev_good_count = qf[jf.MerDNA(seq[j:k+j]).get_canonical()]
                    kmer_count = qf[jf.MerDNA(seq[i:k+i]).get_canonical()]                                                            

                    if j == -1: #even the first kmer is bad                                                                                                                               
                        good_before = -1
                    while kmer_count < threshold and i < len(seq)-k+1: #go forward back to i
                        i+=1
                        kmer_count = qf[jf.MerDNA(seq[i:k+i]).get_canonical()]
                        
                    good_after = i #the first base of the first good kmer after the mismatch
                    #A kmer is bad if it is below the threshold AND its count is less than 1/2 of the good k-mers before
                    too_low_flag = False
                    if (qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()] < threshold/2) and (qf[jf.MerDNA(seq[good_before-k+3:good_before+3]).get_canonical()] < threshold/2):
                        too_low_flag = True #Below absolute thres = (relative) threshold/2, then it is an error kmer even if  it's > previous kemr count/2. Flag just for debugging purpose
                    else:
                        while  qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()] >= prev_good_count/2 and good_before-k+1 < good_after: #kmer count > prev kmer count/2, then not a bad kmer.
                            if good_before == -1:
                                break
                            if (prev_good_count >= threshold/2 and (qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()] < threshold/2) and (qf[jf.MerDNA(seq[good_before-k+3:good_before+3]).get_canonical()] < threshold/2)):
                                too_low_flag = True
                                break
                            #move to the next kmer
                            prev_good_count = qf[jf.MerDNA(seq[good_before-k+2:good_before+2]).get_canonical()]
                            good_before +=1                        
                        if good_before >= len(seq)-1:
                            break
                    #additional check for case 000...high high...000
                    second = seq[max(0,good_before-k+2)+1:max(0,good_before-k+2)+k+1]
                    k_minus_1 = seq[max(0,good_before-k+2)+k-2:max(0,good_before-k+2)+k+k-2]
                    k_th =seq[max(0,good_before-k+2)+k-1:max(0,good_before-k+2)+k+k-1]
                    k_plus_1 = seq[max(0,good_before-k+2)+k:max(0,good_before-k+2)+k+k]
                    if qf[jf.MerDNA(second).get_canonical()] < threshold and qf[jf.MerDNA(k_minus_1).get_canonical()] < threshold and qf[jf.MerDNA(k_th).get_canonical()] < threshold and qf[jf.MerDNA(k_plus_1).get_canonical()]>=threshold:
                        good_after = max(0,good_before-k+2)+k #ie good_before+2
                    to_be_fixed = seq[max(0,good_before-k+2):good_after+k-1]
                    wrong_kmers_list.extend([*range(max(0,good_before-k+2),good_after)])
                    
                    if fix == True:
                        if good_before < 0:
                            continue
                        seq,fixed_base,original,fixed_ind = fixing_sid(seq,to_be_fixed,k,threshold,qf,len([*range(max(0,good_before-k+2),good_after)]),good_before,good_after) #fix simple sub/insert/del cases
                        if fixed_base != "nN":
                            if len(fixed_ind) == 1:
                                fixed_bases_list.append([seqname,fixed_ind[0],fixed_base,original])
                            else:
                                fixed_bases_list.append([seqname,fixed_ind[0],fixed_base[0],original[0]])
                                fixed_bases_list.append([seqname,fixed_ind[1],fixed_base[1],original[1]])
                    
                else: #good kmer                                                                                                                                                      
                    backtracked = False
                    i += k-1
            
                    
            seqs.append(seq)
            total_wrong_kmers += len(wrong_kmers_list)


        if test == True:
            if ite == 0 or ite==num_iter:
                file_name = str(ite)+"qValCalcHelper.csv"
                with open(file_name,'a') as f:
                    f.write("{} {}\n".format(total_wrong_kmers,total_kmers))
                        
            
        
        if fix == True:
            base_fields = ['Contig', 'Base_coord', 'Original',"Mutated"]
            with open(fout,'w') as csvf:
                csvwriter = csv.writer(csvf,delimiter=' ')
                csvwriter.writerow(base_fields)
                csvwriter.writerows(fixed_bases_list)

            i = 0
            with open(fixedout,'w') as of:
                for seqname in seq_names:
                    of.write(">{}\n".format(seqname))
                    new_seq = split_output(seqs[i],60)
                    for l in new_seq:
                        of.write(l+"\n")
                    i+=1
            return fixedout
        
        return 
    except:
        exception_type, exception_object, exception_traceback = sys.exc_info()
        line_number = exception_traceback.tb_lineno
        print(line_number)
        print(sys.exc_info()) #to help debug
        print(jf.MerDNA(seq[j:k+j]).get_canonical())
        print(seq[j:k+j])
        sys.exit(1)




def split_output(seq, num_per_line=60): #make a new line after num_per_line bases
    lines = math.ceil(len(seq)/num_per_line)
    output = []
    for i in range(lines):
        output.append(seq[num_per_line*i:num_per_line*(i+1)])
    return output

         

def jellyfish(contigs,database,k,thre):
    count = math.inf
    threshold = 0
    if database != None:
        db_name = database
        base_name = os.path.splitext(db_name)[0]
    if contigs != None:
        db_name = os.path.splitext(os.path.basename(contigs[0]))[0]+".jf"
        base_name = os.path.splitext(db_name)[0]
        contigs = ' '.join(contig_file for contig_file in contigs)
        os.system("jellyfish count -s 300000000 -t 32 -m {} -C -o {} {}".format(k,db_name,contigs))
    if thre != None:
        return thre,db_name
    else:
        os.system("jellyfish histo -t 32 {}> {}".format(db_name,base_name+".csv"))
        with open(base_name+".csv",'r') as histo:
            csvreader = csv.reader(histo,delimiter=' ')
            for row in csvreader:
                if count >= int(row[-1]):
                    count = int(row[-1])
                    threshold = int(int(row[0])/2)
                else: #found local min
                    return threshold,db_name
                



def fixing_sid(seq,to_be_fixed,k,threshold,qf,num_below_thres_kmers,good_before,good_after): #fix sub and indel
    try:
        fixed_base = "nN"
        original = '-'
        fixed_ind = None

        if num_below_thres_kmers == k: #substitution or insertion
            b = fix_k_case_sub(to_be_fixed,k,threshold,qf)
            if b !=  None:
                original = "s"+seq[good_after-1]
                seq  = seq[:good_after-1] + b + seq[good_after:]
                fixed_base = b
                fixed_ind = [good_after-1]
            else:
                b = fix_insert(to_be_fixed,k,threshold,qf)
                if b != None:
                    original = "i"+seq[good_after-1]
                    seq = seq[:good_after-1] + seq[good_after:]
                    fixed_base = '-'
                    fixed_ind = [good_after-2]
                                    
        elif num_below_thres_kmers == k-1: #deletion or same base insertion or substitution
           removed_base = fix_del(to_be_fixed,k,threshold,qf)
           if removed_base != None: #deletion of a base
               original = "d-"
               fixed_ind = [good_after]
               seq = seq[:good_after]+removed_base+seq[good_after:]
               fixed_base = removed_base 
           elif fix_same_base_insertion(to_be_fixed,k,threshold,qf,num_below_thres_kmers)!= None:
               original = "i"+seq[good_before]
               seq = seq[:good_before] + seq[good_before+1:]
               fixed_base = "-"
               fixed_ind = [good_before] #insertion after this index
           else:
                left,right,l_or_r =  fixdiploid(to_be_fixed,k,threshold,qf,seq,good_before,good_after) #diploidy of two adjacent bases                                
                if l_or_r !=  None:
                    if l_or_r == "s": #lefting base is changed                                                                                                        
                        original = "s"+seq[good_after-1]
                        fixed_base = str(left)
                        fixed_ind = [good_after-1]
                    else:#the righting base is changed                                                                                                                
                        original = "s"+seq[good_before+1]
                        fixed_base = str(right)
                        fixed_ind = [good_before+1]
                    seq = seq[:good_after-1]+left+seq[good_after:good_before+1] + right + seq[good_before+2:]
                    

                  
        elif num_below_thres_kmers < k-1 and num_below_thres_kmers > 1 and len(to_be_fixed)>=k:#skip the good_before = -1 (ie first kmer is bad) case.
            left,right,l_or_r =  fixdiploid(to_be_fixed,k,threshold,qf,seq,good_before,good_after) #diploidy
            if l_or_r !=  None: 
                if l_or_r == "s": #lefting base is changed
                    original = "s"+seq[good_after-1]
                    fixed_base = str(left)
                    fixed_ind = [good_after-1]
                else:#the righting base is changed
                    original = "s"+seq[good_before+1]
                    fixed_base = str(right)
                    fixed_ind = [good_before+1]
                seq = seq[:good_after-1]+left+seq[good_after:good_before+1] + right + seq[good_before+2:]
                
            else: 
                inserted_base = fix_same_base_insertion(to_be_fixed,k,threshold,qf,num_below_thres_kmers)
                if inserted_base != None:
                    original = "i"+inserted_base
                    seq = seq[:good_before] + seq[good_before+len(inserted_base):]     
                    fixed_base = "-"
                    fixed_ind = [good_before] #insertion after this index
                else:
                    removed_base = fix_same_base_del(to_be_fixed,k,threshold,qf,num_below_thres_kmers)
                    if removed_base != None: #deletion of a base
                        original = "d-"
                        fixed_ind = [good_before]
                        seq = seq[:good_before]+removed_base+seq[good_before:]
                        fixed_base = removed_base

        elif num_below_thres_kmers > k: #two or more nearby errors.
            good_kmer_before = seq[good_before-k+1:good_before+1] 
            good_k_mer_after = seq[good_after:good_after+k] 
            fixed_seq = base_extension(to_be_fixed,qf,k,good_kmer_before,good_k_mer_after,threshold)
            if fixed_seq != None:
                seq = seq[:good_before+1]+fixed_seq+seq[good_after:]
                fixed_ind = []
                fixed_base = []
                original = []
                for index, s in enumerate(difflib.ndiff(fixed_seq,to_be_fixed)):
                    if s[0]=='-': #fixed a deletion
                        original.append("d-")
                        fixed_ind.append(index+(good_before+1))
                        fixed_base.append(s[-1])
                    elif s[0]=='+': #fixed an insertion
                        fixed_base.append('-')
                        original.append("i"+s[-1])
                        fixed_ind.append(index+good_before+1)
                        
                        
        return seq,fixed_base, original, fixed_ind
    except:
        exception_type, exception_object, exception_traceback = sys.exc_info()
        line_number = exception_traceback.tb_lineno
        print(line_number)
        print(sys.exc_info()) #to help debug                                  \                                                                   
        sys.exit(1)         
                
def fixdiploid(seq_to_be_fixed,k,threshold,qf,full_seq,good_before,good_after):
    #for fixing the first base of the last k-mer in L we also need to check the kmers before that contains this base.
    # for fixing the last b of first kmer in L we also need to check if the kmer starting at good_after is good still
    try:
        left_bad = seq_to_be_fixed[len(seq_to_be_fixed)-k]
        right_bad = seq_to_be_fixed[k-1]
        left = left_bad
        right = right_bad
        good_before_starting_ind = max(0,good_before-k+1)
        base_after=''
        if good_after+k-1+int((k-1-len(seq_to_be_fixed)+k)/2) < len(full_seq):
            base_after = full_seq[good_after+k-1:good_after+k-1+int((k-1-len(seq_to_be_fixed)+k)/2)] #the end of the kmer starting with the right changed base = right index-left index -1 + good_after+k-1 
        else:
            base_after = full_seq[min(len(full_seq)-1,good_after+k-1):len(full_seq)]
        before_len=len(base_after)
        bases_before=full_seq[max(0,good_before_starting_ind-before_len+1):good_before_starting_ind+1]
        for x in 'ACTG':
            for y in 'ACTG':
                if x ==left_bad and y==right_bad:
                    continue
                if x!=left_bad and y!=right_bad:
                    continue
                trial = seq_to_be_fixed[:len(seq_to_be_fixed)-k]+x+seq_to_be_fixed[len(seq_to_be_fixed)-k+1:k-1]+y+seq_to_be_fixed[k:]
                fixed = True
                check=bases_before+trial+base_after
                for i in  range(len(check)-k+1):
                    if qf[jf.MerDNA(check[i:k+i]).get_canonical()] < threshold:
                        fixed  = False
                        break
                if fixed == True:
                    left = x
                    right = y
                    if x == left_bad:
                        l_or_r = "e" 
                        
                    elif y == right_bad:
                        l_or_r = "s"
                        
                    else: #both changed not acceptable
                        continue
                        
     
                    return(left,right,l_or_r)
        return left,right,None
    except:
         exception_type, exception_object, exception_traceback = sys.exc_info()
         line_number = exception_traceback.tb_lineno
         print(line_number)
         print(sys.exc_info()) #to help debug                                  \                                                                   
         sys.exit(1)    
    


def fix_k_case_sub(seq_to_be_fixed,k,threshold,qf): #when number of conseuctive bad kmers is k
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
                    break
            if fixed == True:
                return b
    return None


    
def fix_nearby_subs(seq_to_be_fixed,k,threshold,qf,num_below_thres_kmers): #a substituion nearby another error
    for x in 'ACTG':
        if x == seq_to_be_fixed[k-1]:
            continue
        for y in 'ACTG':
            fixed = False
            if y == seq_to_be_fixed[len(seq_to_be_fixed)-k]:
                continue
            else:
                trial = seq_to_be_fixed[:k-1]+x+seq_to_be_fixed[k:len(seq_to_be_fixed)-k]+y+seq_to_be_fixed[len(seq_to_be_fixed)-k+1:]
                fixed = True
                for i in range(len(trial)-k+1):
                    if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
                        fixed  = False
                        break
                if fixed == True:
                    return x,y
                else:
                    l = len(seq_to_be_fixed)
                    if qf[jf.MerDNA(seq_to_be_fixed[:k-1]+x).get_canonical()]  >= threshold and qf[jf.MerDNA(seq_to_be_fixed[1:k-1]+x+seq_to_be_fixed[k]).get_canonical()] and qf[jf.MerDNA(seq_to_be_fixed[2:k-1]+x+seq_to_be_fixed[k:k+2]).get_canonical()] >= threshold: # >3 fixed
                        return x,None
                    elif qf[jf.MerDNA(seq_to_be_fixed[l-k-2:l-k]+y+seq_to_be_fixed[l-k+1:l-2]).get_canonical()] >= threshold and  qf[jf.MerDNA(seq_to_be_fixed[l-k-1]+y+seq_to_be_fixed[l-k+1:l-1]).get_canonical()]  >= threshold and qf[jf.MerDNA(y+seq_to_be_fixed[l-k+1:]).get_canonical()] >= threshold:
                       return None,y
    return None, None
               
                

def fix_insert(seq_to_be_fixed,k,threshold,qf):
    ind_to_be_removed = k-1
    base_to_be_removed = seq_to_be_fixed[ind_to_be_removed]
    seq_to_be_fixed = seq_to_be_fixed[:ind_to_be_removed] + seq_to_be_fixed[ind_to_be_removed+1:]
    fixed = True
    for i in range(len(seq_to_be_fixed)-k+1):
        if qf[jf.MerDNA(seq_to_be_fixed[i:k+i]).get_canonical()] < threshold:
            fixed  = False
            break
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
        if fixed == True:
            return alt
    return None

def fix_same_base_del(seq_to_be_fixed,k,threshold,qf,num_below_thres_kmers):
    sb = seq_to_be_fixed[k-2] #sb stands for same base
    fixed = False
    trial=seq_to_be_fixed
    new_bad = 0
    inserted = 0
    original_bad = len(seq_to_be_fixed)-k+1
    current_bad = original_bad
    max_insertions = original_bad
    while inserted < max_insertions:
        new_bad=0
        trial = trial[:k-1]+sb+trial[k-1:]
        fixed = True
        inserted+=1
        for i in range(len(trial)-k+1):
            if qf[jf.MerDNA(trial[i:k+i]).get_canonical()] < threshold:
                fixed  = False
                new_bad +=1
        if fixed == True:
            return sb*inserted
        if (new_bad >= current_bad):
            return None
        else: #added one base may have helped but need more
            current_bad = new_bad
            continue
    return None


def fix_same_base_insertion(seq_to_be_fixed,k,threshold,qf,num_below_thres_kmers):
    ind_to_be_removed = k-1
    sb = seq_to_be_fixed[k-1] #sb stands for same base
    fixed = False
    deleted = 0
    original_bad = len(seq_to_be_fixed)-k+1
    current_bad = original_bad
    max_deletions = original_bad
    while seq_to_be_fixed[k-1] == sb and deleted < max_deletions:
        deleted +=1
        seq_to_be_fixed = seq_to_be_fixed[:ind_to_be_removed] + seq_to_be_fixed[ind_to_be_removed+1:]
        fixed=True
        new_bad = 0
        flag = 0
        for i in  range(len(seq_to_be_fixed)-k+1):
            flag = 1
            if qf[jf.MerDNA(seq_to_be_fixed[i:k+i]).get_canonical()] < threshold:
                fixed=False
                new_bad +=1
        if (flag == 0):
            return None
        if (fixed == True):
            return sb*deleted
        if (new_bad >= current_bad):
            return None
        else: #delete one more base
            current_bad = new_bad
            continue
    return None


def base_extension(seq_to_be_fixed,qf,k,good_kmer_before,good_k_mer_after,threshold):
    if len(good_kmer_before) < k:
        raise Exception("Path not long enough") #remove this later                                                                                                                            
    bases = ["A", "C", "G", "T"]
    paths = []  # array of all possible extensions                                                                                                                                            
    #K+L-1-2(K-1) = L-K+1 deleted                                                                                                                                                             
    L = len(seq_to_be_fixed)+1-k #number of bad kmers                                                                                                                                         
    max_ext = L-k+5 #allow insertion of 4 more bases than before                                                                                                                              
    found_good_path = False
    start_km1 = good_kmer_before[:-1]   # store the k-1 bases in a variable no need to carry these around                                                                                     
    paths.append(good_kmer_before[-1])  # the last base of the initial k-mer makes the first path                                                                                             
    right_end = good_k_mer_after[:-1]
    for i in range(1,1+max_ext):
        paths = [l for l in paths if len(l) > 0]
        if len(paths) > 2000: #may change later
            return None
        last_path = len(paths) - 1  # since the number of paths will be changing we need to record this                                                                                       
        for p in range(last_path+1):
            if not paths[p]:
                continue
            km1 = (start_km1 + paths[p])[-k+1:]  # get the k-1 bases off the end                                                                                                              
            path_ext_count = 0  # this becomes 1 if we find an extension, and if no extension, then we delete the path                                                                        
            for j in range(4):  # try to extend                                                       
                score = qf[jf.MerDNA(km1 + bases[j]).get_canonical()]
                if score >= threshold:
                    if i >= max_ext-8: #the path is <4 bases short of deleted part. Start checking if it is finished.                                                                         
                        found_good_path = True
                        right_half_connected = km1+ bases[j]+right_end
                        for n in range(1,len(right_half_connected)-k+1): #check the remaining kmers               
                            if qf[jf.MerDNA(right_half_connected[n:k+n]).get_canonical()] < threshold:
                                found_good_path = False
                                break
                        if found_good_path:
                            return paths[p][1:]

                    #not a good path YET                                                                                                                                                      
                    if path_ext_count == 0:  # first extension -- extend the current path                                                                                                     
                        paths[p] += bases[j]
                        path_ext_count = 1
                    else:  # another extension -- add a new path                                                                                                                              
                        paths.append(paths[p][:-1] + bases[j])
            if path_ext_count == 0:
                paths[p] = ""  # no extensions, kill this path 
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
    parser.add_argument("-k","--ksize", type=int,help = "The kmer size")
    parser.add_argument("--test", action='store_true',help = "Ouput the total num of bad kmers, and an estimation for Q value")
    parser.add_argument("--fix", action='store_true', help="Output the index of fixed bases and output the new sequence")
    parser.add_argument("--fout",default = "fout.csv",help = "The path to output the index of the fixed bases." )
    parser.add_argument("-ff","--fixedfasta",default = "fixed_seq.fasta",help = "The path to output the fixed assembly sequences")
    parser.add_argument("-p","--num_passes", type=int, default = 2, help = "The number of iterations of fixing.")
    args = parser.parse_args()
    main(args.reads,args.query,args.ksize,args.test,args.fix,args.fout,args.fixedfasta,args.db,args.threshold,args.num_passes)
