## Script to analyse and identify IGHV sequences ##

# Alaine Athenaïs Marsden - 2024: alainem@nicd.ac.za #

################################################################################################################

## Require Pre-Processing Steps ##


#= This script requires several inputs:

1. preprocessed participant sequences - qc'ed, trimmed, dereplicated
2. reference IGHV region sequences: incl -> v-region, RSS, LPART1 and LPART2
3. Tables generated from local BLAST with custom v-region db

script assumes a directory full of participant specific directories, aptly named.

=#

################################################################################################################

using FASTX, StatsBase, OrderedCollections, BioSequences, JLD2, BioAlignments, CSV, DataFrames

################################################################################################################

# load in data #

dir = "../../Data/pacbio_june_2022"
data = readdir(dir)
data = filter(x -> startswith(x, "CAP"), data)

seqs = []

for i in data
    cap_file = dir*"/"*i*"/constant_genes.fa"
    r = FASTA.Reader(open(cap_file, "r"))
    for x in r
        if occursin("IGHA", identifier(x))
            push!(seqs,x)
        end
    end
end

################################################################################################################

# correct orientation #

forward_motif = ExactSearchQuery(dna"GGGCCGCGTCCTCACAGTGCATTCTGTGTTCCAGCATCCCCGACCAGCCCC")
reverse_motif = ExactSearchQuery(dna"CAGCCCCACGCTTCCATCCGGCGCCTGTCTG")

corrected_seqs = []

for i in seqs 
    if occursin(forward_motif, sequence(i))
        if occursin(reverse_motif, sequence(i))
            first_pos = first(findfirst(forward_motif, sequence(i)))
            last_pos = last(findfirst(reverse_motif, sequence(i)))
            new_seq = FASTA.Record(identifier(i), sequence(i)[first_pos:last_pos])
            push!(corrected_seqs, new_seq)
        end
    end
end

################################################################################################################

# further processing #

cap_list = []
compiled_seqs = []
processed_seqs = []
n = 1

for i in corrected_seqs
    if identifier(i) ∈ processed_seqs
        continue
    end

    if isempty(sequence(i))
        continue
    end

    parts = split(identifier(i), "_")
    caps = [parts[1]]
    size = parse(Int64, parts[3])
    query_seq = ExactSearchQuery(sequence(i))
    gene = split(parts[4], "*")[1]
    for x in corrected_seqs
        if isempty(sequence(x))
            continue
        end
        if occursin(query_seq, sequence(x))
            if identifier(x) == identifier(i)
                continue
            elseif identifier(x) ∈ processed_seqs
                continue
            else
                parts_x = split(identifier(x), "_")
                if parts_x ∉ caps
                    push!(caps, parts_x[1])
                end
                size = size + parse(Int64, parts_x[3])
                push!(processed_seqs, identifier(x))
            end
        end
    end
    push!(cap_list, caps)
    new_id = join([gene, string(n), string(size), string(length(caps))], "_")
    new_seq = FASTA.Record(new_id, sequence(i))
    push!(compiled_seqs, new_seq)
    n = n + 1
    push!(processed_seqs, identifier(i))
end


filtered_seqs = []
filtered_caps = []
n = 1

for i in compiled_seqs
    parts = split(identifier(i), "_")
    size = parse(Int64, parts[3])
    count = parse(Int64, parts[4])
    if count > 1
        push!(filtered_seqs, i)
        push!(filtered_caps, cap_list[n])
    elseif size > 4
        push!(filtered_seqs, i)
        push!(filtered_caps, cap_list[n])
    else
    end
    n = n + 1
end

################################################################################################################

# Identifying sequences and extracting CH regions #

w=open(FASTA.Writer, "../../Data/Constant_heavy/filtered_igha.fa")

for i in compiled_seqs
    write(w, i)
end

close(w)

id_list = []

for i in compiled_seqs
    push!(id_list, identifier(i))
end

ref_dict = Dict(zip(id_list, filtered_caps))
save("igha_seq_ref.jld2", "data", ref_dict)

r = FASTA.Reader(open("/Users/alaine/Refs/Constant/IGHA_CH.fa", "r"))

ch1 = []
ch2 = []
ch3 = []

for i in r
    gene = split(identifier(i), "|")[2]
    part = split(description(i), "|")[3]
    new_id = gene*"_"*part
    new_seq = FASTA.Record(new_id, sequence(i))
    if occursin("CH1", identifier(new_seq))
        push!(ch1, new_seq)
    elseif occursin("CH2", identifier(new_seq))
        push!(ch2, new_seq)
    else
        push!(ch3, new_seq)
    end
end

closest_match_ch1 = []
closest_match_ch2 = []
closest_match_ch3 = []
mm_ch1 = []
mm_ch2 = []
mm_ch3 = []
ch1_motif= []
ch2_motif= []
ch3_motif = []

scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-10, gap_extend=-1)

for i in compiled_seqs

    seq_query = sequence(i)
    mm_test_ch1 = 50
    mm_test_ch2 = 50
    mm_test_ch3 = 50
    ch1_seq = nothing
    ch2_seq = nothing
    ch3_seq = nothing
    match_ch1 = nothing
    match_ch2 = nothing
    match_ch3 = nothing

    for x in ch1 
        seq_ref = sequence(x)
        aln = alignment(pairalign(GlobalAlignment(), seq_ref, seq_query, scoremodel))
        gapped = LongDNA{4}([x for (x, _) in aln])
        if count_mismatches(aln) < mm_test_ch1   
            trim_f = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_f = i
                    break
                else
                    continue
                end
            end
            reverse!(gapped)
            trim_r = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_r = i
                    break
                else
                    continue
                end
            end

            ch1_seq = seq_query[trim_f:(length(seq_query)-trim_r+1)]
            mm_test_ch1 = count_mismatches(aln)
            match_ch1 = split(identifier(x), "_")[1]
        end
    end

    push!(mm_ch1, mm_test_ch1)
    push!(ch1_motif, ch1_seq)
    push!(closest_match_ch1, match_ch1)

    for x in ch2 
        seq_ref = sequence(x)
        aln = alignment(pairalign(GlobalAlignment(), seq_ref, seq_query, scoremodel))
        gapped = LongDNA{4}([x for (x, _) in aln])
        if count_mismatches(aln) < mm_test_ch2   
            trim_f = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_f = i
                    break
                else
                    continue
                end
            end
            reverse!(gapped)
            trim_r = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_r = i
                    break
                else
                    continue
                end
            end

            ch2_seq = seq_query[trim_f:(length(seq_query)-trim_r+1)]
            mm_test_ch2 = count_mismatches(aln)
            match_ch2 = split(identifier(x), "_")[1]
        end
    end

    push!(mm_ch2, mm_test_ch2)
    push!(ch2_motif, ch2_seq)
    push!(closest_match_ch2, match_ch2)

    for x in ch3 
        seq_ref = sequence(x)
        aln = alignment(pairalign(GlobalAlignment(), seq_ref, seq_query, scoremodel))
        gapped = LongDNA{4}([x for (x, _) in aln])
        if count_mismatches(aln) < mm_test_ch3  
            trim_f = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_f = i
                    break
                else
                    continue
                end
            end
            reverse!(gapped)
            trim_r = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_r = i
                    break
                else
                    continue
                end
            end

            ch3_seq = seq_query[trim_f:(length(seq_query)-trim_r+1)]
            mm_test_ch3 = count_mismatches(aln)
            match_ch3 = split(identifier(x), "_")[1]
        end
    end

    push!(mm_ch3, mm_test_ch3)
    push!(ch3_motif, ch3_seq)
    push!(closest_match_ch3, match_ch3)

end
        


# read in table

sum_tab = DataFrame(CSV.File("/Users/alaine/Data/Constant_heavy/IGHA_features.csv"))

tab = filter(:reads => x -> x > 9, sum_tab)

new_cap_list = []

for i in tab.cap_list
    ls = split(i, "\"")
    new_list = []
    for x in ls
        if startswith(x, "CAP")
            if x ∉ new_list
                push!(new_list, x)
            end
        end
    end
    new_list = join(new_list, ", ")
    push!(new_cap_list, new_list)
end

tab[!, "cap_list"] = new_cap_list

#read in hinge Sequences

r = FASTA.Reader(open("/Users/alaine/Refs/Constant/igha_hinge.fa", "r"))

hinges = []
for i in r
    push!(hinges, i)
end

scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-20, gap_extend=-5)
closest_match_h = []
mm_h = []
h_motif= []

for i in tab.seq
    seq_query = LongDNA{4}(i)
    mm_test_h = 20
    h_seq = nothing
    match_h = nothing
    for x in hinges 
        seq_ref = sequence(x)
        aln = alignment(pairalign(LocalAlignment(), seq_ref, seq_query, scoremodel))
        gapped = LongDNA{4}([x for (x, _) in aln])
        if count_mismatches(aln) < mm_test_h  
            trim_f = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_f = i
                    break
                else
                    continue
                end
            end
            reverse!(gapped)
            trim_r = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_r = i
                    break
                else
                    continue
                end
            end

            h_seq = seq_query[trim_f:(length(seq_query)-trim_r+1)]
            mm_test_h = count_mismatches(aln)
            match_h = split(identifier(x), "_")[1]
        end
    end

    push!(mm_h, mm_test_h)
    push!(h_motif, h_seq)
    push!(closest_match_h, match_h)
    
end

tab[!, "h_seq"] = h_motif

closest_match_h3 = []
mm_h3 = []
h3_motif= []

for i in tab.seq
    seq_query = LongDNA{4}(i)
    mm_test_h3 = 50
    h3_seq = nothing
    match_h3 = nothing
    for x in ch3
        seq_ref = sequence(x)
        aln = alignment(pairalign(SemiGlobalAlignment(), seq_ref, seq_query, scoremodel))
        gapped = LongDNA{4}([x for (x, _) in aln])
        if count_mismatches(aln) < mm_test_h3  
            trim_f = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_f = i
                    break
                else
                    continue
                end
            end
            reverse!(gapped)
            trim_r = 0
            for i in 1:length(gapped)
                if gapped[i] != DNA_Gap
                    trim_r = i
                    break
                else
                    continue
                end
            end

            h3_seq = seq_query[trim_f:(length(seq_query)-trim_r+1)]
            mm_test_h3 = count_mismatches(aln)
            match_h3 = split(identifier(x), "_")[1]
        end
    end

    push!(mm_h3, mm_test_h3)
    push!(h3_motif, h3_seq)
    push!(closest_match_h3, match_h3)
    
end

tab[!, "ch3_seq"] = h3_motif

CSV.write("/Users/alaine/Data/Constant_heavy/igha_tab.csv", tab)  

spliced_seqs = []

for i in eachrow(tab)
    seq = i.ch1_seq*string(i.h_seq)*i.ch2_seq*string(i.ch3_seq)
    rec = FASTA.Record(i.seq_id,  seq)
    push!(spliced_seqs, rec)
end



x = open(FASTA.Writer, "/Users/alaine/Data/Constant_heavy/spliced_IGHA.fa")

for i in spliced_seqs
    write(x, i)
end

close(x)
                
#amino acid translation

peptides = []
aa_ids = []
for i in spliced_seqs
    transseq = sequence(i)[3:end]
    if isinteger(length(transseq)/3)
        push!(peptides, translate(transseq))
        push!(aa_ids, identifier(i))
    end
end

filter!(:seq_id => n -> n ∈ aa_ids, tab)
tab[!, "spliced_aa_seq"] = peptides



filtered_sseqs = []

for i in spliced_seqs
    if identifier(i) ∈ aa_ids
        push!(filtered_sseqs, i)
    end
end

tab[!, "spliced_nt_seq"] = sequence.(filtered_sseqs)


#checking for AA changes

spliced_refs = []
r = FASTA.Reader(open("/Users/alaine/Refs/Constant/spliced_refs.fa", "r"))
scoremodel = AffineGapScoreModel(PAM30, gap_open=-20, gap_extend=-1)

for i in r
    new_id = split(identifier(i), "|")[2]
    new = FASTA.Record(new_id,translate(sequence(i)[3:end]))
    push!(spliced_refs, new)
end

aa_array = []
aa_mm = []

for row in eachrow(tab)
    if row.novel_call == "novel"
        seq1 = row.spliced_aa_seq
        mm = 10
        aa_changes = []
        gene = split(row.seq_id, "_")[1]
        for x in spliced_refs
            if startswith(identifier(x), gene)
                seq2 = sequence(x)
                aln = alignment(pairalign(SemiGlobalAlignment(), seq1, seq2, scoremodel))
                if count_mismatches(aln) < mm
                    mm = count_mismatches(aln)
                    n = 1
                    del_present = false
                    first_del = nothing
                    last_aa = nothing
                    for y in collect(aln)
                        if del_present
                            if y[1] != AA_Gap
                                mut = first_del*"_"*last_aa*"del"
                                push!(aa_changes, mut)
                                del_present = false
                            end
                        elseif y[1] != y[2]
                            if y[1] == AA_Gap
                                del_present = true
                                first_del = string(y[2])*string(n)
                                last_aa = string(y[2])*string(n)
                                n = n + 1
                                continue
                            end
                            mut = string(y[2])*string(n)*string(y[1])
                            push!(aa_changes, mut)
                        end
                        last_aa = string(y[2])*string(n)
                        n = n + 1
                    end
                end
            end
        end

        if isempty(aa_changes)
            push!(aa_array, nothing)
            push!(aa_mm, 0)
        else
            push!(aa_array, join(aa_changes, ", "))
            push!(aa_mm, mm)
        end


    else
        push!(aa_array, nothing)
        push!(aa_mm, 0)
    end
end

tab[!, "aa_changes"] = aa_array
tab[!, "aa_mismatches"] = aa_mm

CSV.write("/Users/alaine/Data/Constant_heavy/igha_tab.csv", tab, transform=(col, val) -> something(val, missing)) 



x = open(FASTA.Writer, "/Users/alaine/Data/Constant_heavy/spliced_IGHA_nt.fa")

for i in filtered_sseqs
    write(x, i)
end

close(x)

################################################################################################################

# constructing participant genotypes #


# genotypes

caps = []

for i in tab.cap_list
    a = split(i, ", ")
    for x in a
        if x ∉ caps
            push!(caps, x)
        end
    end
end

geno_tab = DataFrame()

for i in caps
    cap_tab = filter(:cap_list => n -> occursin(i, n), tab)
    sort!(cap_tab, 4, rev=true)
    if nrow(cap_tab) > 1
        cap_tab = cap_tab[1:2,:]
    end
    select!(cap_tab, Not(:cap_list))
    new_id = []
    for x in eachrow(cap_tab)
        push!(new_id, i)
    end
    cap_tab[!, "cap_id"] = new_id
    geno_tab = vcat(geno_tab, cap_tab)

end

CSV.write("/Users/alaine/Data/Constant_heavy/geno_tab.csv", geno_tab, transform=(col, val) -> something(val, missing)) 

## getting common phenotypic alleles
# using genotyped Sequences

geno_tab = DataFrame(CSV.File("/Users/alaine/Data/Constant_heavy/geno_tab.csv"))

seqs = unique(geno_tab.seq_id)
select!(geno_tab, Not(:cap_id))
unique!(geno_tab)

novel_count = 0

for i in geno_tab.novel_call
    if i == "novel"
        novel_count = novel_count + 1
    end
end

aa_change_count = 0

for i in geno_tab.aa_changes
    if !ismissing(i)
        aa_change_count = aa_change_count + 1
    end
end
replace!(geno_tab.aa_changes, missing => "none")

aa_types = unique(geno_tab.aa_changes)


counts_aa = []
participants_aa = []
seqs_aa = []
seq_lists = []
seq_count_aa = []

for i in aa_types
    s_list = []
    s_ids = []
    counts = 0
    partis = 0
    for row in eachrow(geno_tab)
        if row.aa_changes == i
            partis = partis + row.participants
            counts = counts + row.reads
            push!(s_ids, row.seq_id)
            push!(s_list, row.spliced_nt_seq)
        end
    end
    push!(seq_count_aa, length(s_ids))
    push!(counts_aa, counts)
    push!(participants_aa, partis)
    push!(seq_lists, s_ids)
    push!(seqs_aa, s_list)
end

aa_type_tab = DataFrame()
aa_type_tab[!, "aa_type"] = aa_types
aa_type_tab[!, "spliced_nt_seq"] = seqs_aa
aa_type_tab[!, "gl_seqids"] = seq_lists
aa_type_tab[!, "participants"] = participants_aa
aa_type_tab[!, "reads"] = counts_aa
aa_type_tab[!, "sequence_counts"] = seq_count_aa



CSV.write("/Users/alaine/Data/Constant_heavy/aa_alleles_tab.csv", aa_type_tab) 

aa_type_tab = DataFrame(CSV.File("/Users/alaine/Data/Constant_heavy/aa_alleles_tab.csv"))

bur_numbered = []

for i in aa_type_tab.aa_type
    if i == "none"
        push!(bur_numbered, "none")
    else
        bur_list = []
        for x in split(i, ", ")
            bur = first(x)*string(parse(Int64, chop(x, head = 1, tail = 1)) + 121)*last(x)
            push!(bur_list, bur)
        end
        push!(bur_numbered, join(bur_list, ", "))
    end
end

aa_type_tab[!, "bur_numbered"] = bur_numbered

capids = []
            
for i in aa_type_tab.gl_seqids
    cs = []
    for x in i
        id_tab = filter(:seq_id => n -> occursin(x, n), geno_tab)
        c_list = unique(id_tab.cap_id)
        for c in c_list
            if c ∉ cs
                push!(cs, c)
            end
        end
    end
    push!(capids, cs)
end

aa_type_tab[!, "participants"] =  length.(capids)
aa_type_tab[!, "cap_ids"] = capids

