
## Script to analyse and identify IGHG sequences ##

# A. Marsden - 2024 #

################################################################################################################

using FASTX, StatsBase, OrderedCollections, BioSequences, JLD2, BioAlignments, CSV, DataFrames

################################################################################################################

# function to get cap IGHG sequences from split gene fastas
function compile_cap_seqs(file)
    seqs = FASTA.Reader(open(file, "r"))
    seq_list = []
    for i in seqs
        if startswith(identifier(i), "CAP")
            mm = parse(Int64, split(identifier(i), "_")[4])
            if mm < 10
                push!(seq_list, i)
            end
        end
    end
    return seq_list
end
# general function to compile sequences from fastas
function compile_seqs(file)
    seqs = FASTA.Reader(open(file, "r"))
    seq_list = []
    for i in seqs
        push!(seq_list, i)
    end
    return seq_list
end
# function to reduce number of sequences per cap participant to 2, using seq_size as the deciding variable
function two_seq_cap(seq_list)
    cap_list = []
    for i in seq_list
        cap = split(identifier(i), "_")[1]
        if cap ∉ cap_list
            push!(cap_list,cap)
        end
    end

    new_seq_list = []
    for cap in cap_list
        cap_seqs = []
        for i in seq_list
            if startswith(identifier(i), cap)
                push!(cap_seqs, i)
            end
        end
        p1 = true
        p2 = true
        for i in cap_seqs
            size = parse(Int64, split(identifier(i), "_")[3])
            if p1 == true
                p1 = i
            elseif p2 == true
                p2 = i
            else
                size_p1 = parse(Int64, split(identifier(p1), "_")[3])
                size_p2 = parse(Int64, split(identifier(p2), "_")[3])
                if size > size_p1
                    p1 = i
                elseif size > size_p2
                    p2 = i
                end
            end
        end
        push!(new_seq_list, p1)
        if p2 == true
            continue
        end
        push!(new_seq_list, p2)
    end
    return new_seq_list
end
#define array for ch designation
struct ch_array
    mm_ch1
    mm_ch2
    mm_ch3
    ch1_motif
    ch2_motif
    ch3_motif
    seq_id
end
#function to identify and analyse ch regions
function ch_assign(seq_list, ch1, ch2, ch3)
    mm_ch1_list = []
    mm_ch2_list = []
    mm_ch3_list = []
    ch1_seq_list= []
    ch2_seq_list= []
    ch3_seq_list = []
    id_list = []

    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-15, gap_extend=-1)

    for i in seq_list

        seq_query = FASTX.sequence(i)
        seq_id = identifier(i)
        gene = split(identifier(i), "_")[1]
        mm_ch1 = nothing
        mm_ch2 = nothing
        mm_ch3 = nothing
        ch1_seq = nothing
        ch2_seq = nothing
        ch3_seq = nothing
        for x in ch1 
            if !startswith(identifier(x), gene)
                continue
            end
            seq_ref = FASTX.sequence(x)
            aln = alignment(pairalign(LocalAlignment(), seq_ref, seq_query, scoremodel))
            ch1_seq = LongDNA{4}([y for (_, y) in aln])
            mm_ch1 = count_mismatches(aln)
        end
        push!(mm_ch1_list, mm_ch1)
        push!(ch1_seq_list, ch1_seq)

        for x in ch2
            if !startswith(identifier(x), gene)
                continue
            end
            seq_ref = FASTX.sequence(x)
            aln = alignment(pairalign(LocalAlignment(), seq_ref, seq_query, scoremodel))
            ch2_seq = LongDNA{4}([y for (_, y) in aln])
            mm_ch2 = count_mismatches(aln)
        end
        push!(mm_ch2_list, mm_ch2)
        push!(ch2_seq_list, ch2_seq)

        for x in ch3
            if !startswith(identifier(x), gene)
                continue
            end
            seq_ref = FASTX.sequence(x)
            aln = alignment(pairalign(LocalAlignment(), seq_ref, seq_query, scoremodel))
            ch3_seq = LongDNA{4}([y for (_, y) in aln])
            mm_ch3= count_mismatches(aln)
        end
        push!(mm_ch3_list, mm_ch3)
        push!(ch3_seq_list, ch3_seq)


        push!(id_list, seq_id)
    end
    return ch_array(mm_ch1_list,  
                    mm_ch2_list, 
                    mm_ch3_list, 
                    ch1_seq_list, 
                    ch2_seq_list, 
                    ch3_seq_list,
                    id_list)
end

function unique_seqs(seq_list)
    cap_dict = Dict()
    unique_s = []
    for i in seq_list
        push!(unique_s, FASTX.sequence(i))
    end
    unique!(unique_s)
    unique_records = []
    seq_num = 1
    for i in unique_s
        cap_list = []
        size = 0
        gene = true
        for x in seq_list
            if i == FASTX.sequence(x)
                if gene == true
                    gene = split(split(identifier(x), "_")[2], "/")[1]
                end
                cap = split(identifier(x), "_")[1]
                seq_size = parse(Int64, split(identifier(x), "_")[3])
                size = size + seq_size
                if cap ∉ cap_list
                    push!(cap_list, cap)
                end
            end
        end
        cap_num = length(cap_list)
        new_id = join([gene, string(seq_num), string(cap_num), string(size)], "_") 
        new_rec = FASTA.Record(new_id, i)
        push!(unique_records, new_rec)
        seq_num = seq_num + 1
        cap_dict[new_id] = cap_list
    end
    return unique_records, cap_dict
end

function get_ch_match(gene, part, id)
    ch_match = "none"
    if gene == "IGHG1"
        for z in 1:length(ighg1_ch_array.seq_id)
            if ighg1_ch_array.seq_id[z] == id
                if part == "CH1"
                    ch_match = ighg1_ch_array.closest_match_ch1[z]
                elseif part == "CH2"
                    ch_match = ighg1_ch_array.closest_match_ch2[z]
                elseif  part == "CH3"
                    ch_match = ighg1_ch_array.closest_match_ch3[z]
                end
            end
        end
    elseif gene == "IGHG2"
        for z in 1:length(ighg2_ch_array.seq_id)
            if ighg2_ch_array.seq_id[z] == id
                if part == "CH1"
                    ch_match = ighg2_ch_array.closest_match_ch1[z]
                elseif part == "CH2"
                    ch_match = ighg2_ch_array.closest_match_ch2[z]
                elseif  part == "CH3"
                    ch_match = ighg2_ch_array.closest_match_ch3[z]
                end
            end
        end
    elseif gene == "IGHG3"
        for z in 1:length(ighg3_ch_array.seq_id)
            if ighg3_ch_array.seq_id[z] == id
                if part == "CH1"
                    ch_match = ighg3_ch_array.closest_match_ch1[z]
                elseif part == "CH2"
                    ch_match = ighg3_ch_array.closest_match_ch2[z]
                elseif  part == "CH3"
                    ch_match = ighg3_ch_array.closest_match_ch3[z]
                end
            end
        end
    elseif gene == "IGHG4"
        for z in 1:length(ighg4_ch_array.seq_id)
            if ighg4_ch_array.seq_id[z] == id
                if part == "CH1"
                    ch_match = ighg4_ch_array.closest_match_ch1[z]
                elseif part == "CH2"
                    ch_match = ighg4_ch_array.closest_match_ch2[z]
                elseif  part == "CH3"
                    ch_match = ighg4_ch_array.closest_match_ch3[z]
                end
            end
        end
    end
    return ch_match
end

function ch_array_to_tab(sample_array)
    num_pp_list = []
    size_list = []
    seq_id_list = []

    for i in sample_array.seq_id
        parts = split(i, "_")
        push!(seq_id_list, parts[1]*"_"*parts[2])
        push!(num_pp_list,parts[3])
        push!(size_list, parts[4])
    end
    array_tab = DataFrame(
        seq_id = seq_id_list,
        cap_num = num_pp_list,
        read_num = size_list,
        ch1_seq = sample_array.ch1_motif,
        ch1_snv = sample_array.mm_ch1,
        ch2_seq = sample_array.ch2_motif,
        ch2_snv = sample_array.mm_ch2,
        ch3_seq = sample_array.ch3_motif,
        ch3_snv = sample_array.mm_ch3,
    )

    allele_call = []
    for row in eachrow(array_tab)
        novel = false
        while novel == false
            if row.ch1_snv > 0
                novel = true
            elseif row.ch2_snv > 0
                novel = true
            elseif row.ch3_snv > 0
                novel = true
            else
                break
            end
        end

        if novel
            push!(allele_call, "novel")
        else
            push!(allele_call, "imgt")
        end
    end
    
    array_tab.allele = allele_call
    return array_tab

end

function fasta_output(file_name, seq_list)
    x = open(FASTA.Writer, file_name)
    for i in seq_list
        write(x, i)
    end
    close(x)
end

function generate_nucleotide_city(seq_list, ref_seq)
    pos_dict = Dict()
    pos = []
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-15, gap_extend=-1)
    ref_seq = FASTA.sequence(ref_seq)
    n = 1
    for i in seq_list
        print("sequence "*string(n))
        seq = FASTA.sequence(i)
        aln = alignment(pairalign(LocalAlignment(), seq, ref_seq, scoremodel))
        for x in 1:length(collect(aln))
            println(x)
            if collect(aln)[x][1] != collect(aln)[x][2] 
                if x ∉ keys(pos_dict)
                    pos_dict[x] = [collect(aln)[x][1]]
                    push!(pos, x)
                end
                if collect(aln)[x][1] ∉ pos_dict[x]
                    push!(pos_dict[x], collect(aln)[x][1])
                end
            else
                if x ∉ keys(pos_dict)
                    pos_dict[x] = []
                    push!(pos, x)
                end
            end
        end
        n = n + 1
    end
    pos_value = []
    for i in values(pos_dict)
        push!(pos_value, length(i))
    end
    pos_tab = DataFrame(position = pos, nucl = pos_value)
    return pos_tab
end

function align_ref(ref_seq, seq_list, ref_name)
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-15, gap_extend=-1)
    aln = alignment(pairalign(LocalAlignment(), FASTA.sequence(seq_list[1]), FASTA.sequence(ref_seq), scoremodel))
    trimmed_ref = LongDNA{4}([y for (_, y) in aln])
    ref_rec = FASTA.Record(ref_name,trimmed_ref)
    out_list = vcat(ref_rec, seq_list)
    return out_list
end

function freq_create(ref_seq,seq_list, ref_name)
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-15, gap_extend=-1)
    aln = alignment(pairalign(LocalAlignment(), FASTA.sequence(seq_list[1]), FASTA.sequence(ref_seq), scoremodel))
    trimmed_ref = LongDNA{4}([y for (_, y) in aln])
    ref_rec = FASTA.Record(ref_name,trimmed_ref)
    new_seq_list = []
    for i in seq_list
        seq_size = parse(Int, split(identifier(i), "_")[4])
        for x in 1:seq_size
            push!(new_seq_list, FASTA.Record(identifier(i)*"_"*string(x), FASTA.sequence(i)))
        end
    end
    out_list = vcat(ref_rec, new_seq_list)
    return out_list
end
################################################################################################################

#compile all participant IGHG seqs and remove pseudogenes
ighg1 = compile_cap_seqs("/Users/alaine/Data/Constant_heavy/holly_test_data/IGHG1_sorted.fa")
ighg2 = compile_cap_seqs("/Users/alaine/Data/Constant_heavy/holly_test_data/IGHG2_sorted.fa")
ighg3 = compile_cap_seqs("/Users/alaine/Data/Constant_heavy/holly_test_data/IGHG3_sorted.fa")
ighg4 = compile_cap_seqs("/Users/alaine/Data/Constant_heavy/holly_test_data/IGHG4_sorted.fa")

#filter the sequences
ighg1_filtered_seqs = two_seq_cap(ighg1)
ighg2_filtered_seqs = two_seq_cap(ighg2)
ighg3_filtered_seqs = two_seq_cap(ighg3)
ighg4_filtered_seqs = two_seq_cap(ighg4)

#collapse unique alleles
ighg1_unique_seqs = unique_seqs(ighg1_filtered_seqs)
ighg2_unique_seqs = unique_seqs(ighg2_filtered_seqs)
ighg3_unique_seqs = unique_seqs(ighg3_filtered_seqs)
ighg4_unique_seqs = unique_seqs(ighg4_filtered_seqs)

fasta_output("/Users/alaine/Data/Constant_heavy/ighg1_unique_seqs.fa", ighg1_unique_seqs[1])
fasta_output("/Users/alaine/Data/Constant_heavy/ighg2_unique_seqs.fa", ighg2_unique_seqs[1])
fasta_output("/Users/alaine/Data/Constant_heavy/ighg3_unique_seqs.fa", ighg3_unique_seqs[1])
fasta_output("/Users/alaine/Data/Constant_heavy/ighg4_unique_seqs.fa", ighg4_unique_seqs[1])

ighg1_seq_caps = ighg1_unique_seqs[2]
ighg2_seq_caps = ighg2_unique_seqs[2]
ighg3_seq_caps = ighg3_unique_seqs[2]
ighg4_seq_caps = ighg4_unique_seqs[2]


#extract and identify novel mutations in CH regions
ch_total = compile_seqs("/Users/alaine/Refs/Constant/ighg_ch_refs.fa")

ch1 = []
ch2 = []
ch3 = []

for i in ch_total
    gene = split(identifier(i), "|")[2]
    part = split(description(i), "|")[5]
    new_id = gene*"_"*part
    new_seq = FASTA.Record(new_id, uppercase(FASTX.sequence(i)))
    if occursin("CH1", identifier(new_seq))
        push!(ch1, new_seq)
    elseif occursin("CH2", identifier(new_seq))
        push!(ch2, new_seq)
    else
        push!(ch3, new_seq)
    end
end

ighg1_ch_array = ch_assign(ighg1_unique_seqs[1], ch1, ch2, ch3)
ighg2_ch_array = ch_assign(ighg2_unique_seqs[1], ch1, ch2, ch3)
ighg3_ch_array = ch_assign(ighg3_unique_seqs[1], ch1, ch2, ch3)
ighg4_ch_array = ch_assign(ighg4_unique_seqs[1], ch1, ch2, ch3)

# save ch_arrays

jldsave("/Users/alaine/Data/Constant_heavy/ighg1_ch_array.jld"; ighg1_ch_array)
jldsave("/Users/alaine/Data/Constant_heavy/ighg2_ch_array.jld"; ighg2_ch_array)
jldsave("/Users/alaine/Data/Constant_heavy/ighg3_ch_array.jld"; ighg3_ch_array)
jldsave("/Users/alaine/Data/Constant_heavy/ighg4_ch_array.jld"; ighg4_ch_array)

# load ch_arrays

# ighg1_ch_array = load("/Users/alaine/Data/Constant_heavy/ighg1_ch_array.jld", "ighg1_ch_array")
# ighg2_ch_array = load("/Users/alaine/Data/Constant_heavy/ighg2_ch_array.jld", "ighg2_ch_array")
# ighg3_ch_array = load("/Users/alaine/Data/Constant_heavy/ighg3_ch_array.jld", "ighg3_ch_array")
# ighg4_ch_array = load("/Users/alaine/Data/Constant_heavy/ighg4_ch_array.jld", "ighg4_ch_array")

# save ch's into separate fastas

ch_seqs = Dict()

for i in [ighg1_ch_array, ighg2_ch_array, ighg3_ch_array, ighg4_ch_array]
    gene = split(i.closest_match_ch1[1], "*")[1]
    num_seqs = length(i.closest_match_ch1)
    ch1_seqs = []
    ch2_seqs = []
    ch3_seqs = []
    seq_ids = []
    for x in 1:num_seqs
        ch1_s = i.ch1_motif[x]
        ch2_s = i.ch2_motif[x]
        ch3_s = i.ch3_motif[x]
        ch1_mm = i.mm_ch1[x]
        ch2_mm = i.mm_ch2[x]
        ch3_mm = i.mm_ch3[x]
        ch1_seq = FASTX.FASTA.Record(i.seq_id[x], ch1_s)
        ch2_seq = FASTX.FASTA.Record(i.seq_id[x], ch2_s)
        ch3_seq = FASTX.FASTA.Record(i.seq_id[x], ch3_s)
        push!(ch1_seqs, ch1_seq)
        push!(ch2_seqs, ch2_seq)
        push!(ch3_seqs, ch3_seq)
    end

    w=open(FASTA.Writer, "/Users/alaine/Data/Constant_heavy/"*gene*"_ch1.fa")
    for i in ch1_seqs
        write(w, i)
    end
    close(w)

    w=open(FASTA.Writer, "/Users/alaine/Data/Constant_heavy/"*gene*"_ch2.fa")
    for i in ch2_seqs
        write(w, i)
    end
    close(w)

    w=open(FASTA.Writer, "/Users/alaine/Data/Constant_heavy/"*gene*"_ch3.fa")
    for i in ch3_seqs
        write(w, i)
    end
    close(w)

    ch_seqs[gene*"_CH1"] = ch1_seqs
    ch_seqs[gene*"_CH2"] = ch2_seqs
    ch_seqs[gene*"_CH3"] = ch3_seqs
end

# translate seqs

ch_seqs_aa = Dict()

for (k,v) in ch_seqs
    gene = k
    trans_list = []
    for i in v
        println(k*" "*identifier(i))
        if identifier(i) == "IGHG2_2_1_3"
            continue
        elseif identifier(i) == "IGHG3_3_1_3"
            continue
        end
        s = FASTA.sequence(i)
        s = replace(s, "-" => "")
        if endswith(gene, "CH3")
            trans = translate(LongDNA{4}(s)[3:length(s)])
        else
            trans = translate(LongDNA{4}(s)[3:length(s)-1])
        end
        push!(trans_list, FASTX.FASTA.Record(identifier(i), trans))
    end
    ch_seqs_aa[gene] = trans_list
end

# load in ch Refs

ch_aa_refs  = compile_seqs("/Users/alaine/Refs/Constant/aa_ighg_refs.fa")

filtered_aa_refs = []
for i in ch_aa_refs
    if occursin("CH",description(i))
        push!(filtered_aa_refs, i)
    end
end

# save translated ch into files

for (k,v) in ch_seqs_aa
    w=open(FASTA.Writer, "/Users/alaine/Data/Constant_heavy/"*k*"_aa.fa")
    gene = split(k, "_")[1]
    part = split(k, "_")[2]
    for i in ch_aa_refs
        if occursin(gene, identifier(i))
            if occursin(part, description(i))
                write(w, i)
            end
        end
    end
    for i in v
        write(w, i)
    end
    close(w)
end

# aa changes



aa_changes_dict = Dict()
scoremodel_prot = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1)
for (k, v) in ch_seqs_aa
    gene_array = []
    gene = split(k, "_")[1]
    part = split(k, "_")[2]
    for i in v
        seq1 = FASTX.FASTA.sequence(i)
        aa_changes = []
        aa_array = []
        id = identifier(i)
        ch_match = get_ch_match(gene, part, id)
        print(ch_match)
        for x in filtered_aa_refs
            if split(identifier(x), "|")[2] == ch_match
                if occursin(part, description(x))
                    seq2 = FASTA.sequence(x)
                    aln = alignment(pairalign(LocalAlignment(), seq1, seq2, scoremodel_prot))
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
                            mut = string(y[2])*string(n+118)*string(y[1])
                            push!(aa_changes, mut)
                        end
                        last_aa = string(y[2])*string(n)
                        n = n + 1
                    end
                end
            end
        end

        if isempty(aa_changes)
            push!(gene_array, [identifier(i)*"_"*part, "none"])
        else
            push!(gene_array, [identifier(i)*"_"*part, join(aa_changes, ", ")])
        end
        push!(gene_array)
    end
    aa_changes_dict[k] = gene_array
end

IGHG1_CH2_IGHG1_1_1_12


################################################################################################################

# construct final tabular outputs: 1 - genotypes table 2 - sequence feature table

cap_list = []
for i in [ighg1_seq_caps, ighg2_seq_caps, ighg3_seq_caps, ighg4_seq_caps]
    for (k,v) in i
        for x in v
            if x ∉ cap_list
                push!(cap_list, x)
            end
        end
    end
end

ighg1_alleles = []
ighg2_alleles = []
ighg3_alleles = []
ighg4_alleles = []
for i in cap_list
    g1_alls = []
    g2_alls = []
    g3_alls = []
    g4_alls = []
    for x in [ighg1_seq_caps, ighg2_seq_caps, ighg3_seq_caps, ighg4_seq_caps]
        for (k,v) in x
            if i ∈ v
                if occursin("IGHG1", k)
                    push!(g1_alls, k)
                elseif occursin("IGHG2", k)
                    push!(g2_alls, k)
                elseif occursin("IGHG3", k)
                    push!(g3_alls, k)
                elseif occursin("IGHG4", k)
                    push!(g4_alls, k)
                end
            end
        end
    end
    push!(ighg1_alleles, g1_alls)
    push!(ighg2_alleles, g2_alls)
    push!(ighg3_alleles, g3_alls)
    push!(ighg4_alleles, g4_alls)
end

genotab = DataFrame()

genotab.cap_participant = cap_list
genotab.ighg1 = ighg1_alleles
genotab.ighg2 = ighg2_alleles
genotab.ighg3 = ighg3_alleles
genotab.ighg4 = ighg4_alleles

CSV.write("/Users/alaine/Data/Constant_heavy/ighg_geno_tab.tsv", genotab, delim="\t")

ighg1_tab = ch_array_to_tab(ighg1_ch_array)
ighg2_tab = ch_array_to_tab(ighg2_ch_array)
ighg3_tab = ch_array_to_tab(ighg3_ch_array)
ighg4_tab = ch_array_to_tab(ighg4_ch_array)

ighg1_tab.full_sequence = FASTA.sequence.(ighg1_unique_seqs[1])
ighg2_tab.full_sequence = FASTA.sequence.(ighg2_unique_seqs[1])
ighg3_tab.full_sequence = FASTA.sequence.(ighg3_unique_seqs[1])
ighg4_tab.full_sequence = FASTA.sequence.(ighg4_unique_seqs[1])


CSV.write("/Users/alaine/Data/Constant_heavy/ighg1_seq_tab.csv", ighg1_tab)
CSV.write("/Users/alaine/Data/Constant_heavy/ighg2_seq_tab.csv", ighg2_tab)
CSV.write("/Users/alaine/Data/Constant_heavy/ighg3_seq_tab.csv", ighg3_tab)
CSV.write("/Users/alaine/Data/Constant_heavy/ighg4_seq_tab.csv", ighg4_tab)


igg_refs = compile_seqs("/Users/alaine/Refs/Constant/IgG_full_refs.fa")

igg01 = []

for i in igg_refs
    if occursin("*01",identifier(i))
        push!(igg01, i)
    end
end

ighg1_pos_tab = generate_nucleotide_city(ighg1_unique_seqs[1], igg01[1])
ighg2_pos_tab = generate_nucleotide_city(ighg2_unique_seqs[1], igg01[2])
ighg3_pos_tab = generate_nucleotide_city(ighg3_unique_seqs[1], igg01[3])
ighg4_pos_tab = generate_nucleotide_city(ighg4_unique_seqs[1], igg01[4])

CSV.write("/Users/alaine/Data/Constant_heavy/ighg1_pos_tab.csv", ighg1_pos_tab)
CSV.write("/Users/alaine/Data/Constant_heavy/ighg2_pos_tab.csv", ighg2_pos_tab)
CSV.write("/Users/alaine/Data/Constant_heavy/ighg3_pos_tab.csv", ighg3_pos_tab)
CSV.write("/Users/alaine/Data/Constant_heavy/ighg4_pos_tab.csv", ighg4_pos_tab)

ighg1_align = align_ref(igg01[1], ighg1_unique_seqs[1], "IGHG1*01_ref")
ighg2_align = align_ref(igg01[2], ighg2_unique_seqs[1], "IGHG2*01_ref")
ighg3_align = align_ref(igg01[3], ighg3_unique_seqs[1], "IGHG3*01_ref")
ighg4_align = align_ref(igg01[4], ighg4_unique_seqs[1], "IGHG4*01_ref")

fasta_output("/Users/alaine/Data/Constant_heavy/ighg1_final_seqs.fa", igh1_align)
fasta_output("/Users/alaine/Data/Constant_heavy/ighg2_final_seqs.fa", ighg2_align)
fasta_output("/Users/alaine/Data/Constant_heavy/ighg3_final_seqs.fa", ighg3_align)
fasta_output("/Users/alaine/Data/Constant_heavy/ighg4_final_seqs.fa", ighg4_align)

ighg1_freq = freq_create(igg01[1], ighg1_unique_seqs[1], "IGHG1*01_ref")
ighg2_freq = freq_create(igg01[2], ighg2_unique_seqs[1], "IGHG2*01_ref")
ighg3_freq = freq_create(igg01[3], ighg3_unique_seqs[1], "IGHG3*01_ref")
ighg4_freq = freq_create(igg01[4], ighg4_unique_seqs[1], "IGHG4*01_ref")

fasta_output("/Users/alaine/Data/Constant_heavy/ighg1_freq_seqs.fa", ighg1_freq)
fasta_output("/Users/alaine/Data/Constant_heavy/ighg2_freq_seqs.fa", ighg2_freq)
fasta_output("/Users/alaine/Data/Constant_heavy/ighg3_freq_seqs.fa", ighg3_freq)
fasta_output("/Users/alaine/Data/Constant_heavy/ighg4_freq_seqs.fa", ighg4_freq)