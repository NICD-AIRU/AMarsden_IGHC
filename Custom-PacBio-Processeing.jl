# Processing of raw PB CSS reads in FASTQ format with Julia
# A. Marsden
#
# Steps:
#   1. Remove small sequences
#   2. Fix orientation
#   3. Trim adapters
#   4. Dereplicate
#   5. Output sequences with size into new FASTA file

using FASTX, StatsBase, OrderedCollections, BioSequences

dir = "../../Data/pacbio_june_2022"

files = readdir(dir)

reverse_primers = FASTA.Reader(open("../../Refs/Primers/germline_reverse.fa", "r"))

primers_r = []
for i in reverse_primers
    push!(primers_r, i)
end


adapter_f = ExactSearchQuery(dna"gcagtcgaacatgtagctgactcaggtcac")
adapter_r = ExactSearchQuery(reverse_complement(dna"tggatcacttgtgcaagcatcacatcgtag"))


for f in files
    if endswith(f, "fastq")
        print("Processing: "*f*"\n")
        r = FASTQ.Reader(open(dir*"/"*f, "r"))
        seqs = []
        cap = split(f, ".")[2]
        mkdir(dir*"/"*cap)
        for i in r
            push!(seqs, sequence(i))
        end
        seq_len = length(seqs)
        print("Initial Sequence Count: "*string(seq_len)*"\n")


        seqs_filter = []
        m = 0
        for i in seqs
            if length(i) < 500     
                m = m + 1
                continue
            else
                push!(seqs_filter, i)
            end

        end
        print(string(m)*" Sequences Removed for being under 500bp \n")

        print("Detecting Sequences in Reverse Complement \n")
        y = 0
        for i in seqs_filter
            for primer_r in primers_r
                primer = ExactSearchQuery(LongDNA{4}(string(sequence(primer_r))))
                if occursin(primer, i)
                    reverse_complement!(i)
                    y = y + 1
                end
            end
        end
        print(string(y)*" Sequence(s) Reverse Complemented \n")

        print("Trimming Adapters from Sequences \n")
        for i in seqs_filter
            if occursin(adapter_f, i)
                deleteat!(i,findfirst(adapter_f, i))
            end
            if occursin(adapter_r,i)
                deleteat!(i,findfirst(adapter_r, i))
            end
        end
        print("Collpasing Sequences \n")
        derep_sort = sort(countmap(seqs_filter), byvalue=true, rev=true)
        count = values(derep_sort)

        print("Sequence Count Statistics \n")
        print("Min: "*string(minimum(count))*"\n")
        print("Max: "*string(maximum(count))*"\n")
        print("Mean: "*string(round(mean(count)))*"\n")
        print("Total number of Unique Sequences: "*string(length(derep_sort))*"\n")

        derep_seqs = filter(kv -> kv.second > 0, derep_sort)
        discard = length(derep_sort) - length(derep_seqs)
        print(string(discard)*" Low Abundance Sequences Discarded \n")

        print("Writing Output File with "*string(length(derep_seqs))*" Sequence(s) \n")

        f_name = join([cap, "Processed.fa"], "_")
        w=open(FASTA.Writer, dir*"/"*f_name)
        
        t = 1
        for (key, value) in derep_seqs
            seq_number = string(t)
            size = string(value)
            seq_name = join([cap, seq_number, size], "_")
            rec = FASTA.Record(seq_name, key)
            write(w, rec)
            t = t + 1
        end

        close(w)
        close(r)
        mv(dir*"/"*f, dir*"/"*cap*"/"*f)
        mv(dir*"/"*f_name, dir*"/"*cap*"/"*f_name)

        print("Processing Complete \n")
        print("--------------------------------------------------- \n")
    end
end


