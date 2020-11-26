#! usr/bin/env julia

#load packages
using ArgParse
using StringDistances
using DataFrames
using CSV

#set up command line arguments
s = ArgParseSettings()
@add_arg_table! s begin
"--Input"
	required = true
	help = "Directory of target sequences"
"--Reference"
	help = "Reference consensus sequence"
"--Sub"
	arg_type = Int
	default = 3
	help = "length of substring"
"--Output"
	required = true
	help = "filename to store distance matrix as a csv"
end

parsed_args = parse_args(ARGS, s)

Input = parsed_args["Input"]
Ref = parsed_args["Reference"]
Subs = parsed_args["Sub"]
Output = parsed_args["Output"]

r = open(FASTX.FASTA.Reader, Ref)

#load input files
df = cd(readdir, Input)

#pairwise comparison of each string using Overlap coeffient
for i in length(df)
	q = open(FASTX.FASTA.Reader, df[i])
	overlap = evaluate(Overlap(Subs), q, r)
	global Zero[1,i] = overlap
	end
end


#convert back to dataframe
Zero_df = DataFrame(Zero)

#create output csv
CSV.write(Output, Zero_df)
