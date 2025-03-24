# %%
import pyranges as pr
import pandas as pd
import numpy as np
import sqlite3
# import gffutils
import time
# import pybedtools
import os
from IPython import get_ipython
from multiprocessing import Pool
import ctypes

# %%
gtf_file_name = "Homo_sapiens.GRCh38.112.chr.gtf"
# gtf_file_name = "/home/U210050044/data/DH607-project/Homo_sapiens.GRCh38.112.chr.gtf"
# sql_db_name = "db-human.sqlite3"
sql_db_name = "file:memdb1?mode=memory&cache=shared"
results_file_name = "results_human.txt"
lib = ctypes.CDLL("./gtf_to_sql.so")
# Define the function signature
lib.run_gtf_to_sql.argtypes = [
    ctypes.c_char_p,  # db_name
    ctypes.c_char_p,  # table_name
    ctypes.c_char_p,  # input_file
    ctypes.c_int,     # num_threads
]

# %%
results_file = open(results_file_name,"w")

# %%
start_time = time.time()
human_gr = pr.read_gtf(gtf_file_name)
end_time = time.time()
results_file.write(f"Finished reading into pyranges object in time {end_time-start_time}\n")

# %%
# start_time = time.time()
# human_db = gffutils.create_db("/home/U210050044/data/DH607-project/Homo_sapiens.GRCh38.112.chr.gtf","/home/U210050044/data/DH607-project/human_db.sqlite3",disable_infer_genes=True,disable_infer_transcripts=True, force=True)
# end_time = time.time()
# results_file.write(f"Finished converting gtf to sql db in time {end_time-start_time}")

# %%
# # Open the database
# human_db = gffutils.FeatureDB("/home/U210050044/data/DH607-project/human_db.sqlite3")

# %%
# Call the function
os.system("g++ -std=c++11 -shared -fPIC -o gtf_to_sql.so gtf_to_sql.cpp -lsqlite3 -pthread")
start_time = time.time()
lib.run_gtf_to_sql(sql_db_name.encode("utf-8"), "human".encode("utf-8"), gtf_file_name.encode("utf-8"), 8)
end_time = time.time()
results_file.write(f"Finished reading into sql db in time {end_time-start_time}\n")

# %%
conn = sqlite3.connect(sql_db_name)
cur = conn.cursor()

# %%
start_time = time.time()
human_gr2 = pr.PyRanges(pd.read_sql("SELECT * FROM human", conn))
end_time = time.time()
results_file.write(f"Finished reading from sql to pyranges in time {end_time-start_time}\n")

# %% [markdown]
# # Aggregate Queries: Counting transcripts, calculating exon lengths, and grouping features
# ### Aggregate Query 1: Count the number of exons for each gene

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o (human_gr[human_gr.Feature == "exon"].df.groupby("gene_id").size().reset_index(name="exon_count"))')
results_file.write(f"Time taken to count exons using pyranges is {time}\n")

# %%
# # Using gffutils
# def exon_counts_gffutls():
#     # Iterate over all genes and count their exons
#     exon_counts = {}
#     for gene in human_db.features_of_type("gene"):
#         # Efficiently count exons for each gene using the `children` method
#         exon_counts[gene.id] = sum(
#             1 for _ in human_db.children(gene, featuretype="exon")
#     )
#     return exon_counts
# %timeit exon_counts_gff = exon_counts_gffutls()

# %%
# Using sql
time = get_ipython().run_line_magic("timeit",'-o exon_counts_sql = pd.read_sql_query("SELECT gene_id, COUNT(*) as exon_count FROM human WHERE Feature = \'exon\' GROUP BY gene_id", conn)')
results_file.write(f"Time taken to count exons using sql is {time}\n")

# %%
# check if exon_counts_pr and exon_counts_sql are identical
exon_counts_pr = (human_gr[human_gr.Feature == "exon"].df.groupby("gene_id").size().reset_index(name="exon_count"))
exon_counts_sql = pd.read_sql_query("SELECT gene_id, COUNT(*) as exon_count FROM human WHERE Feature = 'exon' GROUP BY gene_id", conn)
exon_counts_sql["exon_count"] = exon_counts_sql["exon_count"].astype(int)
exon_counts_pr["exon_count"] = exon_counts_pr["exon_count"].astype(int)
if exon_counts_sql.equals(exon_counts_pr):
    results_file.write("exon counts are identical\n")
else:
    results_file.write("exon counts are not identical\n")

# %% [markdown]
# ### Aggregate Query 2: Calculate the total length of exons for each gene

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o exon_lengths_pr = (human_gr[human_gr.Feature == "exon"].df.assign(length=lambda df: df["End"] - df["Start"]).groupby("gene_id")["length"].sum().reset_index(name="total_exon_length"))')
results_file.write(f"Time taken to calculate total exon length using pyranges is {time}\n")

# %%
# # Using gffutils
# def exon_lengths_gffutls():
#     # Dictionary to store total exon length per gene
#     exon_lengths = {}
#     # Iterate over all genes and calculate total exon length
#     for gene in human_db.features_of_type("gene"):
#         # Compute total exon length for each gene
#         exon_lengths[gene.id] = sum(
#             (child.end - child.start + 1)
#             for child in human_db.children(gene, featuretype="exon")
#         )
#     return exon_lengths
# %timeit exon_lengths_gff = exon_lengths_gffutls()

# %%
# using sql
time = get_ipython().run_line_magic("timeit",'-o exon_lengths_sql = pd.read_sql_query("SELECT gene_id, SUM(End - Start) as total_exon_length FROM human WHERE Feature = \'exon\' GROUP BY gene_id", conn)')
results_file.write(f"Time taken to calculate total exon length using sql is {time}\n")

# %%
# check if exon_lengths_pr and exon_lengths_sql are identical
exon_lengths_pr = (human_gr[human_gr.Feature == "exon"].df.assign(length=lambda df: df["End"] - df["Start"]).groupby("gene_id")["length"].sum().reset_index(name="total_exon_length"))
exon_lengths_sql = pd.read_sql_query("SELECT gene_id, SUM(End - Start) as total_exon_length FROM human WHERE Feature = 'exon' GROUP BY gene_id", conn)
exon_lengths_sql["total_exon_length"] = exon_lengths_sql["total_exon_length"].astype(int)
exon_lengths_pr["total_exon_length"] = exon_lengths_pr["total_exon_length"].astype(int)
if exon_lengths_sql.equals(exon_lengths_pr):
    results_file.write("exon lengths are identical\n")
else:
    results_file.write("exon lengths are not identical\n")

# %% [markdown]
# ### Aggregate Query 3: Identify the chromosome with the highest number of transcripts

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o transcript_counts_pr = human_gr[human_gr.Feature == "transcript"].df["Chromosome"].value_counts().idxmax()')
results_file.write(f"Time taken to find chromosome with most transcripts using pyranges is {time}\n")

# %%
# # Using gffutils
# def transcript_counts_gffutls():
#     # Dictionary to store the number of transcripts per chromosome
#     transcript_counts = {}
#     # Iterate over all transcripts and count the number of transcripts per chromosome
#     for transcript in human_db.features_of_type("transcript"):
#         chrom = transcript.chrom
#         if chrom not in transcript_counts:
#             transcript_counts[chrom] = 0
#         transcript_counts[chrom] += 1
#     return max(transcript_counts, key=transcript_counts.get)
# %timeit transcript_counts_gff = transcript_counts_gffutls()

# %%
# using sql
time = get_ipython().run_line_magic("timeit",'-o transcript_counts_sql = pd.read_sql_query("SELECT Chromosome, COUNT(*) as transcript_count FROM human WHERE Feature = \'transcript\' GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn)')
results_file.write(f"Time taken to find chromosome with most transcripts using sql is {time}\n")

# %%
# verify if transcript_counts_pr and transcript_counts_sql are identical
transcript_counts_pr = human_gr2[human_gr2.Feature == "transcript"].df["Chromosome"].value_counts().idxmax()
transcript_counts_sql = pd.read_sql_query("SELECT Chromosome, COUNT(*) as transcript_count FROM human WHERE Feature = 'transcript' GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn)
if transcript_counts_pr == transcript_counts_sql["Chromosome"].values[0]:
    results_file.write("transcript counts are identical\n")
else:
    results_file.write("transcript counts are not identical\n")

# %% [markdown]
# # Interval Arithmetic
# ### Query 1: Merging Overlapping Exon Intervals

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o exon_intervals_pr = human_gr[human_gr.Feature == "exon"].merge(strand=True)')
results_file.write(f"Time taken to merge exons using pyranges is {time}\n")

# %%
# # using gffutils TODO: run
# def exon_intervals_gffutls():
#     exon_intervals = [
#         gffutils.helpers.asinterval(exon) for exon in human_db.features_of_type("exon")
#     ]
#     exon_intervals_sorted = sorted(exon_intervals, key=lambda x: (x.chrom, x.start))
#     # Use pybedtools to merge overlapping intervals
#     return pybedtools.BedTool(exon_intervals_sorted).merge()
# %timeit exon_intervals_gff = exon_intervals_gffutls()

# %%
# using sql
def get_exon_intervals_sql(chrom_strand):
    chrom, strand = chrom_strand
    exon_intervals_sql = pd.read_sql_query(f"SELECT Start, End FROM human WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn)
    exon_intervals_sql = pr.methods.merge._merge(exon_intervals_sql, chromosome=chrom, count=None, strand=strand)
    return exon_intervals_sql

def get_exon_intervals_sql_multi():
    chrom_strand_tup = pd.read_sql_query("SELECT DISTINCT Chromosome, Strand FROM human", conn)
    chrom_strand_tup = list(zip(chrom_strand_tup["Chromosome"], chrom_strand_tup["Strand"]))
    with Pool(8) as p:
        exon_intervals_sql = p.map(get_exon_intervals_sql, chrom_strand_tup)
    exon_intervals_sql = pd.concat(exon_intervals_sql).sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return exon_intervals_sql

# %%
time = get_ipython().run_line_magic("timeit",'-o exon_intervals_sql = get_exon_intervals_sql_multi()')
results_file.write(f"Time taken to merge exons using sql is {time}\n")

# %%
# check if exon_intervals_pr and exon_intervals_sql are identical
exon_intervals_pr = human_gr[human_gr.Feature == "exon"].merge(strand=True).df.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
exon_intervals_sql = get_exon_intervals_sql_multi()
if exon_intervals_pr.equals(exon_intervals_sql):
    results_file.write("exon intervals are identical\n")
else:
    results_file.write("exon intervals are not identical\n")

# %% [markdown]
# ### Query 2: Finding Overlaps with a Specific Interval
# Find all gene features that overlap a given interval chr1:100000-200000 + strand

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o overlapping_genes_pr = human_gr[human_gr.Feature == "gene"].overlap(pr.from_dict({"Chromosome": ["1"], "Start": [100000], "End": [200000], "Strand": ["+"]}), strandedness="same")')
results_file.write(f"Time taken to find overlapping genes using pyranges is {time}\n")

# %%
# # using gffutils TODO: run
# def overlapping_genes_gffutls():
#     exon_intervals_overlapping = pybedtools.BedTool.from_dataframe(pd.DataFrame({"Chromosome": ["1"], "Start": [100000], "End": [200000]}))
#     exon_intervals = [
#         gffutils.helpers.asinterval(exon) for exon in human_db.features_of_type("exon")
#     ]
#     exon_intervals_sorted = sorted(exon_intervals, key=lambda x: (x.chrom, x.start))
#     # Use pybedtools to merge overlapping intervals
#     return exon_intervals_overlapping.window(pybedtools.BedTool(exon_intervals_sorted).merge(), w=0)
# %timeit overlapping_genes_gff = overlapping_genes_gffutls()

# %%
# using sql
other_genes_all = pd.DataFrame({"Chromosome": ["1"], "Start": [100000], "End": [200000], "Strand": ["+"]})
def get_overlapping_genes_sql(chrom_strand):
    chrom, strand = chrom_strand
    self_genes_sql = pd.read_sql_query(f"SELECT * FROM human WHERE Feature = 'gene' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn)
    other_genes = other_genes_all[(other_genes_all["Chromosome"] == chrom) & (other_genes_all["Strand"] == strand)]
    overlapping_genes_sql = pr.methods.intersection._overlap(self_genes_sql, other_genes, how="first")
    return overlapping_genes_sql

def get_overlapping_genes_sql_multi():
    chrom_strand_tup = pd.read_sql_query("SELECT DISTINCT Chromosome, Strand FROM human", conn)
    chrom_strand_tup = list(zip(chrom_strand_tup["Chromosome"], chrom_strand_tup["Strand"]))
    with Pool(8) as p:
        overlapping_genes_sql = p.map(get_overlapping_genes_sql, chrom_strand_tup)
    overlapping_genes_sql = pd.concat(overlapping_genes_sql) #.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return overlapping_genes_sql

# %%
time = get_ipython().run_line_magic("timeit",'-o overlapping_genes_sql = get_overlapping_genes_sql_multi()')
results_file.write(f"Time taken to find overlapping genes using sql is {time}\n")

# %%
# check if overlapping_genes_pr and overlapping_genes_sql are identical
overlapping_genes_pr = human_gr[human_gr.Feature == "gene"].overlap(pr.from_dict({"Chromosome": ["1"], "Start": [100000], "End": [200000], "Strand": ["+"]}), strandedness="same")
overlapping_genes_sql = get_overlapping_genes_sql_multi()
overlapping_genes_sql = overlapping_genes_sql.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
overlapping_genes_pr = overlapping_genes_pr.df.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
diff = overlapping_genes_pr.compare(overlapping_genes_sql)
if diff.empty:
    results_file.write("overlapping genes are identical\n")
else:
    results_file.write("overlapping genes are not identical\n")

# %% [markdown]
# ### Query 3: Subtracting Intervals
# Remove a set of repetitive regions from the exon features.
# 
# chr1:15M-16M,20M-21M

# %%
# # Using pyranges
# %timeit subtracted_intervals_pr = human_gr[human_gr.Feature == "exon"].merge(strand=False).subtract(pr.from_dict({"Chromosome": ["1", "1"],"Start": [15000000, 20000000],"End": [16000000, 21000000],}))

# %%
# # using gffutils TODO: run
# def subtracted_intervals_gffutls():
#     exon_intervals_subtract = pybedtools.BedTool.from_dataframe(pd.DataFrame({"Chromosome": ["1", "1"], "Start": [15000000, 20000000], "End": [16000000, 21000000]}))
#     exon_intervals = [
#         gffutils.helpers.asinterval(exon) for exon in human_db.features_of_type("exon")
#     ]
#     exon_intervals_sorted = sorted(exon_intervals, key=lambda x: (x.chrom, x.start))
#     # Use pybedtools to merge overlapping intervals
#     return pybedtools.BedTool(exon_intervals_sorted).merge().subtract(exon_intervals_subtract)
# %timeit overlapping_genes_gff = overlapping_genes_gffutls()

# %%
# # using sql
# def find_subtracted_intervals_sql():
#     exon_intervals_subtract = pybedtools.BedTool.from_dataframe(pd.DataFrame({"Chromosome": ["1", "1"], "Start": [15000000, 20000000], "End": [16000000, 21000000]}))
#     exon_intervals = pd.read_sql_query("SELECT Chromosome, Start, End FROM human WHERE Feature = 'exon' ORDER BY Chromosome, Start", conn)
#     # Use pybedtools to merge overlapping intervals
#     return pybedtools.BedTool.from_dataframe(exon_intervals).merge().subtract(exon_intervals_subtract)
# %timeit subtracted_intervals_sql = find_subtracted_intervals_sql()

# %%
results_file.close()
conn.close()


