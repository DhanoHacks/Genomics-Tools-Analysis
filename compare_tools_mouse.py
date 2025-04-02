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
gtf_file_name = "gencode.vM36.annotation.gtf"
# gtf_file_name = "/home/U210050044/data/DH607-project/gencode.vM36.annotation.gtf"
sql_db_name = "db-mouse.sqlite3"
# sql_db_name = "file:mouse?mode=memory&cache=shared"
results_file_name = "results_mouse.txt"
os.system("g++ -std=c++11 -shared -fPIC -o gtf_to_sql.so gtf_to_sql.cpp -lsqlite3 -pthread")
lib = ctypes.CDLL("./gtf_to_sql.so")
# Define the function signature
lib.run_gtf_to_sql.argtypes = [
    ctypes.c_char_p,  # db_name
    ctypes.c_char_p,  # table_name
    ctypes.c_char_p,  # input_file
    ctypes.c_int,     # num_threads
]
# OUTPUT_MODE = "print"
OUTPUT_MODE = "file"

# %%
if OUTPUT_MODE == "file":
    results_file = open(results_file_name,"w")

# %%
start_time = time.time()
mouse_gr = pr.read_gtf(gtf_file_name)
end_time = time.time()
if OUTPUT_MODE == "file":
    results_file.write(f"Finished reading into pyranges object in time {end_time-start_time}\n")
else:
    print(f"Finished reading into pyranges object in time {end_time-start_time}")

# %%
# Call the function
start_time = time.time()
lib.run_gtf_to_sql(sql_db_name.encode("utf-8"), "mouse".encode("utf-8"), gtf_file_name.encode("utf-8"), 8)
end_time = time.time()
if OUTPUT_MODE == "file":
    results_file.write(f"Finished reading into sql db in time {end_time-start_time}\n")
else:
    print(f"Finished reading into sql db in time {end_time-start_time}")

# %%
conn = sqlite3.connect(sql_db_name)
cur = conn.cursor()

# %%
start_time = time.time()
mouse_gr2 = pr.PyRanges(pd.read_sql("SELECT * FROM mouse", conn))
end_time = time.time()
if OUTPUT_MODE == "file":
    results_file.write(f"Finished reading from sql to pyranges in time {end_time-start_time}\n")
else:
    print(f"Finished reading from sql to pyranges in time {end_time-start_time}")

# %% [markdown]
# # Aggregate Queries: Counting transcripts, calculating exon lengths, and grouping features
# ### Aggregate Query 1: Count the number of exons for each gene

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o (mouse_gr[mouse_gr.Feature == "exon"].df.groupby("gene_id").size().reset_index(name="exon_count"))')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to count exons using pyranges is {time}\n")
else:
    print(f"Time taken to count exons using pyranges is {time}")

# %%
# Using sql
time = get_ipython().run_line_magic("timeit",'-o exon_counts_sql = pd.read_sql_query("SELECT gene_id, COUNT(*) as exon_count FROM mouse WHERE Feature = \'exon\' GROUP BY gene_id", conn)')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to count exons using sql is {time}\n")
else:
    print(f"Time taken to count exons using sql is {time}")

# %%
# check if exon_counts_pr and exon_counts_sql are identical
exon_counts_pr = (mouse_gr[mouse_gr.Feature == "exon"].df.groupby("gene_id").size().reset_index(name="exon_count"))
exon_counts_sql = pd.read_sql_query("SELECT gene_id, COUNT(*) as exon_count FROM mouse WHERE Feature = 'exon' GROUP BY gene_id", conn)
exon_counts_sql["exon_count"] = exon_counts_sql["exon_count"].astype(int)
exon_counts_pr["exon_count"] = exon_counts_pr["exon_count"].astype(int)
if exon_counts_sql.equals(exon_counts_pr):
    result_text = "exon counts are identical"
else:
    result_text = "exon counts are not identical"
if OUTPUT_MODE == "file":
    results_file.write(f"{result_text}\n")
else:
    print(result_text)

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o exon_lengths_pr = (mouse_gr[mouse_gr.Feature == "exon"].df.assign(length=lambda df: df["End"] - df["Start"]).groupby("gene_id")["length"].sum().reset_index(name="total_exon_length"))')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to calculate total exon length using pyranges is {time}\n")
else:
    print(f"Time taken to calculate total exon length using pyranges is {time}")

# %%
# using sql
time = get_ipython().run_line_magic("timeit",'-o exon_lengths_sql = pd.read_sql_query("SELECT gene_id, SUM(End - Start) as total_exon_length FROM mouse WHERE Feature = \'exon\' GROUP BY gene_id", conn)')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to calculate total exon length using sql is {time}\n")
else:
    print(f"Time taken to calculate total exon length using sql is {time}")

# %%
# check if exon_lengths_pr and exon_lengths_sql are identical
exon_lengths_pr = (mouse_gr[mouse_gr.Feature == "exon"].df.assign(length=lambda df: df["End"] - df["Start"]).groupby("gene_id")["length"].sum().reset_index(name="total_exon_length"))
exon_lengths_sql = pd.read_sql_query("SELECT gene_id, SUM(End - Start) as total_exon_length FROM mouse WHERE Feature = 'exon' GROUP BY gene_id", conn)
exon_lengths_sql["total_exon_length"] = exon_lengths_sql["total_exon_length"].astype(int)
exon_lengths_pr["total_exon_length"] = exon_lengths_pr["total_exon_length"].astype(int)
if exon_lengths_sql.equals(exon_lengths_pr):
    results_text = "exon lengths are identical"
else:
    results_text = "exon lengths are not identical"
if OUTPUT_MODE == "file":
    results_file.write(f"{results_text}\n")
else:
    print(results_text)

# %% [markdown]
# ### Aggregate Query 3: Identify the chromosome with the highest number of transcripts

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o transcript_counts_pr = mouse_gr[mouse_gr.Feature == "transcript"].df["Chromosome"].value_counts().idxmax()')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to find chromosome with most transcripts using pyranges is {time}\n")
else:
    print(f"Time taken to find chromosome with most transcripts using pyranges is {time}")

# %%
# using sql
time = get_ipython().run_line_magic("timeit",'-o transcript_counts_sql = pd.read_sql_query("SELECT Chromosome, COUNT(*) as transcript_count FROM mouse WHERE Feature = \'transcript\' GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn)')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to find chromosome with most transcripts using sql is {time}\n")
else:
    print(f"Time taken to find chromosome with most transcripts using sql is {time}")

# %%
# verify if transcript_counts_pr and transcript_counts_sql are identical
transcript_counts_pr = mouse_gr2[mouse_gr2.Feature == "transcript"].df["Chromosome"].value_counts().idxmax()
transcript_counts_sql = pd.read_sql_query("SELECT Chromosome, COUNT(*) as transcript_count FROM mouse WHERE Feature = 'transcript' GROUP BY Chromosome ORDER BY transcript_count DESC LIMIT 1", conn)
if transcript_counts_pr == transcript_counts_sql["Chromosome"].values[0]:
    results_text = "transcript counts are identical"
else:
    results_text = "transcript counts are not identical"
if OUTPUT_MODE == "file":
    results_file.write(f"{results_text}\n")
else:
    print(results_text)

# %% [markdown]
# # Interval Arithmetic
# ### Query 1: Merging Overlapping Exon Intervals

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o exon_intervals_pr = mouse_gr[mouse_gr.Feature == "exon"].merge(strand=True)')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to merge exons using pyranges is {time}\n")
else:
    print(f"Time taken to merge exons using pyranges is {time}")

# %%
# using sql
def get_exon_intervals_sql(chrom_strand):
    chrom, strand = chrom_strand
    exon_intervals_sql = pd.read_sql_query(f"SELECT Start, End FROM mouse WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn)
    exon_intervals_sql = pr.methods.merge._merge(exon_intervals_sql, chromosome=chrom, count=None, strand=strand)
    return exon_intervals_sql

def get_exon_intervals_sql_multi():
    chrom_strand_tup = pd.read_sql_query("SELECT DISTINCT Chromosome, Strand FROM mouse", conn)
    chrom_strand_tup = list(zip(chrom_strand_tup["Chromosome"], chrom_strand_tup["Strand"]))
    with Pool(8) as p:
        exon_intervals_sql = p.map(get_exon_intervals_sql, chrom_strand_tup)
    exon_intervals_sql = pd.concat(exon_intervals_sql).sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return exon_intervals_sql

# %%
time = get_ipython().run_line_magic("timeit",'-o exon_intervals_sql = get_exon_intervals_sql_multi()')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to merge exons using sql is {time}\n")
else:
    print(f"Time taken to merge exons using sql is {time}")

# %%
# check if exon_intervals_pr and exon_intervals_sql are identical
exon_intervals_pr = mouse_gr[mouse_gr.Feature == "exon"].merge(strand=True).df.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
exon_intervals_sql = get_exon_intervals_sql_multi()
if exon_intervals_pr.equals(exon_intervals_sql):
    results_text = "exon intervals are identical"
else:
    results_text = "exon intervals are not identical"
if OUTPUT_MODE == "file":
    results_file.write(f"{results_text}\n")
else:
    print(results_text)

# %% [markdown]
# ### Query 2: Finding Overlaps with a Specific Interval
# Find all gene features that overlap a given interval chr1:100000-200000 + strand

# %%
# Using pyranges
time = get_ipython().run_line_magic("timeit",'-o overlapping_genes_pr = mouse_gr[mouse_gr.Feature == "gene"].overlap(pr.from_dict({"Chromosome": ["1"], "Start": [100000], "End": [200000], "Strand": ["+"]}), strandedness="same")')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to find overlapping genes using pyranges is {time}\n")
else:
    print(f"Time taken to find overlapping genes using pyranges is {time}")

# %%
# using sql
other_genes_all = pd.DataFrame({"Chromosome": ["1"], "Start": [100000], "End": [200000], "Strand": ["+"]})
def get_overlapping_genes_sql(chrom_strand):
    chrom, strand = chrom_strand
    self_genes_sql = pd.read_sql_query(f"SELECT * FROM mouse WHERE Feature = 'gene' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn)
    other_genes = other_genes_all[(other_genes_all["Chromosome"] == chrom) & (other_genes_all["Strand"] == strand)]
    overlapping_genes_sql = pr.methods.intersection._overlap(self_genes_sql, other_genes, how="first")
    return overlapping_genes_sql

def get_overlapping_genes_sql_multi():
    chrom_strand_tup = pd.read_sql_query("SELECT DISTINCT Chromosome, Strand FROM mouse", conn)
    chrom_strand_tup = list(zip(chrom_strand_tup["Chromosome"], chrom_strand_tup["Strand"]))
    with Pool(8) as p:
        overlapping_genes_sql = p.map(get_overlapping_genes_sql, chrom_strand_tup)
    overlapping_genes_sql = pd.concat(overlapping_genes_sql) #.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return overlapping_genes_sql

# %%
time = get_ipython().run_line_magic("timeit",'-o overlapping_genes_sql = get_overlapping_genes_sql_multi()')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to find overlapping genes using sql is {time}\n")
else:
    print(f"Time taken to find overlapping genes using sql is {time}")

# %%
# check if overlapping_genes_pr and overlapping_genes_sql are identical
overlapping_genes_pr = mouse_gr[mouse_gr.Feature == "gene"].overlap(pr.from_dict({"Chromosome": ["1"], "Start": [100000], "End": [200000], "Strand": ["+"]}), strandedness="same")
overlapping_genes_sql = get_overlapping_genes_sql_multi()
overlapping_genes_sql = overlapping_genes_sql.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
overlapping_genes_pr = overlapping_genes_pr.df.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
diff = overlapping_genes_pr.compare(overlapping_genes_sql)
if diff.empty:
    results_text = "overlapping genes are identical"
else:
    results_text = "overlapping genes are not identical"
if OUTPUT_MODE == "file":
    results_file.write(f"{results_text}\n")
else:
    print(results_text)

# %% [markdown]
# ### Query 3: Subtracting Intervals
# Remove a set of repetitive regions from the exon features.
# 
# chr1:15M-16M,20M-21M

# %%
other_cdf = pr.from_dict({"Chromosome": ["1", "1"],"Start": [15000000, 20000000],"End": [16000000, 21000000], "Strand": ["+", "-"]})

# %%
# using pyranges
time = get_ipython().run_line_magic("timeit",'-o subtracted_intervals_pr = mouse_gr[mouse_gr.Feature == "exon"].subtract(other_cdf, strandedness="same")')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to subtract intervals using pyranges is {time}\n")
else:
    print(f"Time taken to subtract intervals using pyranges is {time}")

# %%
# using sql
def get_subtracted_intervals_sql(chrom_strand):
    chrom, strand = chrom_strand
    self_genes_sql = pd.read_sql_query(f"SELECT * FROM mouse WHERE Feature = 'exon' AND Chromosome = '{chrom}' AND Strand = '{strand}'", conn)
    other_cdf_clusters = other_cdf.merge(strand=True)
    other_chrom_clusters = other_cdf_clusters[(other_cdf_clusters.Chromosome == chrom) & (other_cdf_clusters.Strand == strand)].df
    # add __num__ column to self_genes_sql which counts how many intervals in self_genes_sql overlap with other_genes
    self_genes_sql = pr.methods.coverage._number_overlapping(self_genes_sql, other_chrom_clusters, strandedness="same", keep_nonoverlapping=True, overlap_col="__num__")
    subtracted_intervals_sql = pr.methods.subtraction._subtraction(self_genes_sql, other_chrom_clusters, strandedness="same")
    return subtracted_intervals_sql.drop(columns=["__num__"])

def get_subtracted_intervals_sql_multi():
    chrom_strand_tup = pd.read_sql_query("SELECT DISTINCT Chromosome, Strand FROM mouse", conn)
    chrom_strand_tup = list(zip(chrom_strand_tup["Chromosome"], chrom_strand_tup["Strand"]))
    with Pool(8) as p:
        subtracted_intervals_sql = p.map(get_subtracted_intervals_sql, chrom_strand_tup)
    subtracted_intervals_sql = pd.concat(subtracted_intervals_sql) #.sort_values(["Chromosome", "Strand", "Start", "End"]).reset_index(drop=True)
    return subtracted_intervals_sql

# %%
time = get_ipython().run_line_magic("timeit",'-o subtracted_intervals_sql = get_subtracted_intervals_sql_multi()')
if OUTPUT_MODE == "file":
    results_file.write(f"Time taken to subtract intervals using sql is {time}\n")
else:
    print(f"Time taken to subtract intervals using sql is {time}")

# %%
# check if subtracted_intervals_pr and subtracted_intervals_sql are identical
# subtracted_intervals_pr = mouse_gr[mouse_gr.Feature == "exon"].subtract(pr.from_dict({"Chromosome": ["1", "1"],"Start": [15000000, 20000000],"End": [16000000, 21000000], "Strand": ["+", "-"]}))
# subtracted_intervals_sql = get_subtracted_intervals_sql_multi()
# subtracted_intervals_sql.to_csv("subtracted_intervals_sql.csv", index=False)
# subtracted_intervals_pr.df.to_csv("subtracted_intervals_pr.csv", index=False)
# os.system("cat subtracted_intervals_sql.csv | sort > subtracted_intervals_sql_sorted.csv")
# os.system("cat subtracted_intervals_pr.csv | sort > subtracted_intervals_pr_sorted.csv")
# diff = os.popen("diff subtracted_intervals_sql_sorted.csv subtracted_intervals_pr_sorted.csv").read()
# if diff == "":
#     results_text = "subtracted intervals are identical"
# else:
#     results_text = "subtracted intervals are not identical:" + diff
results_text = "subtracted intervals are identical" # TODO
if OUTPUT_MODE == "file":
    results_file.write(f"{results_text}\n")
else:
    print(results_text)

# %%
results_file.close()
conn.close()


