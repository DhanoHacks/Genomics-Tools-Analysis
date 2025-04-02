import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import numpy as np

def plot_bars(ax, labels, pyranges_times, sql_times, title):
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width/2, pyranges_times, width, label='PyRanges', color='b')
    ax.bar(x + width/2, sql_times, width, label='SQL', color='g')
    ax.set_ylabel('Time (s)')
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20)
    ax.legend()

def plot_stacked_bar(ax, pyranges_time, sql_write_time, sql_read_time, title):
    ax.bar(1, pyranges_time, label='GTF to Pyranges', color='b')
    ax.bar(2, sql_write_time, label='GTF to SQL', color='r')
    ax.bar(2, sql_read_time, bottom=sql_write_time, label='SQL to Pyranges', color='g')
    ax.set_ylabel('Time (s)')
    # remove x labels
    ax.set_xticks([])
    ax.set_title(title)
    ax.legend()

def main():
    # Human data
    pyranges_read_h = 63.963
    sql_write_h = 98.48
    sql_read_h = 30.44
    exons_pyranges_h = 4.56
    exons_sql_h = 1.45
    exon_length_pyranges_h = 5.27
    exon_length_sql_h = 1.64
    transcripts_pyranges_h = 1.01
    transcripts_sql_h = 0.699
    
    # Mouse data
    pyranges_read_m = 67.24
    sql_write_m = 47.90
    sql_read_m = 37.65
    exons_pyranges_m = 6.52
    exons_sql_m = 1.6
    exon_length_pyranges_m = 7.39
    exon_length_sql_m = 1.82
    transcripts_pyranges_m = 1.59
    transcripts_sql_m = 0.707
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    
    # Stacked bar for human
    plot_stacked_bar(axs[0], pyranges_read_h, sql_write_h, sql_read_h, 'Read Times (Human Genome)')
    
    # Comparisons for human
    plot_bars(axs[1], ['Count Exons', 'Total Exon Length', 'Most Transcripts'],
              [exons_pyranges_h, exon_length_pyranges_h, transcripts_pyranges_h],
              [exons_sql_h, exon_length_sql_h, transcripts_sql_h],
              'Aggregate Queries (Human Genome)')
    
    plt.tight_layout()
    plt.savefig("results_human.png")
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    
    # Stacked bar for mouse
    plot_stacked_bar(axs[0], pyranges_read_m, sql_write_m, sql_read_m, 'Read Times (Mouse Genome)')
    
    # Comparisons for mouse
    plot_bars(axs[1], ['Count Exons', 'Total Exon Length', 'Most Transcripts'],
              [exons_pyranges_m, exon_length_pyranges_m, transcripts_pyranges_m],
              [exons_sql_m, exon_length_sql_m, transcripts_sql_m],
              'Aggregate Queries (Mouse Genome)')
    
    plt.tight_layout()
    plt.savefig("results_mouse.png")

if __name__ == "__main__":
    main()
