import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import numpy as np
import re

def parse_time_with_units(s):
    """Parses a string like '1.2 s' or '685 ms' into seconds (float)."""
    value, unit = re.match(r'([\d.]+)\s*(ms|s)', s).groups()
    value = float(value)
    return value / 1000 if unit == 'ms' else value

def plot_bars(ax, labels, pyranges_times, pyranges_errs, sql_times, sql_errs, title):
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width/2, pyranges_times, width, yerr=pyranges_errs, label='PyRanges', color='b', capsize=5)
    ax.bar(x + width/2, sql_times, width, yerr=sql_errs, label='SQL', color='g', capsize=5)
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
    ax.set_xticks([])
    ax.set_title(title)
    ax.legend()

def extract_float(pattern, text):
    return float(re.findall(pattern, text)[0])

def extract_with_units(pattern, text):
    raw = re.findall(pattern, text)[0]
    return parse_time_with_units(raw)

def parse_results_file(path):
    text = open(path).read()

    results = {}
    results['pyranges_read'] = extract_float(r'Finished reading into pyranges object in time ([\d.]+)', text)
    results['sql_write'] = extract_float(r'Finished writing into sql db in time ([\d.]+)', text)
    results['sql_read'] = extract_float(r'Finished reading from sql to pyranges in time ([\d.]+)', text)

    def get_pair(query, label):
        base = fr'Time taken to {query} using {label} is ([\d.]+|\d+) (ms|s) ± ([\d.]+|\d+) (ms|s)'
        m = re.search(base, text)
        if not m:
            raise ValueError(f"Could not find timing for {query} using {label}")
        time = float(m.group(1)) / 1000 if m.group(2) == 'ms' else float(m.group(1))
        err = float(m.group(3)) / 1000 if m.group(4) == 'ms' else float(m.group(3))
        return time, err

    for query_key, query_phrase in {
        'exons': 'count exons',
        'exon_length': 'calculate total exon length',
        'transcripts': 'find chromosome with most transcripts',
        'merge': 'merge exons',
        'overlap': 'find overlapping genes',
        'subtract': 'subtract intervals'
    }.items():
        results[f'{query_key}_pyr'], results[f'{query_key}_pyr_err'] = get_pair(query_phrase, 'pyranges')
        results[f'{query_key}_sql'], results[f'{query_key}_sql_err'] = get_pair(query_phrase, 'sql')

    return results

def make_plot(results, species_name, filename):
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    plot_stacked_bar(axs[0], results['pyranges_read'], results['sql_write'], results['sql_read'], f'Read Times ({species_name})')

    plot_bars(
        axs[1],
        ['Count Exons', 'Total Exon Length', 'Most Transcripts'],
        [results['exons_pyr'], results['exon_length_pyr'], results['transcripts_pyr']],
        [results['exons_pyr_err'], results['exon_length_pyr_err'], results['transcripts_pyr_err']],
        [results['exons_sql'], results['exon_length_sql'], results['transcripts_sql']],
        [results['exons_sql_err'], results['exon_length_sql_err'], results['transcripts_sql_err']],
        f'Aggregate Queries ({species_name})'
    )

    plot_bars(
        axs[2],
        ['Merge Exons', 'Find Overlaps', 'Subtract Intervals'],
        [results['merge_pyr'], results['overlap_pyr'], results['subtract_pyr']],
        [results['merge_pyr_err'], results['overlap_pyr_err'], results['subtract_pyr_err']],
        [results['merge_sql'], results['overlap_sql'], results['subtract_sql']],
        [results['merge_sql_err'], results['overlap_sql_err'], results['subtract_sql_err']],
        f'Interval Queries ({species_name})'
    )

    plt.tight_layout()
    plt.savefig(filename)

def main():
    human_results = parse_results_file("results_human.txt")
    mouse_results = parse_results_file("results_mouse.txt")
    make_plot(human_results, "Human Genome", "results_human.png")
    make_plot(mouse_results, "Mouse Genome", "results_mouse.png")

if __name__ == "__main__":
    main()
