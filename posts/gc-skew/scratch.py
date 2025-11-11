import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")

with app.setup:
    # Initialization code that runs before all other cells
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['figure.dpi'] = 400
    plt.rcParams['figure.figsize'] = (10, 6)
    plt.rcParams['figure.constrained_layout.use'] = True


@app.cell
def _():
    with open("data/Vibrio_cholerae.txt", "r") as fh:
        genome = fh.read()

    L = len(genome)
    positions = range(0, L+1)
    skew = np.zeros(L+1, dtype=int)

    for pos, base in enumerate(genome, start=1):
        if base == "G":
            skew[pos] = skew[pos-1] + 1
        elif base == "C":
            skew[pos] = skew[pos-1] - 1
        else:
            skew[pos] = skew[pos-1]
    return genome, positions, skew


@app.cell
def _(positions, skew):
    def format_position(x, pos): return f"{x/1_000} Kb"

    fig, ax = plt.subplots(layout="tight")
    ax.plot(positions, skew)
    ax.axvspan(145_000, 155_000, color='lightgreen', alpha=0.3, label='ori')
    ax.axvspan(565_000, 575_000, color='lightcoral', alpha=0.3, label='ter')

    ax.xaxis.set_major_formatter(format_position)
    ax.set_xlabel("Position")
    ax.set_ylabel("Skew \n($count_G - count_C$)")

    plt.show()
    return


app._unparsable_cell(
    r"""
    def kmer_counts(genome: str, k: int):
    

    def most_frequent_kmer(genome: str, k: int):
        L = len(genome)
        kmer_counts = {}
        for i in range(L-k+1):
            kmer = genome[i:i+k]
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1

        max_count = max(kmer_counts.values())
        max_kmers = []
        for kmer, count in kmer_counts.items():
            if count == max_count:
                max_kmers.append(kmer)

        return max_kmers
    """,
    name="_"
)


@app.cell
def _(genome, most_frequent_kmer):
    ks = range(3, 20)
    n_kmers = [len(most_frequent_kmer(genome, k)) for k in ks]
    n_kmers = np.array(n_kmers, dtype=np.int8)
    return ks, n_kmers


@app.cell
def _(ks, n_kmers):
    plt.bar(ks, n_kmers)
    plt.xlim(2, 21)
    plt.show()
    return


@app.cell
def _(genome, most_frequent_kmer):
    mer11 = most_frequent_kmer(genome, 11)
    mer11
    return


@app.cell
def _(genome):
    with open("data/vcholerae.fasta", "w") as out:
        out.write(">Vibrio cholerae genome\n")
        for i in range(0, len(genome), 60):
            out.write(genome[i:i+60])
            out.write("\n")
    return


@app.cell
def _(genome):
    genome[196758:197146].find("TTTCCAGT")
    return


@app.cell
def _(genome, most_frequent_kmer):
    most_frequent_kmer(genome, 3)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
