#import "@preview/oxdraw:0.1.0": *
#import "bench-table.typ": bench-table, duplication-rate-table
#import "comparison-tables.typ": build-comparaison-table, build-ram-comparaison-table
#import "@preview/cetz:0.4.2": canvas, draw
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge, shapes
#import "@preview/fontawesome:0.5.0": *
#import "@preview/codly:1.3.0": *
#import "@preview/codly-languages:0.1.1": *
#show: codly-init.with()

#codly(languages: codly-languages)


#set document(
  title: "FastDedup: A fast and memory-efficient tool for read deduplication",
  date: auto,
)

#set page(
  paper: "a4",
  margin: (x: 2cm, y: 2cm),
  header: context {
    if counter(page).get().first() > 1 [
      #text(size: 9pt, fill: gray)[FastDedup -- A fast and memory-efficient tool for read deduplication #h(1fr) R. Ribes, C. Mandier & A. Baniel]
      #line(length: 100%, stroke: 0.5pt + gray)
    ]
  },
  footer: context {
    let page-num = counter(page).get().first()
    align(center)[#text(size: 9pt)[#page-num]]
  },
)

#set text(size: 11pt, lang: "en")
#set heading(numbering: "1.1")
#set par(justify: true)

// Style des liens
#show link: it => text(fill: blue, it)

// =============================================================================
// Title page
// =============================================================================

#align(center)[
  #v(3cm)
  #text(size: 24pt, weight: "bold")[FastDedup]
  #v(0.5cm)
  #text(size: 18pt, weight: "bold")[A fast and memory-efficient tool for read deduplication]
  #v(3cm)
  #text(size: 12pt)[
    *Raphaël Ribes*#super[1], *Céline Mandier*#super[1], *Alice Baniel*#super[2] \
    #v(0.5em)
    #text(size: 11pt)[
      1: ISDM, Univ Montpellier, CIRAD, INSERM, Montpellier, France \
      2: ISEM, Univ Montpellier, CNRS, IRD, Montpellier, France
    ]
  ]
  
  #v(1cm)
  #text(size: 11pt, fill: gray)[#datetime.today().display("[day] [month repr:long] [year]")]
]

#v(1fr)

#rect(fill: luma(245), inset: 1em, width: 100%)[
  *Abstract* ---
  PCR duplicate removal is a critical first step in high-throughput
  sequencing pipelines, yet existing tools struggle with speed, memory, or
  correctness at modern dataset scales.
  We present `FastDedup`, a Rust-based FASTX deduplicator that transforms each
  read or read pair to a compact `xxh3` hash fingerprint, drastically 
  reducing memory usage and binding most of the execution time to disk I/O.
  Benchmarked against six competing tools on synthetic human WGS datasets
  up to 300 million reads, `FastDedup` consistently leads on paired-end
  data, running more than 10 times faster than `fastp`.
  It also outperforms all tools on uncompressed single-end
  data, deduplicating a million reads in a second.
  We additionally report correctness failures in `prinseq++` and `clumpify`.
  `FastDedup` is available under the MIT License via GitHub, Bioconda, and Cargo.
]

#pagebreak()
// =============================================================================
// Table of contents
// =============================================================================

#outline(
  title: [Table of contents],
  indent: auto,
  depth: 2,
)

#pagebreak()

// =============================================================================
// Context
// =============================================================================
= -- Context

#v(1em) 
With the rise of Next-Generation Sequencing (NGS) technologies, the volume of sequencing data has increased exponentially.
With this surge in data, the development of new tools for processing and analyzing sequencing reads became necessary.
One of the first and most critical steps in the analysis pipeline is read deduplication, which involves identifying and removing duplicate sequences that can arise from PCR amplification or sequencing errors.
This is crucial when working with DNA-seq data, but it can introduce bias in RNA-seq data since the duplication could be biological in nature (e.g., highly expressed genes).

#v(1em)

#align(center)[
  #figure(
    canvas({
      import draw: *

      let draw_read(x, y, h, label_a, label_b) = {
        line((x, y), (x, y + h), stroke: 1.5pt + black)

        // R1
        line((x - 0.2, y), (x - 0.2, y + 1.1), stroke: 1pt + red, mark: (end: ">", fill: red))
        content((x - 0.5, y + 0.3), text(fill: red, size: 8pt)[#label_a])

        // R2
        line((x + 0.2, y + h), (x + 0.2, y + h - 1.1), stroke: 1pt + blue, mark: (end: ">", fill: blue))
        content((x + 0.5, y + h - 0.3), text(fill: blue, size: 8pt)[#label_b])
      }

      // Not a duplication
      content((0.75, 3.5), text(weight: "bold", size: 1.2em)[A)])
      draw_read(0, 0, 3, $"R1"_1$, $"R2"_1$)
      draw_read(1.5, 0, 2, $"R1"_2$, $"R2"_2$)
      
      content((0.75, -0.5), text(fill: red, size: 9pt)[$"R1"_1 = "R1"_2$])
      content((0.75, -1.0), text(fill: blue, size: 9pt)[$"R2"_1 != "R2"_2$])

      // Duplication
      content((5.75, 3.5), text(weight: "bold", size: 1.2em)[B)])
      draw_read(5, 0, 3, $"R1"_1$, $"R2"_1$)
      draw_read(6.5, 0, 3, $"R1"_2$, $"R2"_2$)

      content((5.75, -0.5), text(fill: red, size: 9pt)[$"R1"_1 = "R1"_2$])
      content((5.75, -1.0), text(fill: blue, size: 9pt)[$"R2"_1 = "R2"_2$])
    }),
    caption: [Identifying PCR Duplicates in Paired-End Sequencing.],
  )
]

#v(1em)

#align(center)[
  #block(width: 100%, stroke: .5pt + luma(150), radius: 2pt, inset: 10pt, align(left)[
    #v(0.5em)
    *Solid Black Line:* Represents a single sequenced DNA fragment (height reflects length).
    
    *Red Arrow (#text(fill: red)[R1]):* Forward read (Read 1).
    
    *Blue Arrow (#text(fill: blue)[R2]):* Reverse read (Read 2).
    #v(0.5em)
    *A) Not a duplicate:* Shows two distinct fragments where the start position of Read 1 is shared, but different fragment lengths lead to unique mapping positions for Read 2. \
    #v(0.5em)
    *B) PCR Duplicate:* Shows two exact fragments where the start position of Read 1 and Read 2 is shared.
  ])
]

Given a fragment $F$, we define the 5'->3' sequence as $R 1$ with a length of $r_1$, and the 3'->5' sequence as $R 2$ with a length of $r_2$.
$R 1$ is the prefix of $F$, the reverse complement of $R 2$ is the suffix of $F$, and $r_1 = r_2$.

Therefore, a PCR duplicate occurs when $F_1$ and $F_2$ share the same starting and ending positions in the genome, and thus have the same length.
In that case, $R 1_1$ and $R 1_2$ are identical, and $R 2_1$ and $R 2_2$ are identical.

This step is crucial because removing these artifacts prevents false positives during variant calling, as PCR polymerases can introduce errors during amplification cycles.
If an error occurs early, it is propagated into multiple daughter molecules.
When sequenced, these identical, error-carrying molecules will appear as multiple independent reads supporting a mutation that does not actually exist in the biological sample @Ebbert2016.

Deduplication also helps correct allele frequency bias, as not all DNA molecules amplify with equal efficiency due to factors such as GC content or fragment length.
This uneven amplification can artificially inflate the coverage of certain alleles.
It also improves _de novo_ assembly, especially in whole-genome sequencing.
When paired-end reads are mapped to estimate the order and intervening distance between contigs, duplicates can create false-positive connections between contigs or introduce conflicting connections @Zhang2023.

Before presenting our tool, we will first review existing tools for read deduplication and their limitations.

#pagebreak()
= -- Existing tools and limitations 
== : Fastp
#v(1em)

`fastp` is a widely adopted, all-in-one FASTQ preprocessor written in C++ that performs quality control, adapter trimming, and deduplication within a single pass @Chen2018 @Chen2023.
Unlike alignment-based tools that require computationally expensive mapping and coordinate sorting, `fastp` operates directly on raw sequences.

To achieve high performance, `fastp` implements a probabilistic deduplication mechanism based on Bloom filters and multiple hash arrays.
For each read, or read pair, the sequence is mapped to multiple integer hash values.
If all corresponding positions in the hash arrays are already marked as positive, the read is flagged as a duplicate.
This streaming, single-pass architecture avoids reading and writing intermediate files to disk, drastically reducing input/output (I/O) bottlenecks and resulting in ultra-fast execution speeds.

A significant advantage of `fastp` is its highly configurable memory management.
Users can adjust the `--dup_calc_accuracy` parameter to scale the hash buffer size from 1 GB up to 24 GB, allowing the tool to run on both standard personal computers and High-performance Computing (HPC) clusters.

However, this probabilistic approach involves a minor trade-off in absolute precision, as `fastp` has a theoretical false-positive rate.
Furthermore, `fastp`'s deduplication algorithm relies heavily on the exact matching of sequence lengths and coordinate regions.
Consequently, applying certain quality trimming operations (such as sliding-window trimming via `--cut_front` or `--cut_tail`) prior to deduplication can interfere with the duplicate identification process, as trimmed reads may no longer perfectly match their duplicate counterparts.

We will be using the 1.1.0 version in this study.

== : FastUniq
#v(1em)

`FastUniq` @Xu2012 is a deterministic deduplication tool explicitly designed for paired-end (PE) short reads.
By directly comparing the nucleotide sequences of read pairs, it provides an exact identification of PCR duplicates.

Algorithmically, `FastUniq` operates by loading the entirety of the paired reads into memory, sorting them based on their sequences using a merge sort algorithm, and identifying duplicates by sequentially comparing adjacent read pairs in the sorted list.\
A notable feature of its comparison logic is its ability to handle reads of varying lengths: a shorter sequence is considered identical to a longer one if it perfectly matches its 5' prefix.
If duplicates are found, the pair containing the longest reads is retained.

While `FastUniq` guarantees absolute precision for exact duplicate removal, its architectural design introduces critical resource bottlenecks.
Because it relies on importing the entire dataset into RAM simultaneously using a hierarchical storage structure, its memory footprint scales linearly with the input file size.
Benchmarking studies have highlighted this as a severe limitation for modern high-throughput data.
For instance, processing a moderately sized dataset of approximately 69 million read pairs required over 50 GB of RAM @Dai2020.

While highly accurate, `FastUniq` is often impractical for massive datasets without access to HPC clusters.
Furthermore, its scope is strictly limited to paired-end data, and it does not tolerate sequence mismatches.

We will be using the 1.1 version in this study.

== : Clumpify
#v(1em)

`clumpify`, a tool from the Java-based BBTools suite @Bushnell2014, introduces an algorithmic paradigm to duplicate removal based on $k$-mer grouping.
Rather than relying on exact sequence hashing or global coordinate sorting, `clumpify` identifies reads that share common $k$-mers and groups them into localized clusters or "clumps". 

This spatial reorganization provides a unique secondary benefit. 
By physically placing highly similar sequences adjacent to one another in the output file, `clumpify` dramatically improves standard `gzip` compression efficiency, as the compression algorithm can utilize much closer reference pointers. 

A major functional advantage of `clumpify` is its high sensitivity to optical duplicates.
By parsing the flowcell tile and pixel coordinates embedded within Illumina FASTQ headers, `clumpify` (via parameters such as `dupedist`) can explicitly distinguish between amplification-derived PCR duplicates and sensor-derived optical duplicates.
For paired-end data, `clumpify` predominantly uses the first read (R1) to determine the clump, allowing the paired read (R2) to passively follow its counterpart, which maintains file synchronization while minimizing computational overhead.

To manage memory, `clumpify` implements a multi-phase strategy (KmerSplit and KmerSort).
If the dataset exceeds available RAM, the tool can divide the data into temporary groups on the disk, sorting them individually before merging.

However, its reliance on the Java Virtual Machine (JVM) can lead to instability and high memory overhead on massive datasets.
By default, high-memory tools within the BBTools suite will attempt to allocate all available system RAM @Bushnell2014, and JVM garbage collection overhead can cause severe bottlenecks unless strict multi-threading and heap-size limits are manually tuned.

We will be using the 39.79 bbmap version in this study.

== : SeqKit
#v(1em)

`SeqKit` @Shen2016 @Shen2024 is a comprehensive, cross-platform, and ultra-fast toolkit for FASTA and FASTQ file manipulation.
Written in Go, it benefits from efficient concurrency management and the ability to produce static binaries with no external dependencies, making it highly accessible and stable across different operating systems.

For PCR deduplication, `SeqKit` provides the `rmdup` subcommand.
This tool relies on an exact sequence matching approach using a hashing strategy.
To process massive datasets without exhausting system memory, `seqkit rmdup` offers an MD5 mode (invoked via the `-m` or `--md5` flag).
In this mode, the algorithm computes and stores only the MD5 digest of each sequence in the RAM rather than the full nucleotide sequence.
This extreme memory efficiency allows `SeqKit` to maintain a very low resource footprint, even outperforming newer tools (e.g `prinseq++`) when handling vast amounts of single-end data.

Despite its speed and memory efficiency, `SeqKit rmdup` has a notable architectural limitation: it does not natively support the synchronized deduplication of paired-end reads.
The `rmdup` command evaluates and processes individual files independently.
It checks for duplicates inside each file separately; therefore, if $R 1_1=R 1_2$ but $R 2_1 eq.not R 2_2$, $R 1_2$ will be flagged as a duplicate and removed while $R 2_2$ is retained, leading to data loss and desynchronization of the paired files.
While workarounds exist using other `SeqKit` subcommands (such as `seqkit pair` to resynchronize files), `seqkit rmdup` is generally considered the optimal choice strictly for single-end libraries.
This is why we will not use it for paired-end benchmarking in this study.

We will be using the 2.13.0 version in this study.

== : CD-HIT-DUP
#v(1em)

`CD-HIT-DUP`, a specialized tool within the broader CD-HIT package @Fu2012 @Huang2010, is designed to identify and remove both identical and nearly-identical duplicates.
Unlike tools restricted to exact string matching, `CD-HIT-DUP` is particularly robust against sequencing errors.
It incorporates a greedy incremental clustering algorithm that tolerates a defined threshold of mismatches, including both substitutions and insertions/deletions.
`CD-HIT-DUP` relies on a prefix-suffix comparison approach. The tool first bins reads that share an identical prefix of $k$ nucleotides.
To accelerate this binning and sorting process, it uses a dual numeric encoding strategy, allowing prefixes of up to 27 nucleotides to be compactly represented as 64-bit integers.
Within each bin, sequences are sorted by decreasing length, with the longest sequence acting as the cluster's seed.
The algorithm then evaluates the suffixes of the remaining sequences against this seed, utilizing a $k$-mer word-based index table to statistically predict sequence identity and bypass computationally expensive full pairwise alignments when possible.
For paired-end libraries, `CD-HIT-DUP` conceptually joins the read pairs into a single sequence and checks prefixes at both ends simultaneously to ensure paired-end integrity.

We will be using the 4.8.1 cd-hit-auxtools version in this study.

== : PRINSEQ++
#v(1em)

`prinseq++` @Cantu2019 is an optimized C++ implementation of the Perl-based `prinseq-lite` program.
It is designed to perform a comprehensive suite of quality control, filtering, trimming, and reformatting tasks for genomic and metagenomic sequence data.

For sequence deduplication, `prinseq++` departs from exact string comparison or sorting heuristics by implementing a probabilistic data structure known as a Bloom filter.
During execution, each sequence is transformed using multiple fast, non-cryptographic hash functions, and the corresponding indices in a shared bit-array are flagged.
If a newly evaluated sequence maps entirely to bits that are already set to positive, it is classified as a duplicate.
This probabilistic model is highly advantageous for multi-threaded environments because it enables asynchronous read and write access to the filter with minimal thread blocking.\
Additionally, `prinseq++` natively supports the direct reading and writing of `gzip`-compressed FASTQ and FASTA files without the need for intermediate disk I/O, which significantly reduces hard drive storage requirements.

It also natively supports paired-end reads by processing pairs simultaneously and outputting synchronized R1 and R2 files.
However, the deduplication function (`-derep`) is strictly limited to exact duplicates and does not tolerate sequencing errors.
Furthermore, the probabilistic nature of the Bloom filter introduces a theoretically non-zero rate of false positives.
Finally, when `prinseq++` is executed with multiple threads, the output sequences are not guaranteed to remain in the same order as the input files, which may require consideration if downstream tools expect strict global index synchronization.

We will be using the 1.2.4 version in this study.

== : SAM/BAM ecosystem tools
#v(1em)

In addition to tools that operate directly on raw FASTQ sequences, the bioinformatics ecosystem relies heavily on utilities designed for SAM and BAM formats.
These include PCR deduplication tools such as `samtools` (via `markdup`) @Li2009 and Picard's `MarkDuplicates` @Picard.\
While these tools are typically used downstream after mapping reads to a reference genome, it is technically possible to convert raw, unaligned FASTQ files directly into unmapped SAM/BAM files to process them through these pipelines.
However, even without the computationally heavy step of genomic alignment, the sheer format conversion from FASTQ to SAM/BAM introduces significant and unnecessary I/O and processing overhead.
Because our primary objective is to evaluate and develop a highly optimized, memory-efficient tool strictly for direct FASTQ-to-FASTQ deduplication, these format-dependent tools fall outside the scope of our comparative analysis.

#pagebreak()
= -- FastDedup: High-Performance FASTX Deduplication
#v(1em)

The goal behind the development of `FastDedup` (or `FDedup`; @SE_arch @PE_arch) was to create an exact FASTX PCR deduplication tool that prioritizes maximum speed and memory efficiency.
To achieve this, the tool relies on `xxh3` @Collet_xxHash @DoumanAsh_xxhash_rust, a non-cryptographic hash function, to compute a unique fingerprint for each read.
These fingerprints are cached in memory using `fxhash` @cbreeden_fxhash, which provides a low memory footprint.

`FastDedup` provides users with precise control over the hash collision rate while maintaining high performance.
It automatically scales its hashing strategy, seamlessly switching between 64-bit and 128-bit hashes based on the estimated input sequence count and a user-defined collision probability threshold.

Specifically designed for seamless integration into OMICS pipelines, the tool features robust incremental deduplication and auto-recovery mechanisms.
In the event of an interruption, `FastDedup` can safely preload existing hashes from the output file to prevent duplicate processing and smoothly resume operations.
If an uncompressed output file becomes corrupted due to a crash, the tool automatically detects the issue, calculates a fail-safe truncation point, and truncates the file to the last valid sequence before continuing.

Finally, a dry-run mode (`--dry-run` or `-s`) is available to calculate the duplication rate without writing any output files.
Since file I/O operations account for the majority of execution time, this mode is exceptionally fast and serves as a highly efficient feature for pipeline planning and data analyses.

== : Algorithmic Complexity and Hardware-Level Optimizations<AlgComplexity>

The architectural principle of `FastDedup`'s performance lies in shifting the deduplication problem from a string-matching paradigm to an integer-matching paradigm.
Transforming raw paired-end sequences (e.g., $2 times 150$ bytes) into a single 64-bit hash reduces the data by 37.5 times in both space and time complexity.

=== Space Complexity and Memory Footprint

In a deterministic approach, storing a read pair requires retaining sequences of length $L$ (where $L = r_1 + r_2$).
The space complexity for $N$ unique sequences scales linearly as $O(N times L)$.
For a standard 150 bp paired-end dataset, each pair consumes at least 300 bytes of raw character data, plus substantial data structure overhead (e.g., string pointers).
Storing $10^8$ sequences typically exceeds 32 GB of RAM.
By computing a 64-bit hash using `xxh3`, `FastDedup` reduces the sequence identity into exactly 8 bytes.
This decouples the memory requirement from the read length, reducing the space complexity to $O(N)$ (since $8 << L$).

=== Time Complexity and Arithmetic Logic Unit (ALU) Efficiency

String-based deduplication necessitates comparing sequences letter by letter, meaning the time it takes scales with the read length.
In the worst-case scenario, comparing two sequences requires $O(L)$ time.
`FastDedup` changes this dynamic. While it still takes $O(L)$ time to read a sequence and calculate its initial hash, every subsequent comparison uses only that hash.
Because the hash is a simple 64-bit integer, the processor's Arithmetic Logic Unit can compare two of them natively in a single clock cycle, instantly dropping the lookup and comparison time complexity to $O(1)$.

=== Cache Locality and Zero-Allocation Pipeline

Beyond theoretical Big-$O$ improvements, the 8-byte hash architecture heavily exploits modern CPU hardware design.
Processors fetch data from RAM in discrete 64-byte cache lines.
A 300-byte sequence spans multiple cache lines, leading to frequent cache misses and severe memory bandwidth bottlenecks.
In contrast, exactly eight 64-bit hashes fit perfectly into a single cache line.
As `fxhash` probes the internal hash set for collisions, it operates almost exclusively within the ultra-fast CPU cache (L1/L2), minimizing latency and maximizing throughput.

Furthermore, `FastDedup` processes byte slices directly from the `needletail` @Needletail parser into the `xxh3` hasher without ever allocating the raw string to the heap.
This zero-allocation pipeline completely bypasses system-level memory allocation (`malloc`) overheads, ensuring that the execution time is strictly bound by disk I/O rather than RAM latency.

== : The Mathematics of Dynamic Hashing
#v(1em)

While `FastDedup` defaults to the highly efficient 64-bit hashing strategy of `xxh3`, it dynamically assesses the risk of hash collisions to guarantee data integrity.
The theoretical risk of a hash collision probability $p$ for $x$ sequences follows the standard Birthday Problem approximation @feller1968probability @cormen2009algorithms:

$ p = x^2 / (2^(N+1)) $

Where $N in {64, 128}$ represents the bit-length of the hash.
This probability expresses the chance that *at least one collision exists anywhere in the entire dataset*, not a per-read false-duplicate rate.
When a collision does occur, it affects exactly one read; since $p$ is small, the expected number of false duplicates is also approximately $p$.

At the default 64-bit configuration and a probability threshold of $0.01$, it requires approximately $0.6 times 10^9$ sequences before the probability of any single collision occurring reaches $1%$.
In practical terms, even at this boundary, a run would produce fewer than one spurious duplicate on average.
If the estimated dataset size pushes the collision probability beyond this threshold, `FastDedup` automatically scales to 128-bit hashing, extending the same guarantee to approximately $2.6 times 10^18$ sequences --- far beyond any current sequencing experiment.

In everyday use cases, the 64-bit hashing strategy is more than sufficient.
The 128-bit option serves as a safeguard for extreme cases, ensuring that `FastDedup` can maintain its performance and integrity even as sequencing datasets continue to grow exponentially.
Even at the theoretical boundary, the expected number of false duplicates remains below one, and any dataset large enough to push beyond that boundary would contain hundreds of millions of reads, making a handful of removals statistically insignificant relative to the overall data volume.

For paired-end reads, maintaining this theoretical probability requires careful bitwise operations.
Let $S_1$ and $S_2$ represent the nucleotide sequences of $R 1$ and $R 2$, and $H$ be our non-cryptographic uniform hash function mapping to ${0,1}^N$.
We define the independent hashes for each read as:

$ h_1 = H(S_1) $
$ h_2 = H(S_2) $

If the algorithm combined these hashes using a simple bitwise XOR ($plus.o$) to create a unified fingerprint $h_u$, a critical failure occurs when $R 1$ and $R 2$ are identical sequences ($S_1 = S_2$).
In this case, $h_1 = h_2$, resulting in:

$ h_u = h_1 plus.o h_1 = 0 $

Every paired-end read with identical forward and reverse sequences would map to $0$.
This would cause massive, deterministic hash collisions, meaning the collision probability would no longer be bounded by the theoretical uniform distribution, but rather by the biological frequency of identical sequence pairs in the dataset.

To prevent this symmetrical cancellation, the algorithm introduces a bitwise left-rotation by $N\/2$ bits, denoted as $R(x, k)$ where $x$ is the hash and $k$ is the number of bits to rotate.\
The unified fingerprint is instead calculated as:

$ h_u = h_1 plus.o R(h_2, N/2) $

Under this operation, even if $S_1 = S_2$:

$ h_u = h_1 plus.o R(h_1, N/2) eq.not 0 $

Because $R(h_1, N/2)$ shifts the lower half of the bits to the upper half and vice versa, $h_1 plus.o R(h_1, N/2)$ will almost never equal $0$.
It would only equal $0$ if $h_1$ consisted of two perfectly identical halves, an event with an astronomically low probability of $2^(-N/2)$.

Since a bitwise rotation is a bijective function over the $N$-bit vector space, it perfectly preserves the original entropy of $h_2$.
When two independent, uniformly distributed random variables are combined via XOR, the resulting variable remains uniformly distributed.
Therefore, the unified fingerprint $h_u$ safely maintains a strict uniform distribution, ensuring the overall collision risk remains uncorrupted and strictly bounded by the initial probability formula.

This is the mathematical foundation that allows `FastDedup` to guarantee a collision probability of less than $1%$ for datasets containing up to approximately 600 million paired-end reads (1.2 billion total sequences) when using 64-bit hashes.

== : FxHash: Algorithmic Mechanics and ALU Optimizations
#v(1em)

While `FastDedup` utilizes `xxh3` to compute sequence fingerprints, it relies on `fxhash` @cbreeden_fxhash to manage these fingerprints within its internal hash set.
Concretely, `FastDedup` uses a `FxHashSet`, a hash set data structure that uses `fxhash` as its internal bucket-routing function: the fixed-size fingerprint produced by `xxh3` is re-hashed by `fxhash` to determine its slot in the table.
This two-level design is intentional, `xxh3` provides high-quality, collision-resistant fingerprints of biological sequences, while `fxhash` provides minimal-overhead addressing into the hash table itself.

`fxhash` is a high-performance, non-cryptographic hashing algorithm derived from the internal hasher utilized in the Rust compiler (`rustc`) and Firefox.
Its architectural strength lies not in collision resistance, but rather in raw computational throughput and extremely low latency.

For every word extracted from the input, `fxhash` applies a highly optimized state transition function.
Let $H$ represent the internal state of the hash (initialized to $0$), and $W$ represent an incoming data word:

$ H_{i+1} = (R_L (H_i, 5) xor W_i) times K \pmod{2^N} $

Where:
- $N$ is the machine word size (typically 64 bits).
- $R_L (x, k)$ denotes a bitwise left-rotation of the state $x$ by $k$ bits (`ROTATE = 5`).
- $K_p$ is a predefined $p=N$-bit magic constant (e.g. $K_(64) = "0x517cc1b727220a95"$).

The bitwise left-rotation ensures that previously hashed bits propagate rapidly across the state vector.
Being a bijective function, it preserves the entropy of $H_i$ while preventing localized clustering.
The XOR then introduces the new word $W_i$ in a single clock cycle, and the wrapping multiplication by the large prime-like constant $K$ acts as an avalanche multiplier, diffusing the newly mixed entropy across all 64 bits.
The algorithm explicitly relies on modulo $2^N$ arithmetic, leveraging the CPU's native handling of integer overflow without any computational penalty.

Within `FastDedup` specifically, `fxhash` receives only the fixed-size fingerprints already produced by `xxh3`, either 64-bit or 128-bit values, corresponding to $L in {8, 16}$ bytes.
Because $L$ is bounded to at most 2 machine words, the state transition function executes in exactly 1 or 2 iterations regardless of the input dataset size, making the bucket-routing step a strict $O(1)$ operation.
This stands in contrast to the general case, where a byte-oriented hasher processing an arbitrary sequence of length $L$ would require $O(L)$ state transitions.

The algorithm operates entirely in-place on a single internal state variable and requires no auxiliary data structures, lookup tables, or intermediate buffers, yielding an auxiliary space complexity of $O(1)$.
Combined with the zero-allocation and cache-locality properties described in @AlgComplexity, this ensures the memory footprint remains flat during hash table operations, preventing fragmentation and avoiding OS-level memory management overhead.

#page(flipped: true)[
    == : Architecture
    
    #figure(
      image("SE.mmd.svg", width: 110%),
      caption: [
        Architectural design of `FastDedup` for single-end reads.
      ],
    ) <SE_arch>
    #v(1fr)
    #figure(
      image("PE.mmd.svg", width: 110%),
      caption: [
        Architectural design of `FastDedup` for Paired-End reads.
      ],
    ) <PE_arch>

]

= -- Benchmarking and performance evaluation

To test the performance of the different tools mentioned in the previous sections (@tools_table), we benchmarked them on a series of synthetic datasets of varying sizes, ranging from 100,000 to 300 million paired-end reads, to which we added 10% duplicates.
These reads were generated based on a real human whole genome (GRCh38).
This allows us to evaluate the scalability of each tool and its suitability for different types of sequencing projects, from small-scale experiments to large population genomics studies.
To run these benchmarks in a controlled environment, we set up a Snakemake pipeline where each tool had 1 core and 33 GB of RAM available.
This pipeline ran on the High-Performance Computing (HPC) ISDM-MESO cluster with AMD EPYC 9654 processors.
If a tool exceeds 32 GB of RAM, it is killed and documented as OOM (Out Of Memory) in the results table.
At a time when RAM and electricity are becoming increasingly costly and environmentally impactful, it is crucial to consider the resource efficiency of bioinformatics tools.
Since these tools represent the first step in data preprocessing, it is inefficient to run them using excessive resources.

== : Compressed Paired-End FASTQ

Here we will evaluate the performances of the tools in a Paired-End Compressed in Compressed out (PE-cc) scenario.
Before discussing the results, it should be noted that `FastUniq` and `CD-HIT-DUP` are excluded from this benchmark because they do not support compressed files.
The default compression level of `fastp` is 4, whereas all other tools use a compression level of 6, which is the default for `gzip`.

#figure(
    image("results/pe_cc_runtime_and_ram.svg", width: 100%),
    caption: [
    Execution time and memory usage for Compressed Paired-End FASTQ deduplication.
    ],
)<PE_cc_bench>

After setting `fastp` to use the same compression level, we observe that `FastDedup` is around 50% faster than `fastp` while using significantly less RAM up to 300M reads, at which point it uses roughly the same amount as `fastp` for the 64-bit version (@PE_cc_bench, @PE_cc_table).
We also observed no significant differences in execution time between the 64-bit and 128-bit hashing strategies.
Furthermore, the RAM usage for the 128-bit strategy is strictly doubled compared to the 64-bit strategy, perfectly reflecting the twofold increase in fingerprint size.
Looking closely at the memory scaling, `clumpify`'s RAM usage increases steeply, exhausting the 32 GB limit before reaching 100M reads (consuming over 23 GB at 50M reads alone).
In contrast, `FastDedup` maintains a very flat memory profile, using only around 4.6GB of RAM at 300M reads.
Execution time for `FastDedup` and other tools scales almost linearly with volume, demonstrating their respective I/O-bound bottlenecks, but `FastDedup` remains faster than `fastp` by approximately 50%.

== : Uncompressed Paired-End FASTQ

Here we will evaluate the performances of the tools in a Paired-End Uncompressed in Uncompressed out (PE-dd) scenario.
For uncompressed FASTQ files as inputs and outputs, the raw computational speed of the tools is fully exposed.
`FastDedup` consistently outpaces all other tools in execution time (@PE_dd_bench, @PE_dd_table), completing 300M reads in under 10 minutes (approx. 450s).
In comparison, `prinseq++` requires approximately 2700s and `fastp` roughly 4500s for the same volume---meaning `FastDedup` is around 6× faster than `prinseq++` and 10× faster than `fastp`.
We can also see that `FastDedup` maintains a highly stable memory footprint (around 4.6GB and 8.7GB of RAM for 64-bit and 128-bit hashing respectively at 300M reads), while `prinseq++` edges it out in strict memory efficiency starting at 5M reads.
Deterministic string-matching tools struggle significantly with memory capacity here: `CD-HIT-DUP` and `FastUniq` load full sequences into RAM, causing vertical spikes in their memory usage leading to OOM crashes early on (at 30M and 50M reads respectively). `clumpify` also succumbed to OOM errors at 300M reads.


#figure(
    image("results/pe_dd_runtime_and_ram.svg", width: 100%),
    caption: [
    Execution time and memory usage for Uncompressed Paired-End FASTQ deduplication.
    ],
)<PE_dd_bench>

#v(1fr)

== : Compressed in, uncompressed out Paired-End FASTQ

Here we will evaluate the performances of the tools in a Paired-End Compressed in Uncompressed out (PE-cd) scenario.
#figure(
    image("results/pe_cd_runtime_and_ram.svg", width: 100%),
    caption: [
    Execution time and memory usage for Compressed-input to Uncompressed-output Paired-End FASTQ deduplication.
    ],
)<PE_cd_bench>

We observe very similar results to those obtained with uncompressed files (@PE_cd_bench, @PE_cd_table).
`FastDedup` being the fastest tool and using the least amount of RAM up until 5M reads, where `prinseq++` starts to use less memory than `FastDedup`.
We can also note very similar results for relative performance between each tool, which means that all the tools scale in a similar way whether the input files are compressed or not.

== : Compressed Single-End FASTQ

Here we will evaluate the performances of the tools in a Single-End Compressed in Compressed out (SE-cc) scenario.
#figure(
    image("results/se_cc_runtime_and_ram.svg", width: 100%),
    caption: [
    Execution time and memory usage for Compressed Single-End FASTQ deduplication.
    ],
)<SE_cc_bench>

For compressed single-end data, `SeqKit` and `fastp` hold a slight advantage in execution time (@SE_cc_bench, @SE_cc_table) over `FastDedup` due to highly optimized compression pipelines.
However, pulling attention to the memory graph reveals a dramatic difference in scalability: while `SeqKit`'s memory scales aggressively up to approximately 19 GB at 300M reads, `FastDedup`'s memory footprint strictly remains stable around 4.6 GB (a flat fourfold reduction).
This showcases `FastDedup`'s superior suitability for memory-restricted environments, even when it is roughly 25% slower in pure elapsed time.

== : Uncompressed Single-End FASTQ

Here we will evaluate the performances of the tools in a Single-End Uncompressed in Uncompressed out (SE-dd) scenario.
#figure(
    image("results/se_dd_runtime_and_ram.svg", width: 100%),
    caption: [
    Execution time and memory usage for Uncompressed Single-End FASTQ deduplication.
    ],
)<SE_dd_bench>

Removing the `gzip` bottleneck entirely, `FastDedup` processes 300M single-end reads in roughly 363 seconds, demonstrating its fundamental algorithmic advantage in parsing and deduplicating (@SE_dd_bench, @SE_dd_table).
It stands around 3 times faster than both `SeqKit` and `fastp` (which hover around 900-1500s). In addition to speed, it uses a quarter of the RAM of `SeqKit` (4.6 GB vs 19 GB) and roughly the same amount as `fastp`.

== : Compressed in, uncompressed out Single-End FASTQ

Here we will evaluate the performances of the tools in a Single-End Compressed in Uncompressed out (SE-cd) scenario.
We got very similar results to those for uncompressed single-end data, with `FastDedup` being the fastest tool for all volumes and using the least amount of RAM up until 5M reads, where `prinseq++` starts to use less memory than `FastDedup` (@SE_cd_bench, @SE_cd_table).

#figure(
    image("results/se_cd_runtime_and_ram.svg", width: 100%),
    caption: [
    Execution time and memory usage for Compressed-input to Uncompressed-output Single-End FASTQ deduplication.
    ],
)<SE_cd_bench>

#pagebreak()
= -- Discussion and future perspectives

Now that we have significant data from isolated benchmarks, we can discuss the implications of these results for real-world applications and the future development of `FastDedup`.

== : Output Quality and Deduplication Accuracy<outqual>

Despite each tool employing very different approaches to deduplication, they all produced nearly identical outputs regarding the number of unique sequences and the duplication rate.
The differences between tools are marginal, staying entirely within a 0.1% margin.
However, we were unable to estimate this metric for `CD-HIT-DUP` because it does not output a log.
Moreover, a critical problem emerges with `prinseq++` starting at 10M reads: its reported duplication rate exceeds the expected 10.1%, climbing to over 90% at 300M reads (@PE_dd_dedup, @pe_rate).
Because of this severe issue, we strongly advise against using `prinseq++` in production pipelines; although it becomes more memory-efficient than `FastDedup` at 5M reads, its false-positive rate becomes unacceptably high by 10M reads.

#figure(
    image("results/pe_dd_duplication_rate.svg", width: 100%),
    caption: [
    Duplication rate observed for Uncompressed Paired-End FASTQ deduplication.
    ],
)<PE_dd_dedup>

#duplication-rate-table(
  "results/pe_dd_duplication.csv", 
  caption: [Paired-end duplication rate analysis.],
  label: <pe_rate>
)

A similar pattern is observed for `clumpify` on single-end data, but starting at 3M reads.
It rises above 10.01% at 10M reads, continuing to climb steadily up to 100M reads.
While this still exceeds the false-positive threshold, it remains significantly better than the 90% false-positive rate exhibited by `prinseq++` at 300M reads (@SE_dd_dedup, @se_rate).

#figure(
    image("results/se_dd_duplication_rate.svg", width: 100%),
    caption: [
    Duplication rate observed for Uncompressed Single-End FASTQ deduplication.
    ],
)<SE_dd_dedup>

#duplication-rate-table(
  "results/se_dd_duplication.csv", 
  caption: [Single-end duplication rate analysis.],
  label: <se_rate>
)

The actual duplication rate is not strictly constant across different sequencing volumes because `wgsim` acts as a random read simulator.
We can calculate the expected number of natural coincidental duplicates $C$ with $N$ randomly generated reads from a uniform distribution of $M$ possible start positions as $C approx N^2 / (2M)$.
However, the human genome is not uniform, as it contains high-frequency repeat families @treangen2012repetitive, so we have to approximate $M$ with the number of unique start positions in the simulated dataset.
We can approximate $k=(2G)^(-1)$ using a known data point from your `wgsim` simulation.
Looking at the generation of $270,000,000$ single-end reads, the simulation naturally produced roughly $60,000$ coincidental duplicates.

We can solve $k$ like:
$ 
60,000 &= k times (270,000,000)^2 \
k &= (60,000) / (270,000,000)^2
$
We can now express the number of natural duplicates $C$ as $C=k times N^2$ and add this value to the artificial 10% duplication rate we intentionally introduced, yielding the theoretical expected duplication rate for each dataset size.

Before leaving this, we have to talk about why `prinseq++`'s false-positive rate explodes so dramatically with increasing dataset size.
We can calculate the expected false-positive rate $F$ for a Bloom filter with $m$ bits and $n$ inserted elements using $k$ hash functions as:
$ F = (1 - e^(-(k n)/m))^k $
We have multiple data points for `prinseq++`'s false-positive rate at different volumes, so we can solve for $m$ using the observed false-positive rates like so:
$ m approx - (k n) / ln(1 - F^(1/k)) $
However, we also need $k$, the number of hash functions used by `prinseq++`, which is not documented.

However I know how to read code. In `main.cpp`, we can find this snippet of code starting at line 497.

```CPP
    if (derep) {
        parameters.projected_element_count = 10000000;
        parameters.false_positive_probability = 0.000001; // 1 in 10000
```

The algorithm begins with two fixed, hardcoded constraints:
- $n$ (Projected elements): $10,000,000$
- $p$ (Target false positive probability): $10^(-6)$

Because a Bloom filter must use a whole number of hash functions, the algorithm performs a discrete optimization.
It iterates over integers $k in ZZ^+$ to find the argument $k$ that minimizes the size $m(k)$:

$ k_"opt" = op("arg min")_(k in ZZ^+) ( (-k dot n) / ln(1 - p^(1/k)) ) $

The continuous global minimum for this function occurs when $p = (1/2)^k$, we can directly solve for the ideal $k$:

$ k approx -ln(p) / ln(2) <=> k approx -ln(10^(-6)) / ln(2) approx 13.8155 / 0.6931 approx 19.931 $

When evaluating the function $m(k)$ across integer values, the discrete minimum strictly aligns with the closest integer to this theoretical continuous minimum. Therefore, the search loop locks in on: $k_"opt" = 20$
With the optimal number of hash functions established, the algorithm substitutes this value back into the objective function to find the raw bit size:

$ m_"raw" = (-20 dot 10,000,000) / ln(1 - (10^(-6))^(1/20))  approx 287,551,586.5 " bits" $

The code then uses `static_cast` to truncate the decimal, resulting in an integer raw size of *287,551,586 bits*.
We saw experimentally that `prinseq++` had a 45MB RAM footprint at 5M reads, which is consistent with this theoretical calculation since $287,551,586$ bits is approximately 34.3 MB.
The rest of the RAM usage can be attributed to the overhead of the Bloom filter data structure.

Now that we have explicitly determined the internal parameters, we can mathematically model how the false-positive rate degrades when the filter is fed a dataset of 100 million paired-end reads.

Crucially, because `prinseq++` was executed using exclusively the `-derep` flag, no prior quality trimming or length filtering occurred.
This guarantees that 100% of the dataset, all 200 million individual sequences, was forced into the Bloom filter.
Assuming an estimated true biological duplication rate of 10%, the approximate number of unique sequences inserted into the filter is $N approx 180,000,000$.

When `prinseq++` outputs the final deduplication percentage, it reports the fraction of reads flagged as duplicates across the entire run.
As the filter fills up, the instantaneous false-positive rate after $x$ unique elements have been inserted, $F(x)$, increases monotonically toward saturation:

$ F(x) = (1 - e^(-(20 x)/m))^20 $

The expected average false-positive rate $overline(F)(N)$ for a dataset containing $N$ unique elements is the integral of this instantaneous rate over the domain of $N$, divided by $N$:

$ overline(F)(1.8 dot 10^8) = 1/(1.8 dot 10^8) integral_0^(1.8 dot 10^8) (1 - e^(-(20 x) / 287551586))^20 upright(d)x $

Evaluating this integral yields the theoretical average false-positive rate among unique reads:

$ overline(F)(1.8 dot 10^8) approx 0.7126 $

This indicates that 71.26% of the unique reads are expected to be incorrectly flagged as duplicates due to filter saturation.
To obtain the total reported duplication rate, we must account for both the true duplicates and the false positives relative to the total number of reads:

$ D_"total" = (D_"true" + overline(F)(N) dot N) / N_"total" = (0.10 times 2 dot 10^8 + 0.7126 times 1.8 dot 10^8) / (2 dot 10^8) approx 0.7414 $

This theoretical prediction of 74.14% is in strikingly close of the empirical observation of 74.04%, differing by only 0.1 percentage points.
The exclusive use of the `-derep` flag is what allows this mathematical model to map so cleanly to reality, as there are no "lost" reads from upstream filters skewing the input volume $N$.
It is worth noting that this model assumes all $N$ unique reads are sequentially inserted into the filter.
In practice, reads that trigger a false positive may not advance the filter state, which would slow saturation.
The close agreement between the theoretical and empirical values suggests this effect is negligible at this scale.
Regardless, the formulation strongly demonstrates that the 34.3 MB hardcoded memory footprint guarantees a catastrophic failure for a 200-million-read insertion volume.

#page(flipped: true)[
    == : In short
    #figure(
        align(center)[
            #table(
                columns: (auto, auto, auto, auto, auto, 1fr),
                align: (col, row) => if row != 0 and col == 5 { left } else { center + horizon },
                stroke: 0.5pt + luma(200),
                fill: (col, row) => if row == 0 { luma(240) } else { none },
                
                [*Tool*], [*Algorithm / Approach*], [*PE \ Support*], [*Error \ Tolerance*], [*Memory \ Usage*], [*Key Limitations*],
                
                [`fastp` @Chen2018 @Chen2023], [Probabilistic (Bloom filter/hashes)], [Yes], [No], [Configurable], [False-positive rate; long runtime],
                
                [`FastUniq` @Xu2012], [Deterministic (Merge sort)], [Yes], [No], [Very High], [Severe RAM bottleneck; sensitive to input file order.],
                
                [`clumpify` @Bushnell2014], [k-mer grouping], [Yes], [Yes], [Variable], [Variable stability on massive datasets; high memory overhead.],
                
                [`seqkit rmdup` @Shen2016 @Shen2024], [Exact matching (Hashing/MD5)], [No], [No], [Very Low], [Desynchronizes paired-end files.],
                
                [`CD-HIT-DUP` @Fu2012 @Huang2010], [Greedy incremental clustering], [Yes], [Yes], [Very High], [Stability issues in containerized environments.],
                
                [`prinseq++` @Cantu2019], [Probabilistic (Bloom filter)], [Yes], [No], [Low], [False-positives rate],
                
                [`FastDedup`], [Exact matching (Hashing/xxh3)], [Yes], [No], [Very Low], [Theoretical negligible collision risk]
            )
        ],
        caption: [Summary of existing tools for reads deduplication and their limitations.],
    ) <tools_table>
    
    #v(1fr)
    #figure(
        align(center)[
            #build-comparaison-table((
            "results/pe_cc.csv", 
            "results/pe_cd.csv", 
            "results/pe_dd.csv", 
            "results/se_cc.csv", 
            "results/se_cd.csv", 
            "results/se_dd.csv"
            ))
        ],
        caption: [Fastest speed of each tool across different dataset sizes and file formats.],
    )
    #figure(
        align(center)[
            #build-ram-comparaison-table((
            "results/pe_cc.csv", 
            "results/pe_cd.csv", 
            "results/pe_dd.csv", 
            "results/se_cc.csv", 
            "results/se_cd.csv", 
            "results/se_dd.csv"
            ))
        ],
        caption: [Lowest RAM of each tool across different dataset sizes and file formats. Prinseq++ is excluded from this table because of the duplication rate issue.],
    )
]
= -- Conclusion and Recommendations

`FastDedup` represents a significant improvement for FASTX deduplication, combining unmatched speed with a remarkably low memory footprint across all tested scenarios.
Its design, shifting the deduplication problem from string-matching to integer-matching through non-cryptographic hashing, allows it to scale from small experimental datasets to the largest population-scale whole-genome sequencing projects without changing configuration or requiring HPC-class hardware.

== : Summary of findings

For paired-end data, `FastDedup` is the clear winner on both dimensions.
It consistently outperforms all competitors in execution time, being up to 6× faster than `prinseq++` and more than 10× faster than `fastp` on uncompressed 300M-read datasets, while keeping its memory footprint well below that of every other tool that survived to large volumes except for `fastp` at 300M reads.
`FastUniq` and `CD-HIT-DUP` are eliminated by memory constraints before reaching 100M reads, and `clumpify` crashes at 100M reads under the 32GB RAM budget regardless of format for paired-end cc and cd.
`prinseq++`, despite being the most memory-efficient competitor, disqualifies itself on correctness grounds: its Bloom filter false-positive rate explodes beyond 10% starting at 10M reads, reaching above 90% at 300M reads, making it unsuitable for any analysis where duplicate removal must be reliable.
For single-end data, the picture is more nuanced.

`SeqKit` outperforms `FastDedup` on speed when both input and output are compressed, owing to its mature multi-threaded compression pipeline.
However, `FastDedup` reverses this advantage completely on uncompressed I/O, where it runs around 3 times faster than both `SeqKit` and `fastp` at 300M reads, and maintains a drastically lower memory footprint throughout: at 300M reads, `SeqKit` consumes approximately 19 GB of RAM compared to `FastDedup`'s 4.6 GB, a fourfold difference.
Users working on memory-constrained nodes or seeking maximum throughput through uncompressed intermediate files should therefore prefer `FastDedup` even for single-end libraries.

== : Practical recommendations

For maximum performance in any pipeline, we recommend the following strategy: run `FastDedup` immediately after quality control, and direct its output to an uncompressed temporary file.
Writing uncompressed output eliminates the compression bottleneck that otherwise dominates execution time, and the temporary files can be recompressed or deleted once downstream steps are complete.
For users whose storage constraints require compressed outputs throughout, the compressed-output mode of `FastDedup` remains faster than all alternatives for paired-end data, and is still competitive for single-end data except for compressed outputs where `SeqKit` leads on speed alone.

The 64-bit hashing mode is appropriate for the overwhelming majority of datasets: the collision probability remains below 1‰ for libraries of up to approximately 600 million paired-end reads (1.2 billion total sequences).
For larger datasets, the automatic upgrade to 128-bit hashing provides the same pipeline interface with negligible runtime overhead and approximately double the memory cost, a trade-off that remains highly favourable compared to any deterministic alternative.

We advise against `prinseq++` in any production pipeline handling more than 5M reads.
Their false-positive rates are not a theoretical concern but an empirically observed failure: they remove correct unique sequences at increasing rates as dataset size grows, silently corrupting downstream variant calling and assembly results.
`CD-HIT-DUP`, `FastUniq` and `clumpify` remain valid choices only for small, uncompressed datasets where exact error tolerance or mismatch tolerance are required and where memory is not a constraint.\
We have to warn about `fastp`, since it use the same bloom filter approach as `prinseq++`, it is likely to suffer from the same false-positive rate explosion at larger volumes, though we did not observe this within the tested range.
== : Limitations and future work

The current version of `FastDedup` is strictly single-threaded.
While this simplifies the zero-allocation pipeline and guarantees deterministic output ordering, it leaves performance on the table for users with access to multiple cores.
The most impactful near-term improvement would be parallelizing the compression and decompression stages, which account for more than 90% of wall-clock time in compressed-output modes.
Parallelizing the hashing and hash-set lookup stages would offer additional gains, though this requires careful handling of concurrent access to the shared hash set.

A second avenue for development is abundance tracking: rather than discarding duplicate records entirely, `FastDedup` could optionally record the copy number of each unique sequence in the output headers.
This would extend its utility to metagenomics and transcriptomics workflows where abundance information carries biological meaning.

We can also note that the current benchmarking suite covers synthetic data derived from a single human reference genome.
Extending the evaluation to metagenomic, amplicon, and low-complexity datasets would provide a more complete picture of performance across the full range of sequencing applications.
It should be noted that all these benchmarks were performed with `FastDedup` version 1.0.0; improvements and optimizations have been made since then, so the results documented in this paper may not fully reflect the current, likely enhanced performance of `FastDedup`.

Finally, `FastDedup` could continue to expand its feature set to include additional pre-processing steps such as quality filtering, adapter and sequence trimming for example.

`FastDedup` is available under the MIT License at:
- GitHub #link("https://github.com/RaphaelRibes/FastDedup")
- Bioconda #link("https://anaconda.org/bioconda/fdedup")
- Cargo #link("https://crates.io/crates/fastdedup")

The code of the benchmark can be found at: #link("https://github.com/RaphaelRibes/bench-fdedup.git")

= -- Acknowledgements
#v(1em)

This paper was written by Raphaël Ribes, reviewed and edited by Céline Mandier, with the assistance of Gemini 3.1 for grammar, and language correction.
Code architecture and design was developed by Raphaël Ribes, with autocompletion assisted by Gemini 3.1.
Feedback on the algorithm design and ideas for features were provided by Céline Mandier and Alban Mancheron, the implementation was done by Raphaël Ribes.
The testsuite was mostly developed by Claude Sonnet 4.6, reviewed by Raphaël Ribes.

Raphaël Ribes was supported through an internship grant awarded to Alice Baniel by the University of Montpellier (I-SITE Excellence Programme 2022 - SOFIT).

Computations were performed on the ISDM-MESO HPC platform, funded in the framework of State-region planning contracts (Contrat de plan État-région – CPER) by the French Government, the Occitanie/Pyrénées-Méditerranée Region, Montpellier Méditerranée Métropole, and the University of Montpellier.

We would like to thank the Data Science Institute of Montpellier (ISDM) funded by the ANR and France 2030 under the UM2030 programme (ANR-21-IDES-0005), for their expertise and their valuable support.
#pagebreak()
#bibliography("works.bib", style: "american-physics-society")


#page(flipped: true)[
    = -- Annexes
    Note -- In all Annex tables, values highlighted in #text(fill: green)[green] denote tools whose performance is not significantly worse than the best tool in that column.
    The comparison uses a one-sided Welch's _t_-test (H#sub[0]: the tool is no slower than the column-best) with a family-wise error rate of #sym.alpha = 0.05 across all pairwise comparisons within each table.
    `prinseq++` is excluded from this comparison due to the duplication accuracy issues discussed in @outqual
    Its values are never highlighted regardless of performance. A tool appearing without green highlighting is significantly slower or more memory-intensive than the column-best at the corrected significance level.
    #bench-table(
        "results/pe_cc.csv",
        label: <PE_cc_table>,
        caption: [
            Execution time and memory usage for compressed input and compressed output paired-end FASTQ files.
        ],
    )
    #bench-table(
        "results/pe_dd.csv",
        label: <PE_dd_table>,
        caption: [
            Execution time and memory usage for uncompressed input and uncompressed output paired-end FASTQ files.
        ],
    )

    #bench-table(
        "results/pe_cd.csv",
        label: <PE_cd_table>,
        caption: [
            Execution time and memory usage for compressed input and uncompressed output paired-end FASTQ files.
        ],
    )
    
    #bench-table(
        "results/se_cc.csv",
        label: <SE_cc_table>,
        caption: [
            Execution time and memory usage for compressed input and compressed output single-end FASTQ files.
        ],
    )

    #bench-table(
        "results/se_dd.csv",
        label: <SE_dd_table>,
        caption: [
            Execution time and memory usage for uncompressed input and uncompressed output single-end FASTQ files.
],
    )
    
    #bench-table(
        "results/se_cd.csv",
        label: <SE_cd_table>,
        caption: [
            Execution time and memory usage for compressed input and uncompressed output single-end FASTQ files.
        ],
    )
]