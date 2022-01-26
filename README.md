# OCTOPUS
[![pytest-integration](https://github.com/octantbio/octopus/actions/workflows/pytest-integration.yml/badge.svg)](https://github.com/octantbio/octopus/actions/workflows/pytest-integration.yml)

OCTOPUS (Optimized Cloning Through Open-source Pipelined Unabridged Sequencing), is a light-weight, cost-effective, and robust method for full-plasmid sequence verification using next-generation sequencing. Importantly, OCTOPUS can be performed at scale with common lab equipment and uses crude E. coli lysate (rather than purified DNA) as the input. This reduces up-front capital investment by eliminating the need for expensive nano-volume liquid handlers and automated plasmid purification pipelines. Using OCTOPUS, a single molecular biology researcher can sequence verify 96-1152 colonies in three days for ~$5 a sample.

## How OCTOPUS Works

![OCTOPUS overview](./img/overview.jpg)

OCTOPUS is divided into six main steps.
1. colony picking and overnight growth
2. lysate and glycerol stock generation
3. iGenomX Riptide kit library preparation
4. Illumina sequencing
5. turn-key data analysis
6. pick perfect sequences

## Why we Love It

We've used OCTOPUS to sequence ~10,000 plasmids in a 6 month period.

![Plasmids over time](./img/wells-over-time.png)


## Getting Started

- [Installation](https://github.com/octantbio/octopus/wiki/Installation)
- [Bench Protocol](https://github.com/octantbio/octopus/wiki/Bench-Protocol)
- [Running the Analysis Pipeline](https://github.com/octantbio/octopus/wiki/Running-the-Analysis-Pipeline)

## Diving Deeper

- [Pipeline Details](https://github.com/octantbio/octopus/wiki/Pipeline-Details)
- [Experimental Validation](https://github.com/octantbio/octopus/wiki/Experimental-Validation)

## Contributing

Please feel free to open an issue or pull request.

## License

This project is licensed under the Apache 2.0 License - see the [LICENSE](LICENSE) file for details. Additional licensing information:

- [BBTools](docker/bbtools-license) - custom
- [mlr](docker/mlr-license) - custom
- [starcode](docker/starcode-license) - GPL-3.0
- [SPAdes](docker/spades-license) - GPL-2.0
- [fgbio](https://github.com/fulcrumgenomics/fgbio/blob/master/LICENSE) - MIT
- [minimap2](https://github.com/lh3/minimap2/blob/master/LICENSE.txt) - MIT
- [samtools](https://github.com/samtools/samtools/blob/develop/LICENSE) - MIT
- [bcftools](https://github.com/samtools/bcftools/blob/develop/LICENSE) - MIT
- [htslib](https://github.com/samtools/htslib/blob/develop/LICENSE) - MIT
- [freebayes](https://github.com/ekg/freebayes/blob/master/LICENSE) - MIT

