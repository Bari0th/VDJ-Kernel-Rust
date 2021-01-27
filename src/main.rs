#![allow(dead_code)]

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
    Unknown,
}

use std::convert::From;

impl From<char> for Nucleotide {
    fn from(ch: char) -> Nucleotide {
        use Nucleotide::*;
        match ch {
            'A' | 'a' => A,
            'C' | 'c' => C,
            'G' | 'g' => G,
            'T' | 't' => T,
            'n' => Unknown,
            _ => panic!("Invalid character : {:?}", ch),
        }
    }
}

impl Into<char> for Nucleotide {
    fn into(self) -> char {
        use Nucleotide::*;
        match self {
            A => 'a',
            C => 'c',
            G => 'g',
            T => 't',
            Unknown => 'n',
        }
    }
}

#[derive(Debug, Clone)]
pub struct Sequence(Vec<Nucleotide>);

impl Sequence {
    pub fn new(s: &str) -> Sequence {
        Sequence(
            s.chars()
                .filter(|&c| c != '\r' && c != '\n' && c != ' ')
                .filter_map(|ch| {
                    match ch {
                        '.' => None, //ignore .
                        _ => Some(Nucleotide::from(ch)),
                    }
                })
                .collect(),
        )
    }

    pub fn compare(gene: &[Nucleotide], sequence: &[Nucleotide]) -> i32 {
        if gene.is_empty() {
            -(sequence.len() as i32)
        } else if sequence.is_empty() {
            -(gene.len() as i32)
        } else {
            let mut result = vec![None; sequence.len() * gene.len()];

            Sequence::compare_rec_first(&gene, &sequence, &mut result)
        }
    }

    fn compare_rec_first(
        gene: &[Nucleotide],
        sequence: &[Nucleotide],
        result: &mut [Option<i32>],
    ) -> i32 {
        if sequence.is_empty() {
            -(gene.len() as i32) * 3
        } else {
            let skip_seq = Sequence::compare_rec_first(gene, &sequence[1..], result);

            let step = Sequence::compare_rec(&gene[1..], &sequence[1..], gene.len(), result)
                + if gene[0] == sequence[0] { 1 } else { -1 };

            let skip_gene = Sequence::compare_rec(&gene[1..], sequence, gene.len(), result) - 3;

            skip_seq.max(step).max(skip_gene)
        }
    }

    /*
        CGT  .CGT
        ACGT

        AGCTACGT..... | ACGTA
        AGCG..CGCAGC. | ACTA

        CGT  CGT
        GTAC .GTAC

        AGTC  AG....TC
        .....AGTC..
        AGAAAA.TCAA
    */

    fn compare_rec(
        gene: &[Nucleotide],
        sequence: &[Nucleotide],
        gene_size: usize,
        result: &mut [Option<i32>],
    ) -> i32 {
        if gene.is_empty() {
            0
        } else if sequence.is_empty() {
            -(gene.len() as i32) * 3
        } else {
            let idx = gene_size * (sequence.len() - 1) + gene.len() - 1;
            match result[idx] {
                None => {
                    if gene[0] == sequence[0] {
                        let step =
                            Sequence::compare_rec(&gene[1..], &sequence[1..], gene_size, result)
                                + 1;
                        result[idx] = Some(step);
                        step
                    } else {
                        let skip_seq =
                            Sequence::compare_rec(gene, &sequence[1..], gene_size, result) - 2;

                        let step =
                            Sequence::compare_rec(&gene[1..], &sequence[1..], gene_size, result)
                                - 1;

                        let skip_gene =
                            Sequence::compare_rec(&gene[1..], sequence, gene_size, result) - 2;

                        let max = skip_seq.max(step).max(skip_gene);
                        result[idx] = Some(max);
                        max
                    }
                }
                Some(val) => val,
            }
        }
    }
}

use std::ops::Deref;

impl Deref for Sequence {
    type Target = [Nucleotide];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

#[derive(Debug)]
struct Genes(Vec<(String, Sequence)>);

impl Genes {
    pub fn get_best<'a>(&'a self, sequence: &Sequence) -> &'a str {
        self.0
            .iter()
            .fold(("No genes", i32::MIN), |res, gene| {
                let score = Sequence::compare(&gene.1, sequence);
                if score > res.1 {
                    (&gene.0, score)
                } else {
                    res
                }
            })
            .0
    }

    pub fn load_from(path: &str) -> Self {
        use std::fs::File;
        use std::io::Read;
        let mut file = File::open(path).expect("Cannot open file");

        let mut buffer = String::new();

        file.read_to_string(&mut buffer).expect("Cannot read file");

        Genes(
            buffer
                .split('>')
                .skip(1)
                .map(|s| {
                    let mut format_gene = s.split('|');
                    (
                        format_gene.nth(1).unwrap().to_string(),
                        Sequence::new(format_gene.last().unwrap()),
                    )
                })
                .collect(),
        )
    }
}

#[derive(Debug)]
struct Sequences(Vec<(String, Sequence)>);

impl Sequences {
    pub fn load_from(path: &str) -> Self {
        use std::fs::File;
        use std::io::Read;
        let mut file = File::open(path).expect("Cannot open file");

        let mut buffer = String::new();

        file.read_to_string(&mut buffer).expect("Cannot read file");

        Sequences(
            buffer
                .split('>')
                .skip(1)
                .map(|s| {
                    let mut format_sequence = s.split('\n');
                    (
                        format_sequence.next().unwrap().trim().to_string(),
                        Sequence::new(format_sequence.next().unwrap()),
                    )
                })
                .collect(),
        )
    }
}

use std::collections::HashMap;

#[derive(Debug)]
struct SequenceResult(HashMap<String, (String, String, String)>);

impl SequenceResult {
    pub fn load_from(path: &str) -> Self {
        use std::fs::File;
        use std::io::Read;
        let mut file = File::open(path).expect("Cannot open file");

        let mut buffer = String::new();

        file.read_to_string(&mut buffer).expect("Cannot read file");

        SequenceResult(
            buffer
                .split('\n')
                .rev()
                .skip(1)
                .map(|s| {
                    let mut format_result = s.trim().split('\t');
                    (
                        format_result
                            .next()
                            .expect("Problem reading formated date !")
                            .to_string(),
                        (
                            format_result
                                .next()
                                .expect("Problem reading formated date !")
                                .to_string(),
                            format_result
                                .next()
                                .expect("Problem reading formated date !")
                                .to_string(),
                            format_result
                                .next()
                                .expect("Problem reading formated date !")
                                .to_string(),
                        ),
                    )
                })
                .collect(),
        )
    }

    pub fn save_to(&self, path: &str) {
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create(path).expect("Cannot create file");
        file.write_all(
            self.0
                .iter()
                .map(|(key, value)| format!("{}\t{}\t{}\t{}\n", key, value.0, value.1, value.2))
                .collect::<String>()
                .as_bytes(),
        )
        .expect("Cannot write to file");
    }

    pub fn calcul_from(
        sequences: &Sequences,
        genes_v: &Genes,
        genes_d: &Genes,
        genes_j: &Genes,
    ) -> Self {
        SequenceResult(
            sequences
                .0
                .iter()
                .map(|sequence| {
                    (
                        sequence.0.clone(),
                        (
                            genes_v.get_best(&sequence.1).to_string(),
                            genes_d.get_best(&sequence.1).to_string(),
                            genes_j.get_best(&sequence.1).to_string(),
                        ),
                    )
                })
                .collect(),
        )
    }

    pub fn calcul_from_and_compare(
        sequences: &Sequences,
        genes_v: &Genes,
        genes_d: &Genes,
        genes_j: &Genes,
        resultats: &SequenceResult,
    ) -> Self {
        SequenceResult(
            sequences
                .0
                .iter()
                .map(|sequence| {
                    let tmp = (
                        sequence.0.clone(),
                        (
                            genes_v.get_best(&sequence.1).to_string(),
                            genes_d.get_best(&sequence.1).to_string(),
                            genes_j.get_best(&sequence.1).to_string(),
                        ),
                    );
                    let resultat = resultats.0.get(&tmp.0.to_string()).unwrap();
                    println!(
                        "{} :\t{} {} {}\n\t{} {} {}\n\n",
                        tmp.0, resultat.0, resultat.1, resultat.2, tmp.1 .0, tmp.1 .1, tmp.1 .2
                    );
                    tmp
                })
                .collect(),
        )
    }

    pub fn get_confusion_matrixes(
        predictions: &SequenceResult,
        results: &SequenceResult,
        genes_v: &Genes,
        genes_d: &Genes,
        genes_j: &Genes,
    ) -> (ConfusionMatrix, ConfusionMatrix, ConfusionMatrix) {
        let mut confusion_v = ConfusionMatrix::new(genes_v);
        let mut confusion_d = ConfusionMatrix::new(genes_d);
        let mut confusion_j = ConfusionMatrix::new(genes_j);
        predictions.0.iter().for_each(|(name, prediction)| {
            let result = results
                .0
                .get(name)
                .expect("Cannot get result for this sequence");
            confusion_v.increment(&prediction.0, &result.0);
            confusion_d.increment(&prediction.1, &result.1);
            confusion_j.increment(&prediction.2, &result.2);
        });
        (confusion_v, confusion_d, confusion_j)
    }
}

struct ConfusionMatrix(Vec<Vec<usize>>, HashMap<String, usize>);

impl ConfusionMatrix {
    pub fn new(genes: &Genes) -> ConfusionMatrix {
        let confusion = vec![vec![0; genes.0.len()]; genes.0.len()];
        let hash_index = genes
            .0
            .iter()
            .enumerate()
            .map(|(idx, gene)| (gene.0.clone(), idx))
            .collect();
        ConfusionMatrix(confusion, hash_index)
    }

    pub fn increment(&mut self, prediction: &str, result: &str) {
        if let Some(&idx) = self.1.get(prediction) {
            let i = idx;
            if let Some(&idx) = self.1.get(result) {
                let j = idx;
                self.0[i][j] += 1;
            }
        }
    }
}

use std::fmt;

impl fmt::Display for ConfusionMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            self.0
                .iter()
                .map(|line| format!("{:?}\n", line))
                .collect::<String>()
        )
    }
}

fn main() {
    let genes_j = Genes::load_from("../data/database/IGHJgenes.txt");
    let genes_d = Genes::load_from("../data/database/IGHDgenes.txt");
    let genes_v = Genes::load_from("../data/database/IGHVgenes.txt");

    let sequences = Sequences::load_from("../data/datasets/simulations/sim1.fasta");

    let results = SequenceResult::load_from("../data/datasets/simulations/simTrueGenes.txt");

    let predictions =
        SequenceResult::calcul_from_and_compare(&sequences, &genes_v, &genes_d, &genes_j, &results);
    // let predictions = SequenceResult::load_from("../data/datasets/simulations/simTrueGenes.txt");

    predictions.save_to("first_predictions");

    let confusion_matrixes = SequenceResult::get_confusion_matrixes(
        &predictions,
        &results,
        &genes_v,
        &genes_d,
        &genes_j,
    );

    println!(
        "{}\n\n{}\n\n{}",
        confusion_matrixes.0, confusion_matrixes.1, confusion_matrixes.2
    );
}
