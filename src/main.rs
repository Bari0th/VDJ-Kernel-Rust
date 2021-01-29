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
        let mut result = vec![0; (sequence.len() + 1) * (gene.len() + 1)];
        let size = gene.len() + 1;
        for (j, res) in result.iter_mut().enumerate().take(size) {
            *res = -2 * j as i32;
        }

        for i in 1..=sequence.len() {
            for j in 1..gene.len() {
                let skip = result[(i - 1) * size + j] - 3;
                let skip_seq = result[i * size + j - 1] - 3;
                let step = result[(i - 1) * size + j - 1]
                    + if sequence[sequence.len() - i] == gene[gene.len() - j] {
                        2
                    } else {
                        -1
                    };
                result[i * size + j] = skip.max(skip_seq.max(step));
            }
            let j = gene.len();

            let skip = result[(i - 1) * size + j];
            let skip_seq = result[i * size + j - 1] - 3;
            let step = result[(i - 1) * size + j - 1]
                + if sequence[sequence.len() - i] == gene[gene.len() - j] {
                    2
                } else {
                    -1
                };
            result[i * size + j] = skip.max(skip_seq.max(step));
        }

        result[result.len() - 1]
    }

    pub fn get_best_alignment(gene: &[Nucleotide], sequence: &[Nucleotide]) -> (i32, String) {
        let mut result = vec![(0, 2); (sequence.len() + 1) * (gene.len() + 1)];
        let size = gene.len() + 1;
        for (j, res) in result.iter_mut().enumerate().take(size) {
            *res = (-2 * j as i32, 0);
        }

        for i in 1..=sequence.len() {
            for j in 1..gene.len() {
                let skip = (result[(i - 1) * size + j].0 - 3, 2);
                let skip_seq = (result[i * size + j - 1].0 - 3, 0);
                let step = (
                    result[(i - 1) * size + j - 1].0
                        + if sequence[sequence.len() - i] == gene[gene.len() - j] {
                            2
                        } else {
                            -4
                        },
                    1,
                );
                result[i * size + j] = skip.max(skip_seq.max(step));
            }
            let j = gene.len();

            let skip = (result[(i - 1) * size + j].0, 2);
            let skip_seq = (result[i * size + j - 1].0 - 3, 0);
            let step = (
                result[(i - 1) * size + j - 1].0
                    + if sequence[sequence.len() - i] == gene[gene.len() - j] {
                        2
                    } else {
                        -4
                    },
                1,
            );
            result[i * size + j] = skip.max(skip_seq.max(step));
        }

        let mut i = sequence.len();
        let mut j = gene.len();
        let mut s = String::new();

        while i != 0 || j != 0 {
            match result[i * size + j].1 {
                2 => {
                    s.push('.');
                    i -= 1;
                }
                1 => {
                    s.push(gene[gene.len() - j].into());
                    j -= 1;
                    i -= 1;
                }
                0 => {
                    j -= 1;
                }
                _ => unreachable!("Should not be reached !"),
            }
        }
        (result[result.len() - 1].0, s)
    }
}

use std::ops::Deref;

impl Deref for Sequence {
    type Target = [Nucleotide];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Into<String> for &Sequence {
    fn into(self) -> String {
        self.0.iter().map(|&n| Into::<char>::into(n)).collect()
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&Into::<String>::into(self))
    }
}

#[derive(Debug)]
struct Genes(Vec<(String, Sequence)>);

impl Genes {
    pub fn get_best<'a>(&'a self, sequence: &[Nucleotide]) -> &'a str {
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

    pub fn get_gene<'a>(&'a self, name: &str) -> Option<&'a Sequence> {
        self.0
            .iter()
            .filter_map(|(s, seq)| if s == name { Some(seq) } else { None })
            .next()
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

    pub fn calcul_and_save_to(
        path: &str,
        sequences: &Sequences,
        genes_v: &Genes,
        genes_d: &Genes,
        genes_j: &Genes,
        results: &SequenceResult,
    ) -> Self {
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create(path).expect("Cannot create file");
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
                    println!("{}", sequence.0);
                    let result = results.0.get(&tmp.0.to_string()).unwrap();
                    file.write_all(
                        (format!(
                            "{} :\n{}\n\n{}\n{}\n{}\n\n",
                            tmp.0,
                            sequence.1.to_string(),
                            Sequence::get_best_alignment(
                                genes_v.get_gene(&tmp.1 .0).unwrap(),
                                &sequence.1
                            )
                            .1,
                            Sequence::get_best_alignment(
                                genes_d.get_gene(&tmp.1 .1).unwrap(),
                                &sequence.1
                            )
                            .1,
                            Sequence::get_best_alignment(
                                genes_j.get_gene(&tmp.1 .2).unwrap(),
                                &sequence.1
                            )
                            .1,
                        ) + &format!(
                            "{}\n{}\n{}\n\n",
                            if let Some(gene) = genes_v.get_gene(&result.0) {
                                Sequence::get_best_alignment(gene, &sequence.1).1
                            } else {
                                ('.'..='.')
                                    .cycle()
                                    .take(sequence.1.len())
                                    .collect::<String>()
                            },
                            if let Some(gene) = genes_d.get_gene(&result.1) {
                                Sequence::get_best_alignment(gene, &sequence.1).1
                            } else {
                                ('.'..='.')
                                    .cycle()
                                    .take(sequence.1.len())
                                    .collect::<String>()
                            },
                            if let Some(gene) = genes_j.get_gene(&result.2) {
                                Sequence::get_best_alignment(gene, &sequence.1).1
                            } else {
                                ('.'..='.')
                                    .cycle()
                                    .take(sequence.1.len())
                                    .collect::<String>()
                            },
                        ))
                            .as_bytes(),
                    )
                    .expect("Cannot write to file");
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
        let (size, hash_index) =
            genes
                .0
                .iter()
                .fold((0, HashMap::new()), |(idx, mut hash_map), gene| {
                    if hash_map.contains_key(&gene.0) {
                        (idx, hash_map)
                    } else {
                        hash_map.insert(gene.0.clone(), idx);
                        (idx + 1, hash_map)
                    }
                });
        let confusion = vec![vec![0; size]; size];
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
    let genes_j = Genes::load_from("data/database/IGHJgenes.txt");
    let genes_d = Genes::load_from("data/database/IGHDgenes.txt");
    let genes_v = Genes::load_from("data/database/IGHVgenes.txt");

    let sequences = Sequences::load_from("data/datasets/simulations/sim1.fasta");

    let results = SequenceResult::load_from("data/datasets/simulations/simTrueGenes.txt");

    let predictions = SequenceResult::calcul_and_save_to(
        "test_alignment/genes_alignment.txt",
        &sequences,
        &genes_v,
        &genes_d,
        &genes_j,
        &results,
    );
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
