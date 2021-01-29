use super::prelude::*;
use std::fmt;

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
