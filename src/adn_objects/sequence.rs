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

    pub fn compare<const N: usize>(gene: &[Nucleotide], sequence: &[Nucleotide]) -> i32 {
        const GENE_DELETION: i32 = 3;
        const GENE_INSERTION: i32 = 3;
        const GENE_MUTATION: i32 = 1;
        const MATCHING_GENE: i32 = 2;
        let mut result = [0; N];
        let offset = N / 2;
        let size = gene.len().min(sequence.len());
        result[0] = (gene.len() - size) as i32 * (-GENE_DELETION);

        for idx_gen in 1..N - offset {
            result[idx_gen] = result[idx_gen - 1] - GENE_INSERTION;
        }

        for idx_seq in 0..=offset {
            result[offset - idx_seq] = result[offset - idx_seq + 1] - GENE_DELETION;
            for idx_gen in offset - idx_seq + 1..N - 1 {
                let insertion = result[idx_gen + 1] - GENE_INSERTION;
                let deletion = result[idx_gen - 1] - GENE_DELETION;
                let matching = result[idx_gen]
                    + if gene[idx_seq + idx_gen - offset] == sequence[idx_seq] {
                        MATCHING_GENE
                    } else {
                        -GENE_MUTATION
                    };
                result[idx_gen] = insertion.max(deletion.max(matching));
            }

            let deletion = result[N - 2] - GENE_DELETION;
            let matching = result[N - 1]
                + if gene[idx_seq + N - 1 - offset] == sequence[idx_seq] {
                    MATCHING_GENE
                } else {
                    -GENE_MUTATION
                };
            result[N - 1] = deletion.max(matching);
        }

        for idx_seq in offset + 1..=size + offset - N {
            let insertion = result[1] - GENE_INSERTION;
            let matching = result[0]
                + if gene[idx_seq - offset] == sequence[idx_seq] {
                    MATCHING_GENE
                } else {
                    -GENE_MUTATION
                };
            result[0] = insertion.max(matching);
            for idx_gen in 1..N - 1 {
                let insertion = result[idx_gen + 1] - GENE_INSERTION;
                let deletion = result[idx_gen - 1] - GENE_DELETION;
                let matching = result[idx_gen]
                    + if gene[idx_seq + idx_gen - offset] == sequence[idx_seq] {
                        MATCHING_GENE
                    } else {
                        -GENE_MUTATION
                    };
                result[idx_gen] = insertion.max(deletion.max(matching));
            }

            let deletion = result[N - 2] - GENE_DELETION;
            let matching = result[N - 1]
                + if gene[idx_seq + N - 1 - offset] == sequence[idx_seq] {
                    MATCHING_GENE
                } else {
                    -GENE_MUTATION
                };
            result[N - 1] = deletion.max(matching);
        }

        for idx_seq in size + offset - N + 1..size {
            let insertion = result[1] - GENE_INSERTION;
            let matching = result[0]
                + if gene[idx_seq - offset] == sequence[idx_seq] {
                    MATCHING_GENE
                } else {
                    -GENE_MUTATION
                };
            result[0] = insertion.max(matching);
            for idx_gen in 1..size - idx_seq {
                let insertion = result[idx_gen + 1] - GENE_INSERTION;
                let deletion = result[idx_gen - 1] - GENE_DELETION;
                let matching = result[idx_gen]
                    + if gene[idx_seq + idx_gen - offset] == sequence[idx_seq] {
                        MATCHING_GENE
                    } else {
                        -GENE_MUTATION
                    };
                result[idx_gen] = insertion.max(deletion.max(matching));
            }
        }

        result[0]
    }

    pub fn compare_all_sequence(gene: &[Nucleotide], sequence: &[Nucleotide]) -> i32 {
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
                    s.push(gene[gene.len() - j] as isize as u8 as char);
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
        self.0.iter().map(|&n| n as isize as u8 as char).collect()
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&Into::<String>::into(self))
    }
}
