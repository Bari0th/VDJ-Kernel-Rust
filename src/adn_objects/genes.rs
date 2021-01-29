use super::prelude::*;
use std::ops::Deref;

#[derive(Debug)]
pub struct Genes(Vec<(String, Sequence)>);

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
        self.iter()
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

impl Deref for Genes {
    type Target = Vec<(String, Sequence)>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
